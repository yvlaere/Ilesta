use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Default)]
pub struct RescueStats {
    pub input_read_count: usize,
    pub mapped_read_count: usize,
    pub rescued_read_count: usize,
    pub rescued_read_bases: u64,
    pub rescued_unitig_count: usize,
    pub rescued_assembly_length: u64,
    pub rescued_node_ids: HashSet<String>,
    pub rescue_primary_reads: Option<PathBuf>,
}

pub fn run_rescue_stage(
    filtered_reads: &Path,
    improved_assembly: &Path,
    existing_mapping_paf: Option<&Path>,
    output_dir: &Path,
    config: &crate::configs::AssembleConfig,
    graph: &mut crate::create_overlap_graph::OverlapGraph,
    overlaps: &mut std::collections::HashMap<(usize, usize), crate::alignment_filtering::Overlap>,
) -> std::io::Result<RescueStats> {
    let generated_paf = output_dir.join("rescue_mapping.paf");
    let mapping_paf = if let Some(paf) = existing_mapping_paf {
        paf
    } else {
        map_reads(
            filtered_reads,
            improved_assembly,
            &generated_paf,
            config.threads,
            config.read_type,
        )?;
        &generated_paf
    };
    let mapped_reads = mapped_read_names(mapping_paf)?;
    let rescue_path = output_dir.join("rescue_reads.fq");
    let mut stats = write_rescue_reads(filtered_reads, &rescue_path, &mapped_reads)?;
    println!("Reads considered for rescue: {}", stats.input_read_count);
    println!(
        "Reads mapped to primary assembly (excluded): {}",
        stats.mapped_read_count
    );
    println!(
        "Reads written to rescue_reads.fq: {}",
        stats.rescued_read_count
    );
    println!("Rescue read bases: {}", stats.rescued_read_bases);

    if stats.rescued_read_count == 0 {
        println!("No poorly represented reads found; skipping rescue assembly.");
        return Ok(stats);
    }

    let assembly_dir = output_dir.join("rescue_graph");
    std::fs::create_dir_all(&assembly_dir)?;
    let (_, improved_assembly_length) = fasta_stats(improved_assembly)?;
    let genome_size = improved_assembly_length.clamp(1, u32::MAX as u64) as u32;
    let rescue_target_coverage =
        rescue_target_coverage(stats.rescued_read_bases, improved_assembly_length);
    let rescue_bases_for_assembly = improved_assembly_length
        .saturating_mul(u64::from(rescue_target_coverage))
        .min(stats.rescued_read_bases);
    println!(
        "Rescue assembly target: {}× primary assembly length ({} bases) to include all {} unmapped bases",
        rescue_target_coverage, rescue_bases_for_assembly, stats.rescued_read_bases
    );
    let rescue_paf = assembly_dir.join("rescue.paf");
    let rescue_reads = crate::align_reads::align_reads(
        &rescue_path,
        config.threads.clamp(1, 4),
        &rescue_paf,
        &assembly_dir,
        config.min_read_length.saturating_sub(1).clamp(1, 500),
        config.min_base_quality,
        Some(genome_size),
        rescue_target_coverage,
        &config.minimap_batch_size,
        Some("0.001"),
        Some("10,5000"),
        config.read_type,
    )?;
    stats.rescue_primary_reads = Some(rescue_reads.primary_reads);
    let rescue_overlaps = crate::alignment_filtering::run_alignment_filtering(
        &rescue_paf,
        &config.min_overlap_length,
        &config.min_overlap_count,
        &config.min_percent_identity,
        &config.overhang_ratio,
    )
    .map_err(|error| std::io::Error::new(std::io::ErrorKind::Other, error.to_string()))?
    .overlaps;
    let rescue_graph = crate::create_overlap_graph::run_create_overlap_graph(&rescue_overlaps)?;
    stats.rescued_node_ids = rescue_graph.nodes.keys().cloned().collect();

    let key_offset = overlaps
        .keys()
        .flat_map(|(query, target)| [*query, *target])
        .max()
        .unwrap_or(0)
        .saturating_add(1);
    for ((query, target), overlap) in rescue_overlaps {
        overlaps.insert(
            (
                query.saturating_add(key_offset),
                target.saturating_add(key_offset),
            ),
            overlap,
        );
    }

    for node_id in rescue_graph.nodes.keys() {
        graph.add_node(node_id.clone());
    }
    for (node_id, node) in rescue_graph.nodes {
        for edge in node.edges {
            graph.add_edge(
                &node_id,
                &edge.target_id,
                edge.edge_len,
                edge.overlap_len,
                edge.identity,
            );
        }
    }
    Ok(stats)
}

fn rescue_target_coverage(rescued_bases: u64, assembly_length: u64) -> u32 {
    if assembly_length == 0 {
        return u32::MAX;
    }
    rescued_bases
        .div_ceil(assembly_length)
        .max(1)
        .min(u32::MAX as u64) as u32
}

fn map_reads(
    reads: &Path,
    assembly: &Path,
    paf: &Path,
    threads: usize,
    read_type: crate::cli::ReadType,
) -> std::io::Result<()> {
    let mut command = Command::new("minimap2");
    for arg in read_type.mapping_args() {
        command.arg(arg);
    }
    let status = command
        .arg("-t")
        .arg(threads.to_string())
        .arg(assembly)
        .arg(reads)
        .arg("-o")
        .arg(paf)
        .status()?;
    if !status.success() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::Other,
            "minimap2 failed during plasmid rescue read mapping",
        ));
    }
    Ok(())
}

fn mapped_read_names(paf_path: &Path) -> std::io::Result<HashSet<String>> {
    let reader = BufReader::new(File::open(paf_path)?);
    let mut mapped_reads = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        if let Some(read_name) = fields.first().filter(|name| !name.is_empty()) {
            mapped_reads.insert((*read_name).to_string());
        }
    }
    Ok(mapped_reads)
}

fn write_rescue_reads(
    filtered_reads: &Path,
    output_path: &Path,
    mapped_reads: &HashSet<String>,
) -> std::io::Result<RescueStats> {
    let mut reader = crate::utils::open_fastq_reader(filtered_reads)?;
    let mut writer = BufWriter::new(File::create(output_path)?);
    let mut stats = RescueStats::default();

    while let Some((header, sequence, plus, quality)) =
        crate::utils::read_fastq_record(reader.as_mut())?
    {
        let read_name = header
            .strip_prefix('@')
            .unwrap_or(&header)
            .split_whitespace()
            .next()
            .unwrap_or_default();
        if !mapped_reads.contains(read_name) {
            writeln!(writer, "{}\n{}\n{}\n{}", header, sequence, plus, quality)?;
            stats.rescued_read_count += 1;
            stats.rescued_read_bases += sequence.len() as u64;
        } else {
            stats.mapped_read_count += 1;
        }
        stats.input_read_count += 1;
    }
    writer.flush()?;
    Ok(stats)
}

fn fasta_stats(path: &Path) -> std::io::Result<(usize, u64)> {
    let reader = BufReader::new(File::open(path)?);
    let mut unitig_count = 0;
    let mut assembly_length = 0;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            unitig_count += 1;
        } else {
            assembly_length += line.trim().len() as u64;
        }
    }
    Ok((unitig_count, assembly_length))
}

#[cfg(test)]
mod tests {
    use super::{mapped_read_names, rescue_target_coverage, write_rescue_reads};
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn rescues_only_reads_without_any_alignment() {
        let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("target")
            .join(format!("rescue_reads_test_{}", std::process::id()));
        fs::create_dir_all(&directory).unwrap();
        let input = directory.join("reads.fq");
        let output = directory.join("rescue_reads.fq");
        let paf = directory.join("mapping.paf");
        fs::write(
            &input,
            "@mapped_full\nACGTACGT\n+\nIIIIIIII\n@mapped_partial\nACGTACGT\n+\nIIIIIIII\n@unmapped\nACGTACGT\n+\nIIIIIIII\n",
        )
        .unwrap();
        fs::write(
            &paf,
            "mapped_full\t8\t0\t8\t+\tunitig_0\t100\t0\t8\t8\t8\t60\n\
             mapped_partial\t8\t0\t1\t+\tunitig_1\t100\t0\t1\t1\t1\t60\n",
        )
        .unwrap();
        let mapped = mapped_read_names(&paf).unwrap();

        let stats = write_rescue_reads(&input, &output, &mapped).unwrap();
        assert_eq!(stats.rescued_read_count, 1);
        assert_eq!(stats.rescued_read_bases, 8);
        assert_eq!(stats.input_read_count, 3);
        assert_eq!(stats.mapped_read_count, 2);
        let result = fs::read_to_string(output).unwrap();
        assert!(result.contains("@unmapped"));
        assert!(!result.contains("@mapped_full"));
        assert!(!result.contains("@mapped_partial"));
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn rescue_coverage_budget_includes_all_rescued_bases() {
        assert_eq!(rescue_target_coverage(381_937_170, 4_787_666), 80);
        assert_eq!(rescue_target_coverage(20_000_000, 4_787_415), 5);
        assert_eq!(rescue_target_coverage(100, 0), u32::MAX);
    }
}
