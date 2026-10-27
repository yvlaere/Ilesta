//! Rescue: reads that do not map to the primary assembly are assembled separately
//! (they are rarely covered well enough to survive the main graph) and the resulting graph is
//! merged into the main one.

use crate::alignment_filtering::Overlap;
use crate::cli::{OverlapOpts, ReadOpts};
use crate::create_overlap_graph::OverlapGraph;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

pub struct RescueResult {
    /// Rescued reads (down-sampled for assembly), to be appended to the assembly reads.
    pub reads: PathBuf,
    /// Nodes of the rescue graph, used to report which unitigs come from the rescue.
    pub node_ids: HashSet<String>,
}

/// Names of all reads with at least one alignment in a PAF.
fn mapped_read_names(paf: &Path) -> std::io::Result<HashSet<String>> {
    let mut names = HashSet::new();
    for line in BufReader::new(File::open(paf)?).lines() {
        if let Some(name) = line?.split('\t').next().filter(|n| !n.is_empty()) {
            names.insert(name.to_string());
        }
    }
    Ok(names)
}

/// Write the reads that are not in `mapped`; returns (reads written, bases written).
fn write_unmapped_reads(reads: &Path, out: &Path, mapped: &HashSet<String>) -> std::io::Result<(usize, u64)> {
    let mut reader = crate::utils::open_fastq_reader(reads)?;
    let mut writer = BufWriter::new(File::create(out)?);
    let (mut count, mut bases) = (0, 0);
    while let Some((header, seq, plus, qual)) = crate::utils::read_fastq_record(reader.as_mut())? {
        let name = header.trim_start_matches('@').split_whitespace().next().unwrap_or_default();
        if !mapped.contains(name) {
            write!(writer, "{header}\n{seq}\n{plus}\n{qual}\n")?;
            count += 1;
            bases += seq.len() as u64;
        }
    }
    writer.flush()?;
    Ok((count, bases))
}

/// Assemble the reads of `all_reads` that do not map to `assembly_fa` and add the result to
/// `graph` / `overlaps`. Returns `None` if there are no such reads.
pub fn run_rescue_stage(
    all_reads: &Path,
    assembly_fa: &Path,
    out_dir: &Path,
    reads: &ReadOpts,
    overlap: &OverlapOpts,
    graph: &mut OverlapGraph,
    overlaps: &mut HashMap<(usize, usize), Overlap>,
) -> Result<Option<RescueResult>, Box<dyn std::error::Error>> {
    let mapping_paf = out_dir.join("rescue_mapping.paf");
    crate::utils::run_minimap2(reads.read_type.mapping_args(), reads.threads, assembly_fa, all_reads, &mapping_paf)?;
    let rescue_fq = out_dir.join("rescue_reads.fq");
    let (count, bases) = write_unmapped_reads(all_reads, &rescue_fq, &mapped_read_names(&mapping_paf)?)?;
    println!("Reads not mapping to the primary assembly: {count} ({bases} bases)");
    if count == 0 {
        println!("No poorly represented reads found; skipping rescue assembly.");
        return Ok(None);
    }

    // assemble the unmapped reads with a coverage budget that includes all of them
    let assembly_len: u64 = BufReader::new(File::open(assembly_fa)?)
        .lines()
        .map_while(Result::ok)
        .filter(|l| !l.starts_with('>'))
        .map(|l| l.trim().len() as u64)
        .sum();
    let target_coverage = rescue_target_coverage(bases, assembly_len);
    println!("Rescue assembly target: {target_coverage}x primary assembly length");
    let dir = out_dir.join("rescue_graph");
    std::fs::create_dir_all(&dir)?;
    let rescue_opts = ReadOpts {
        reads_fq: rescue_fq.to_string_lossy().into_owned(),
        threads: reads.threads.clamp(1, 4),
        paf: "rescue.paf".into(),
        min_read_length: reads.min_read_length.saturating_sub(1).clamp(1, 500),
        genome_size: Some(assembly_len.clamp(1, u32::MAX as u64) as u32),
        target_coverage,
        minimap2_f: Some("0.001".into()),
        minimap2_K: reads.minimap2_K.clone(),
        minimap2_U: reads.minimap2_U.clone(),
        ..reads.clone()
    };
    let rescued = crate::align_reads::align_reads(&rescue_opts, &dir)?;
    let rescue_overlaps = crate::alignment_filtering::run_alignment_filtering(&dir.join("rescue.paf"), overlap)?.overlaps;
    let rescue_graph = crate::create_overlap_graph::run_create_overlap_graph(&rescue_overlaps)?;

    // merge: overlap keys are shifted so that they cannot collide with the main ones
    let key_offset = overlaps.keys().flat_map(|&(q, t)| [q, t]).max().map_or(0, |m| m + 1);
    overlaps.extend(rescue_overlaps.into_iter().map(|((q, t), o)| ((q + key_offset, t + key_offset), o)));
    let node_ids = rescue_graph.nodes.keys().cloned().collect();
    for id in rescue_graph.nodes.keys() {
        graph.add_node(id.clone());
    }
    for (id, node) in rescue_graph.nodes {
        for e in node.edges {
            graph.add_edge(&id, &e.target_id, e.edge_len, e.overlap_len, e.identity);
        }
    }
    Ok(Some(RescueResult { reads: rescued.primary_reads, node_ids }))
}

/// Coverage (x the primary assembly length) needed to include all rescued bases.
fn rescue_target_coverage(rescued_bases: u64, assembly_len: u64) -> u32 {
    if assembly_len == 0 {
        return u32::MAX;
    }
    rescued_bases.div_ceil(assembly_len).clamp(1, u32::MAX as u64) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rescues_only_reads_without_any_alignment() {
        let dir = std::env::temp_dir().join(format!("rescue_reads_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (input, output, paf) = (dir.join("reads.fq"), dir.join("rescue_reads.fq"), dir.join("mapping.paf"));
        let rec = |n: &str| format!("@{n}\nACGTACGT\n+\nIIIIIIII\n");
        std::fs::write(&input, rec("mapped_full") + &rec("mapped_partial") + &rec("unmapped")).unwrap();
        std::fs::write(
            &paf,
            "mapped_full\t8\t0\t8\t+\tunitig_0\t100\t0\t8\t8\t8\t60\nmapped_partial\t8\t0\t1\t+\tunitig_1\t100\t0\t1\t1\t1\t60\n",
        )
        .unwrap();
        let mapped = mapped_read_names(&paf).unwrap();
        assert_eq!(write_unmapped_reads(&input, &output, &mapped).unwrap(), (1, 8));
        let result = std::fs::read_to_string(output).unwrap();
        assert!(result.contains("@unmapped") && !result.contains("@mapped"));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn rescue_coverage_budget_includes_all_rescued_bases() {
        assert_eq!(rescue_target_coverage(381_937_170, 4_787_666), 80);
        assert_eq!(rescue_target_coverage(20_000_000, 4_787_415), 5);
        assert_eq!(rescue_target_coverage(100, 0), u32::MAX);
    }
}