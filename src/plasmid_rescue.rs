use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;

#[derive(Debug, Default)]
pub struct RescueStats {
    pub input_read_count: usize,
    pub mapped_read_count: usize,
    pub rescued_read_count: usize,
    pub rescued_read_bases: u64,
    pub rescued_unitig_count: usize,
    pub rescued_assembly_length: u64,
}

pub fn run_rescue_stage(
    filtered_reads: &Path,
    improved_assembly: &Path,
    output_dir: &Path,
    threads: usize,
    read_type: crate::cli::ReadType,
    min_base_quality: f32,
    rescue_min_read_length: u32,
) -> std::io::Result<RescueStats> {
    let paf_path = output_dir.join("rescue_mapping.paf");
    map_reads(
        filtered_reads,
        improved_assembly,
        &paf_path,
        threads,
        read_type,
    )?;
    let mapped_reads = mapped_read_names(&paf_path)?;
    let rescue_path = output_dir.join("rescue_reads.fq");
    let stats = write_rescue_reads(filtered_reads, &rescue_path, &mapped_reads)?;
    println!("Reads considered for rescue: {}", stats.input_read_count);
    println!(
        "Reads mapped to primary assembly (excluded): {}",
        stats.mapped_read_count
    );
    println!("Reads written to rescue_reads.fq: {}", stats.rescued_read_count);
    println!("Rescue read bases: {}", stats.rescued_read_bases);

    if stats.rescued_read_count == 0 {
        println!("No poorly represented reads found; skipping rescue assembly.");
        return Ok(stats);
    }

    let assembly_dir = output_dir.join("rescue_assembly");
    std::fs::create_dir_all(&assembly_dir)?;
    let genome_size = stats
        .rescued_read_bases
        .div_ceil(50)
        .max(1)
        .min(u32::MAX as u64) as u32;
    let output = Command::new(std::env::current_exe()?)
        .arg("assemble")
        .arg("--reads-fq")
        .arg(&rescue_path)
        .arg("--output-dir")
        .arg(&assembly_dir)
        .arg("--output-prefix")
        .arg("rescued")
        .arg("--threads")
        .arg(threads.to_string())
        .arg("--read-type")
        .arg(read_type_name(read_type))
        .arg("--min-read-length")
        .arg(rescue_min_read_length.to_string())
        .arg("--min-base-quality")
        .arg(min_base_quality.to_string())
        .arg("--genome-size")
        .arg(genome_size.to_string())
        .output()?;
    if !output.status.success() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!(
                "rescue assembly subprocess failed with {}: {}{}",
                output.status,
                String::from_utf8_lossy(&output.stderr),
                String::from_utf8_lossy(&output.stdout)
            ),
        ));
    }

    let rescued_fasta = assembly_dir.join("rescued.fa");
    (stats.rescued_unitig_count, stats.rescued_assembly_length) = fasta_stats(&rescued_fasta)?;
    Ok(stats)
}

fn read_type_name(read_type: crate::cli::ReadType) -> &'static str {
    match read_type {
        crate::cli::ReadType::Ont => "ont",
        crate::cli::ReadType::PbClr => "pb-clr",
        crate::cli::ReadType::PbHifi => "pb-hifi",
    }
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
    use super::{mapped_read_names, write_rescue_reads};
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn rescues_reads_below_half_aligned_fraction_and_unmapped_reads() {
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
}
