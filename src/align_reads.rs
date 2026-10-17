// Read filtering, genome size estimation and all-vs-all alignment before assembly

use crate::cli::ReadOpts;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// Bases in the first subsample used to estimate the genome size.  The subsample
/// is doubled until the estimated coverage is high enough, so small genomes
/// never pay for an all-vs-all alignment of the full data set.
const ESTIMATION_START_BASES: u64 = 50_000_000;
/// Minimum coverage of the subsample for a reliable genome size estimate.
const ESTIMATION_MIN_COVERAGE: f64 = 10.0;

pub struct AlignReadsResult {
    pub primary_reads: PathBuf,
    pub filtered_reads: PathBuf,
    pub filtered_read_count: usize,
    pub primary_read_count: usize,
}

/// Write reads that pass the length and quality filters; returns (name, length) per kept read.
fn filter_fastq(
    input: &Path,
    output: &Path,
    min_len: u32,
    min_q: f32,
) -> std::io::Result<Vec<(String, u32)>> {
    let mut reader = crate::utils::open_fastq_reader(input)?;
    let mut writer = BufWriter::new(File::create(output)?);
    let mut kept = Vec::new();
    while let Some((header, seq, plus, qual)) = crate::utils::read_fastq_record(reader.as_mut())? {
        let len = seq.len() as u32;
        let sum_q: u64 = qual.bytes().map(|b| b.saturating_sub(33) as u64).sum();
        if len >= min_len && sum_q as f32 / len as f32 >= min_q {
            kept.push((header.trim_start_matches('@').to_string(), len));
            write!(writer, "{header}\n{seq}\n{plus}\n{qual}\n")?;
        }
    }
    writer.flush()?;
    println!("Processed {} reads", kept.len());
    Ok(kept)
}

/// Keep the longest reads until `nr_bases` is reached; returns (bases, reads) written.
fn subsample_fastq(
    input: &Path,
    output: &Path,
    reads: &[(String, u32)],
    nr_bases: u64,
) -> std::io::Result<(u64, usize)> {
    let mut by_length: Vec<&(String, u32)> = reads.iter().collect();
    by_length.sort_by_key(|r| std::cmp::Reverse(r.1));
    let mut selected = HashSet::new();
    let mut total: u64 = 0;
    for (name, len) in by_length {
        selected.insert(name.as_str());
        total += *len as u64;
        if total >= nr_bases {
            break;
        }
    }
    println!("Selected {} reads with total length {} bases", selected.len(), total);

    let mut reader = crate::utils::open_fastq_reader(input)?;
    let mut writer = BufWriter::new(File::create(output)?);
    while let Some((header, seq, plus, qual)) = crate::utils::read_fastq_record(reader.as_mut())? {
        if selected.contains(header.trim_start_matches('@')) {
            write!(writer, "{header}\n{seq}\n{plus}\n{qual}\n")?;
        }
    }
    writer.flush()?;
    Ok((total, selected.len()))
}

/// All-vs-all alignment of a read file against itself.
fn align_all_vs_all(reads: &Path, paf: &Path, opts: &ReadOpts) -> std::io::Result<()> {
    println!("Running minimap2...");
    let mut args: Vec<&str> = opts.read_type.minimap2_args().to_vec();
    if let Some(f) = &opts.minimap2_f {
        args.extend(["-f", f]);
    }
    if let Some(u) = &opts.minimap2_u {
        args.extend(["-U", u]);
    }
    args.extend(["-K", &opts.minimap_batch_size]);
    crate::utils::run_minimap2(&args, opts.threads, reads, reads, paf)
}

/// Estimate coverage and genome size from a minimap2 all-vs-all PAF.
///
/// minimap2's `ava` presets report every overlapping read pair only once, so an
/// overlap adds aligned bases to *both* reads.
fn estimate_genome_size(paf: &Path) -> std::io::Result<Option<(f64, u64)>> {
    let mut seen = HashSet::new();
    let mut total_bases: u64 = 0;
    let mut total_aligned: u64 = 0;
    for line in BufReader::new(File::open(paf)?).lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        let n = |i: usize| f[i].parse::<u64>().ok();
        let (Some(qlen), Some(qs), Some(qe), Some(tlen), Some(ts), Some(te)) =
            (n(1), n(2), n(3), n(6), n(7), n(8))
        else {
            continue;
        };
        if seen.insert(f[0].to_string()) {
            total_bases += qlen;
        }
        if seen.insert(f[5].to_string()) {
            total_bases += tlen;
        }
        total_aligned += qe.saturating_sub(qs) + te.saturating_sub(ts);
    }
    if total_bases == 0 || total_aligned == 0 {
        return Ok(None);
    }
    let coverage = total_aligned as f64 / total_bases as f64;
    Ok(Some((coverage, (total_bases as f64 / coverage) as u64)))
}

pub fn align_reads(opts: &ReadOpts, out_dir: &Path) -> std::io::Result<AlignReadsResult> {
    println!("Computing read stats...");
    let paf = out_dir.join(&opts.paf);
    let filtered_all = out_dir.join("filtered_all.fq");
    let reads = filter_fastq(
        Path::new(&opts.reads_fq),
        &filtered_all,
        opts.min_read_length,
        opts.min_base_quality,
    )?;

    // `aligned_all`: the estimation alignment already covers every read, so it can double as the final one.
    let mut aligned_all = false;
    let genome_size = match opts.genome_size {
        Some(size) => {
            println!("Using provided genome size: {size} bases");
            size as u64
        }
        None => {
            println!("No genome size provided, will estimate from data");
            let subsample = out_dir.join("subsampled.fq");
            let mut bases = ESTIMATION_START_BASES;
            let estimate = loop {
                let (got, _) = subsample_fastq(&filtered_all, &subsample, &reads, bases)?;
                align_all_vs_all(&subsample, &paf, opts)?;
                let estimate = estimate_genome_size(&paf)?;
                aligned_all = got < bases;
                if aligned_all || estimate.is_some_and(|(cov, _)| cov >= ESTIMATION_MIN_COVERAGE) {
                    break estimate;
                }
                bases *= 2;
            };
            std::fs::remove_file(subsample)?;
            let (coverage, size) = estimate.ok_or_else(|| {
                std::io::Error::other("no read overlaps found, cannot estimate the genome size; provide --genome-size")
            })?;
            println!("Coverage: {coverage:.2}\nEstimated genome size: {size}");
            size
        }
    };

    let target_bases = genome_size.saturating_mul(opts.target_coverage as u64);
    let total_bases: u64 = reads.iter().map(|r| r.1 as u64).sum();
    let (primary, primary_read_count) = if aligned_all && target_bases >= total_bases {
        (filtered_all.clone(), reads.len())
    } else {
        let primary = out_dir.join("filtered.fq");
        let (_, n) = subsample_fastq(&filtered_all, &primary, &reads, target_bases)?;
        align_all_vs_all(&primary, &paf, opts)?;
        (primary, n)
    };
    println!("Final alignment finished. Alignments written to {}", paf.display());

    Ok(AlignReadsResult {
        primary_reads: primary,
        filtered_reads: filtered_all,
        filtered_read_count: reads.len(),
        primary_read_count,
    })
}