// Read filtering, genome size estimation and all-vs-all alignment before assembly

use crate::cli::ReadOpts;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// Bases in the first subsample used to estimate the genome size.  The subsample
/// is doubled until the estimated coverage is high enough, so small genomes
/// never pay for an all-vs-all alignment of the full data set.
const ESTIMATION_START_BASES: u64 = 200_000_000;
/// Minimum coverage of the subsample for a reliable genome size estimate.
const ESTIMATION_MIN_COVERAGE: f64 = 20.0;
/// Alignments shorter than this (on either read) do not count towards the coverage estimate.
const ESTIMATION_MIN_ALIGNMENT: u64 = 2_000;

pub struct AlignReadsResult {
    pub primary_reads: PathBuf,
    pub filtered_reads: PathBuf,
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
    let mut header = String::new();
    let mut seq = String::new();
    let mut plus = String::new();
    let mut qual = String::new();
    let mut line = String::new();
    while crate::utils::read_fastq_record_into(
        reader.as_mut(),
        &mut header,
        &mut seq,
        &mut plus,
        &mut qual,
        &mut line,
    )? {
        let len = seq.len() as u32;
        if len < min_len {
            continue;
        }
        let sum_q: u64 = qual.bytes().map(|b| b.saturating_sub(33) as u64).sum();
        if sum_q as f32 / len as f32 >= min_q {
            kept.push((header.trim_start_matches('@').to_string(), len));
            writer.write_all(header.as_bytes())?;
            writer.write_all(b"\n")?;
            writer.write_all(seq.as_bytes())?;
            writer.write_all(b"\n")?;
            writer.write_all(plus.as_bytes())?;
            writer.write_all(b"\n")?;
            writer.write_all(qual.as_bytes())?;
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
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
    write_selected_reads(input, output, &selected)?;
    Ok((total, selected.len()))
}

/// Copy the reads whose name is in `selected` from `input` to `output`.
fn write_selected_reads(input: &Path, output: &Path, selected: &HashSet<&str>) -> std::io::Result<()> {
    let mut reader = crate::utils::open_fastq_reader(input)?;
    let mut writer = BufWriter::new(File::create(output)?);
    while let Some((header, seq, plus, qual)) = crate::utils::read_fastq_record(reader.as_mut())? {
        if selected.contains(header.trim_start_matches('@')) {
            write!(writer, "{header}\n{seq}\n{plus}\n{qual}\n")?;
        }
    }
    writer.flush()
}

/// A read is removed as fold-back when its hairpin (both arms) covers at least this fraction of
/// it. A read with a small hairpin is mostly genuine sequence and does no harm.
const MIN_FOLD_BACK_FRACTION: f64 = 0.25;
/// Reads whose hairpins are looked for are processed in chunks of about this many bases.
const FOLD_BACK_CHUNK_BASES: usize = 64_000_000;

/// For the reads in `wanted`: the fraction of the read covered by its hairpin, for those that have one.
fn find_hairpins(input: &Path, wanted: &HashSet<&str>, threads: usize) -> std::io::Result<HashMap<String, f64>> {
    let mut found = HashMap::new();
    let mut chunk: Vec<(String, String)> = Vec::new();
    let mut chunk_bases = 0;
    let process = |chunk: &mut Vec<(String, String)>, found: &mut HashMap<String, f64>| {
        let per_thread = chunk.len().div_ceil(threads.max(1)).max(1);
        std::thread::scope(|scope| {
            let handles: Vec<_> = chunk
                .chunks(per_thread)
                .map(|part| {
                    scope.spawn(move || {
                        part.iter()
                            .filter_map(|(name, seq)| {
                                crate::fold_back::find_hairpin(seq.as_bytes())
                                    .map(|h| (name.clone(), h.span as f64 / seq.len() as f64))
                            })
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            for handle in handles {
                found.extend(handle.join().expect("hairpin search panicked"));
            }
        });
        chunk.clear();
    };
    let mut reader = crate::utils::open_fastq_reader(input)?;
    let (mut header, mut seq, mut plus, mut qual, mut line) =
        (String::new(), String::new(), String::new(), String::new(), String::new());
    while crate::utils::read_fastq_record_into(reader.as_mut(), &mut header, &mut seq, &mut plus, &mut qual, &mut line)? {
        let name = header.trim_start_matches('@');
        if wanted.contains(name) {
            chunk_bases += seq.len();
            chunk.push((name.to_string(), std::mem::take(&mut seq)));
            if chunk_bases >= FOLD_BACK_CHUNK_BASES {
                process(&mut chunk, &mut found);
                chunk_bases = 0;
            }
        }
    }
    process(&mut chunk, &mut found);
    Ok(found)
}

/// Keep the longest reads until `target_bases` is reached, leaving out fold-back reads. Fold-back
/// reads are replaced by the next longest reads, so the coverage of the selection does not drop.
/// Only the reads that can end up in the selection are searched for hairpins. Returns the number
/// of reads selected.
fn select_longest_reads(
    input: &Path,
    output: &Path,
    reads: &[(String, u32)],
    target_bases: u64,
    threads: usize,
    remove_fold_back: bool,
) -> std::io::Result<usize> {
    if !remove_fold_back {
        return Ok(subsample_fastq(input, output, reads, target_bases)?.1);
    }
    let mut by_length: Vec<&(String, u32)> = reads.iter().collect();
    by_length.sort_by_key(|r| std::cmp::Reverse(r.1));

    let mut checked = 0; // by_length[..checked] has been searched
    let mut hairpins: HashMap<String, f64> = HashMap::new();
    let mut budget = target_bases.saturating_add(target_bases / 4);
    let (selected, removed, removed_bases, small) = loop {
        // search the next batch of candidates
        let mut end = checked;
        let mut batch_bases = 0u64;
        while end < by_length.len() && batch_bases < budget {
            batch_bases += by_length[end].1 as u64;
            end += 1;
        }
        let wanted: HashSet<&str> = by_length[checked..end].iter().map(|r| r.0.as_str()).collect();
        if !wanted.is_empty() {
            hairpins.extend(find_hairpins(input, &wanted, threads)?);
        }
        checked = end;

        // walk the reads by length, skipping fold-back reads
        let (mut selected, mut bases, mut removed, mut removed_bases, mut small) = (HashSet::new(), 0u64, 0, 0u64, 0);
        let mut reached = false;
        for r in &by_length[..checked] {
            match hairpins.get(&r.0) {
                Some(&f) if f >= MIN_FOLD_BACK_FRACTION => {
                    removed += 1;
                    removed_bases += r.1 as u64;
                    continue;
                }
                Some(_) => small += 1,
                None => {}
            }
            selected.insert(r.0.as_str());
            bases += r.1 as u64;
            if bases >= target_bases {
                reached = true;
                break;
            }
        }
        if reached || checked == by_length.len() {
            break (selected, removed, removed_bases, small);
        }
        budget = budget.saturating_add(budget / 2);
    };
    println!(
        "Fold-back reads left out of the selection: {removed} ({removed_bases} bases); {small} reads with smaller hairpins kept"
    );
    write_selected_reads(input, output, &selected)?;
    Ok(selected.len())
}

/// All-vs-all alignment of a read file against itself.
fn align_all_vs_all(reads: &Path, paf: &Path, opts: &ReadOpts) -> std::io::Result<()> {
    let mut args: Vec<&str> = opts.read_type.minimap2_args().to_vec();
    if let Some(f) = &opts.minimap2_f {
        args.extend(["-f", f]);
    }
    if let Some(u) = &opts.minimap2_u {
        args.extend(["-U", u]);
    }
    args.extend(["-K", &opts.minimap_batch_size]);
    crate::utils::timed_step("Run minimap2 all-vs-all alignment", || {
        crate::utils::run_minimap2(&args, opts.threads, reads, reads, paf)
    })
}

/// Estimate coverage and genome size from a minimap2 all-vs-all PAF.
///
/// minimap2's `ava` presets report every overlapping read pair only once, so an overlap adds
/// aligned bases to *both* reads. The depth of a read is its aligned bases divided by its
/// length (about the coverage minus one); the coverage is the **median** over the reads.
/// A mean is dominated by high-copy elements (multi-copy plasmids, concatemer reads) and
/// repeats: in one benchmark dataset 10% of the longest reads had a depth of 120x against 11x
/// for the chromosome, which made the mean-based genome size 4x too small.
fn estimate_genome_size(paf: &Path) -> std::io::Result<Option<(f64, u64)>> {
    Ok(estimate_from_paf(BufReader::new(File::open(paf)?)))
}

fn estimate_from_paf(paf: impl BufRead) -> Option<(f64, u64)> {
    // read name -> (length, aligned bases)
    let mut reads: HashMap<String, (u64, u64)> = HashMap::new();
    let mut add = |name: &str, len: u64, aligned: u64| {
        if let Some(r) = reads.get_mut(name) {
            r.1 += aligned;
        } else {
            reads.insert(name.to_string(), (len, aligned));
        }
    };
    for line in paf.lines().map_while(Result::ok) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        let n = |i: usize| f[i].parse::<u64>().ok();
        let (Some(qlen), Some(qs), Some(qe), Some(tlen), Some(ts), Some(te)) = (n(1), n(2), n(3), n(6), n(7), n(8))
        else {
            continue;
        };
        let (qa, ta) = (qe.saturating_sub(qs), te.saturating_sub(ts));
        // short alignments are mostly repeats and spurious hits
        if qa.min(ta) < ESTIMATION_MIN_ALIGNMENT {
            continue;
        }
        add(f[0], qlen, qa);
        add(f[5], tlen, ta);
    }
    if reads.is_empty() {
        return None;
    }
    let total_bases: u64 = reads.values().map(|r| r.0).sum();
    let mut depths: Vec<f64> = reads.values().map(|&(len, aligned)| aligned as f64 / len.max(1) as f64).collect();
    depths.sort_by(f64::total_cmp);
    // a read at coverage c overlaps (c - 1) others at each position
    let coverage = depths[depths.len() / 2] + 1.0;
    Some((coverage, (total_bases as f64 / coverage) as u64))
}

/// Names and lengths of all reads of a FASTQ file, without writing a copy.
fn scan_fastq(input: &Path) -> std::io::Result<Vec<(String, u32)>> {
    let mut reader = crate::utils::open_fastq_reader(input)?;
    let mut reads = Vec::new();
    let (mut header, mut seq, mut plus, mut qual, mut line) =
        (String::new(), String::new(), String::new(), String::new(), String::new());
    while crate::utils::read_fastq_record_into(reader.as_mut(), &mut header, &mut seq, &mut plus, &mut qual, &mut line)? {
        if !seq.is_empty() {
            reads.push((header.trim_start_matches('@').to_string(), seq.len() as u32));
        }
    }
    Ok(reads)
}

/// Read filtering. Returns the file that holds the reads that passed, and their names and lengths.
///
/// With filtering disabled a plain FASTQ is used as it is (it is only scanned for names and
/// lengths). A gzipped input is always written out as a plain copy: every later step reads the
/// reads again (genome size estimation, selection, rescue), and decompressing 3-7 Gb of FASTQ each
/// time takes far longer than reading a plain copy.
fn prepare_reads(opts: &ReadOpts, out_dir: &Path) -> std::io::Result<(PathBuf, Vec<(String, u32)>)> {
    let input = Path::new(&opts.reads_fq);
    if opts.min_read_length == 0 && opts.min_base_quality <= 0.0 && !crate::utils::is_gzip(input)? {
        return Ok((input.to_path_buf(), scan_fastq(input)?));
    }
    let filtered_all = out_dir.join("filtered_all.fq");
    let reads = filter_fastq(input, &filtered_all, opts.min_read_length, opts.min_base_quality)?;
    Ok((filtered_all, reads))
}

pub fn align_reads(opts: &ReadOpts, out_dir: &Path) -> std::io::Result<AlignReadsResult> {
    let paf = out_dir.join(&opts.paf);
    let (filtered_all, reads) = crate::utils::timed_step("Filter reads", || prepare_reads(opts, out_dir))?;
    let filtered_bases: u64 = reads.iter().map(|r| r.1 as u64).sum();
    println!("Filtered reads: {} ({} bases)", reads.len(), filtered_bases);

    // `aligned_all`: the estimation alignment already covers every read, so it can double as the final one.
    let mut aligned_all = false;
    let genome_size = crate::utils::timed_step(
        "Determine genome size",
        || -> std::io::Result<u64> {
            match opts.genome_size {
                Some(size) => Ok(size as u64),
                None => {
                    let subsample = out_dir.join("subsampled.fq");
                    let mut bases = ESTIMATION_START_BASES;
                    let estimate = loop {
                        let (got, _) = subsample_fastq(&filtered_all, &subsample, &reads, bases)?;
                        align_all_vs_all(&subsample, &paf, opts)?;
                        let estimate = estimate_genome_size(&paf)?;
                        aligned_all = got < bases;
                        if aligned_all
                            || estimate.is_some_and(|(cov, _)| cov >= ESTIMATION_MIN_COVERAGE)
                        {
                            break estimate;
                        }
                        bases *= 2;
                    };
                    std::fs::remove_file(subsample)?;
                    let (coverage, size) = estimate.ok_or_else(|| {
                    std::io::Error::other(
                        "no read overlaps found, cannot estimate the genome size; provide --genome-size",
                    )
                })?;
                    println!("Estimated coverage: {coverage:.2}x");
                    Ok(size)
                }
            }
        },
    )?;
    println!(
        "Genome size: {} bases{}",
        genome_size,
        if opts.genome_size.is_some() {
            " (provided)"
        } else {
            " (estimated)"
        }
    );

    let target_bases = genome_size.saturating_mul(opts.target_coverage as u64);
    let (primary, primary_read_count) = crate::utils::timed_step(
        "Select primary reads and create final alignment",
        || -> std::io::Result<(PathBuf, usize)> {
            let primary = out_dir.join("filtered.fq");
            let selected = select_longest_reads(
                &filtered_all,
                &primary,
                &reads,
                target_bases,
                opts.threads,
                !opts.keep_fold_back_reads,
            )?;
            if aligned_all && selected == reads.len() {
                // the estimation alignment already covers exactly these reads
                std::fs::remove_file(&primary)?;
                Ok((filtered_all.clone(), reads.len()))
            } else {
                align_all_vs_all(&primary, &paf, opts)?;
                Ok((primary, selected))
            }
        },
    )?;
    println!("Primary assembly reads: {primary_read_count}");
    println!("All-vs-all alignments: {}", paf.display());

    Ok(AlignReadsResult {
        primary_reads: primary,
        filtered_reads: filtered_all,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn rec(q: &str, t: &str, len: u64) -> String {
        format!("{q}\t10000\t0\t{len}\t+\t{t}\t10000\t0\t{len}\t{len}\t{len}\t60\n")
    }

    /// 40 reads of 10 kb tiling a genome with 4 neighbours each side (depth ~ 9), plus
    /// high-copy "concatemer" reads that overlap each other all the time.
    fn paf(with_hub: bool) -> String {
        let mut s = String::new();
        for i in 0..40 {
            for j in 1..=4 {
                if i + j < 40 {
                    s += &rec(&format!("r{i}"), &format!("r{}", i + j), 5_000);
                }
            }
        }
        if with_hub {
            for i in 0..8 {
                for j in 0..8 {
                    for _ in 0..30 {
                        if i != j {
                            s += &rec(&format!("h{i}"), &format!("h{j}"), 9_000);
                        }
                    }
                }
            }
        }
        s
    }

    #[test]
    fn genome_size_estimate_ignores_high_copy_reads() {
        let (cov, size) = estimate_from_paf(Cursor::new(paf(false))).unwrap();
        let (hub_cov, hub_size) = estimate_from_paf(Cursor::new(paf(true))).unwrap();
        assert!((cov - hub_cov).abs() < 1.0, "{cov} vs {hub_cov}");
        // 40 reads of 10 kb at ~4.5x -> ~90 kb; the hub reads add 80 kb of reads, not of genome
        assert!((size as f64 - 90_000.0).abs() < 25_000.0, "{size}");
        assert!(hub_size < 200_000, "{hub_size}");
    }

    #[test]
    fn short_alignments_are_ignored() {
        assert!(estimate_from_paf(Cursor::new(rec("a", "b", 500))).is_none());
    }
}