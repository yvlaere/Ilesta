// module for filtering and anigning reads before assembly

use std::process::{Command};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write, Read};
use std::collections::HashSet;
use std::time::Instant;

#[derive(Clone)]
struct ReadStats {
    name: String,
    len: usize,
    mean_q: f32,
}

fn filter_fastq(input_path: &str, output_path: &str, min_len: usize, min_q: f32) -> std::io::Result<Vec<ReadStats>> {
    
    // time for debugging
    let start_time = Instant::now();

    // read the file
    let file = File::open(input_path)?;
    let reader = BufReader::new(file);

    let output = File::create(output_path)?;
    let mut writer = BufWriter::new(output);

    // buffers reused for every read
    let mut header = String::with_capacity(100);
    let mut seq = String::with_capacity(10_000);
    let mut plus = String::with_capacity(10);
    let mut qual = String::with_capacity(10_000);

    let mut results = Vec::new();
    let mut reader_lines = reader.lines();

    let duration = start_time.elapsed();
    println!("File {} opened and reader initialized in {:?}", input_path, duration);

    loop {
        // clear buffers
        header.clear();
        seq.clear();
        plus.clear();
        qual.clear();

        // Read 4 lines per FASTQ record
        if reader_lines.next().map(|l| { header.push_str(&l.unwrap()); }).is_none() {
            break; // EOF
        }
        reader_lines.next().map(|l| seq.push_str(&l.unwrap()));
        reader_lines.next().map(|l| plus.push_str(&l.unwrap()));
        reader_lines.next().map(|l| qual.push_str(&l.unwrap()));

        // parse read name
        let name = if header.starts_with('@') { &header[1..] } else { &header };

        // get read length
        let len = seq.len();

        let sum_q: u64 = qual.bytes().map(|b| (b.saturating_sub(33)) as u64).sum();
        let mean_q = sum_q as f32 / len as f32;

        if len >= min_len && mean_q >= min_q {
            results.push(ReadStats {
                name: name.to_string(),
                len,
                mean_q,
            });

            // Write FASTQ record
            writeln!(writer, "{}", header)?;
            writeln!(writer, "{}", seq)?;
            writeln!(writer, "{}", plus)?;
            writeln!(writer, "{}", qual)?;
        }
    }

    writer.flush()?;
    println!(
        "Processed {} reads in {:?}",
        results.len(),
        start_time.elapsed()
    );

    Ok(results)
}

fn subsample_fastq(input_path: &str, output_path: &str, stats: &Vec<ReadStats>, nr_bases: usize) -> std::io::Result<()> {
    
    // time for debugging
    let start_time = Instant::now();

    // keep the n longest reads untill we reach nr_bases
    let mut sorted_stats = stats.clone();
    sorted_stats.sort_by_key(|r| std::cmp::Reverse(r.len));
    let mut selected_reads = HashSet::new();
    let mut total_len = 0;

    for r in sorted_stats {
        if total_len + r.len >= nr_bases {
            break;
        }
        total_len += r.len;
        selected_reads.insert(r.name);
    }

    println!(
        "Selected {} reads with total length {} bases in {:?}",
        selected_reads.len(),
        total_len,
        start_time.elapsed()
    );

    let file = File::open(input_path)?;
    let reader = BufReader::new(file);

    let output = File::create(output_path)?;
    let mut writer = BufWriter::new(output);

    let mut header = String::with_capacity(100);
    let mut seq = String::with_capacity(10_000);
    let mut plus = String::with_capacity(10);
    let mut qual = String::with_capacity(10_000);

    let mut reader_lines = reader.lines();

    let duration = start_time.elapsed();
    println!("File {} opened and reads identified in {:?}", input_path, duration);

    loop {
        header.clear();
        seq.clear();
        plus.clear();
        qual.clear();

        if reader_lines.next().map(|l| { header.push_str(&l.unwrap()); }).is_none() {
            break; // EOF
        }

        let name = if header.starts_with('@') { &header[1..] } else { &header };

        if !selected_reads.contains(name) {
            // skip this read
            reader_lines.next();
            reader_lines.next();
            reader_lines.next();
            continue;
        }
        reader_lines.next().map(|l| seq.push_str(&l.unwrap()));
        reader_lines.next().map(|l| plus.push_str(&l.unwrap()));
        reader_lines.next().map(|l| qual.push_str(&l.unwrap()));

        writeln!(writer, "{}", header)?;
        writeln!(writer, "{}", seq)?;
        writeln!(writer, "{}", plus)?;
        writeln!(writer, "{}", qual)?;
    }

    writer.flush()?;
    println!(
        "Processed reads in {:?}",
        start_time.elapsed()
    );
    Ok(())
}

fn run_minimap2(query: &str, threads: usize, out_path: &str) -> std::io::Result<()> {

    // time for debugging
    let start_time = Instant::now();

    let mut child = Command::new("minimap2")
        .args(&["-x", "ava-ont", "-t", &threads.to_string(), query, query, "-o", out_path])
        .spawn()?;

    let status = child.wait()?;
    assert!(status.success());

    println!("Minimap2 finished in {:?}", start_time.elapsed());

    Ok(())
}

fn estimate_genome_size(paf_path: &str) -> std::io::Result<usize> {
    let file = File::open(paf_path)?;
    let reader = BufReader::new(file);

    let mut seen = HashSet::new();
    let mut total_bases: u64 = 0;
    let mut total_aligned: u64 = 0;

    for line in reader.lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();

        let qname = fields[0];
        let qlen: u64 = fields[1].parse().unwrap();
        let qstart: u64 = fields[2].parse().unwrap();
        let qend: u64 = fields[3].parse().unwrap();

        // count each read once
        if seen.insert(qname.to_string()) {
            total_bases += qlen;
        }

        total_aligned += qend - qstart;
    }

    let coverage = total_aligned as f64 / total_bases as f64;
    let genome_size = total_bases as f64 / coverage;

    println!("Total bases: {}", total_bases);
    println!("Coverage: {:.2}", coverage);
    println!("Estimated genome size: {:.0}", genome_size);

    Ok(genome_size as usize)
}

pub fn align_reads(reads_fq: &str, threads: usize, output_paf: &str) -> std::io::Result<()> {
    println!("Computing read stats...");
    let basic_filtering = "temp_filtered.fq";
    let stats = filter_fastq(reads_fq, basic_filtering, 1000, 10.0)?;

    let subsampled_output = "temp_subsampled.fq";
    subsample_fastq(basic_filtering, subsampled_output, &stats, 500_000_000)?;

    println!("Running minimap2...");
    run_minimap2(subsampled_output, threads, output_paf)?;
    println!("Minimap2 finished. Alignments written to {}", output_paf);

    let genome_size = estimate_genome_size(output_paf)?;

    // second round of subsampling based on the estimated genome size
    let subsampled_output_2 = "temp_subsampled_2.fq";
    subsample_fastq(reads_fq, subsampled_output_2, &stats, genome_size * 50)?;
    run_minimap2(subsampled_output_2, threads, output_paf)?;
    println!("Second round of minimap2 finished. Alignments written to {}", output_paf);

    Ok(())
}