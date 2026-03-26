// module for filtering and anigning reads before assembly

use std::process::{Command};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write, Read};
use std::collections::HashSet;
use std::time::Instant;

#[derive(Clone)]
struct ReadStats {
    name: String,
    len: u32,
    mean_q: f32,
}

fn filter_fastq(input_path: &std::path::Path, output_path: &std::path::Path, min_len: u32, min_q: f32) -> std::io::Result<Vec<ReadStats>> {
    
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
        let len = seq.len() as u32;

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
        "Processed {} reads",
        results.len()
    );

    Ok(results)
}

fn subsample_fastq(input_path: &std::path::Path, output_path: &std::path::Path, stats: &Vec<ReadStats>, nr_bases: u32) -> std::io::Result<u32> {
    
    // time for debugging
    let start_time = Instant::now();

    // keep the n longest reads untill we reach nr_bases
    let mut sorted_stats = stats.clone();
    sorted_stats.sort_by_key(|r| std::cmp::Reverse(r.len));
    let mut selected_reads = HashSet::new();
    let mut total_len = 0;

    for r in sorted_stats {
        if total_len + r.len >= nr_bases {
            total_len += r.len;
            selected_reads.insert(r.name);
            break;
        }
        total_len += r.len;
        selected_reads.insert(r.name);
    }

    println!(
        "Selected {} reads with total length {} bases",
        selected_reads.len(),
        total_len
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
    Ok(total_len)
}

fn run_minimap2(query: &std::path::Path, threads: usize, out_path: &std::path::Path) -> std::io::Result<()> {

    // time for debugging
    let start_time = Instant::now();

    let mut child = Command::new("minimap2")
        .arg("-x").arg("ava-ont").arg("-t").arg(&threads.to_string()).arg(query).arg(query).arg("-o").arg(out_path)
        .spawn()?;

    let status = child.wait()?;
    assert!(status.success());

    //println!("Minimap2 finished in {:?}", start_time.elapsed());

    Ok(())
}

fn estimate_genome_size(paf_path: &std::path::Path) -> std::io::Result<(f64, u32)> {
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

    Ok((coverage as f64, genome_size as u32))
}

pub fn align_reads(reads_fq: &std::path::Path, threads: usize, output_paf: &std::path::Path, out_dir: &std::path::Path, min_read_length: u32, min_base_quality: f32, input_genome_size: Option<u32>) -> std::io::Result<std::path::PathBuf> {
    
    println!("Computing read stats...");

    // compute read stats and filter reads
    let basic_filtering = "basic_filtered.fq";
    let basic_filtering_path = out_dir.join(basic_filtering);
    let stats = filter_fastq(reads_fq, &basic_filtering_path, min_read_length, min_base_quality)?;

    let mut genome_size = 0;

    match input_genome_size {
        Some(size) => {
            println!("Using provided genome size: {} bases", size);
            genome_size = size;
        },
        None => {
            println!("No genome size provided, will estimate from data");
            // subsample reads for genome size estimation
            let subsampled = "subsampled.fq";
            let subsampled_path = out_dir.join(subsampled);
            let subsampled_bases = subsample_fastq(&basic_filtering_path, &subsampled_path, &stats, 500_000_000)?;

            // align reads for genome size estimation
            println!("Running minimap2...");
            run_minimap2(&subsampled_path, threads, output_paf)?;
            println!("Genome size estimation alignment finished. Alignments written to {}", output_paf.display());

            // estimate genome size
            let mut coverage = 0 as f64;
            (coverage, genome_size) = estimate_genome_size(output_paf)?;

            let mut subsampling_rounds = 1;
            let mut previous_subsampled_bases = 0;

            while coverage < 10.0 {
                // subsample reads for genome size estimation
                let subsampled = "subsampled.fq";
                let subsampled_path = out_dir.join(subsampled);
                let subsampled_bases = subsample_fastq(&basic_filtering_path, &subsampled_path, &stats, subsampling_rounds * 500_000_000)?;

                if subsampled_bases == previous_subsampled_bases {
                    // coverage too low, throw an error
                    eprintln!("Coverage is too low. Stopping.");
                    break;
                }

                // align reads for genome size estimation
                println!("Running minimap2...");
                run_minimap2(&subsampled_path, threads, output_paf)?;
                println!("Genome size estimation alignment finished. Alignments written to {}", output_paf.display());

                // estimate genome size
                (coverage, genome_size) = estimate_genome_size(output_paf)?;

                previous_subsampled_bases = subsampled_bases;
                subsampling_rounds += 1;
            }
        }
    }
    
    // final round of subsampling based on the estimated genome size
    let subsampled_output = "filtered.fq";
    let subsampled_output_path = out_dir.join(subsampled_output);
    subsample_fastq(reads_fq, &subsampled_output_path, &stats, genome_size * 50)?;
    run_minimap2(&subsampled_output_path, threads, output_paf)?;
    println!("Final alignment finished. Alignments written to {}", output_paf.display());

    let path = std::path::PathBuf::from(subsampled_output_path);

    Ok(path)
}