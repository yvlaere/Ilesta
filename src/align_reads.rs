// module for filtering and anigning reads before assembly

use std::process::{Command};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write, Read};
use std::collections::HashSet;
use std::time::Instant;
use flate2::read::GzDecoder;

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

fn run_minimap2(query: &str, threads: usize, out_path: &str) -> std::io::Result<()> {
    let mut child = Command::new("minimap2")
        .args(&["-x", "ava-ont", "-t", &threads.to_string(), query, query, "-o", out_path])
        .spawn()?;

    let status = child.wait()?;
    assert!(status.success());

    Ok(())
}

pub fn align_reads(reads_fq: &str, threads: usize, output_paf: &str) -> std::io::Result<()> {
    println!("Computing read stats...");
    let output = File::create("temp_filtered.fq")?;
    let stats = filter_fastq(reads_fq, "temp_filtered.fq", 1000, 10.0)?;

    println!("Running minimap2...");
    run_minimap2("temp_filtered.fq", threads, output_paf)?;
    println!("Minimap2 finished. Alignments written to {}", output_paf);

    Ok(())
}