use crate::create_overlap_graph::OverlapGraph;
use flate2::read::MultiGzDecoder;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

/// Run a pipeline step with timestamped start, completion, and elapsed-time output.
pub fn timed_step<T, E: std::fmt::Display>(
    label: &str,
    operation: impl FnOnce() -> Result<T, E>,
) -> Result<T, E> {
    let started = Instant::now();
    println!("[{}] === START: {label} ===", utc_timestamp());
    match operation() {
        Ok(value) => {
            println!(
                "[{}] === FINISHED: {label} ({:.2}s) ===",
                utc_timestamp(),
                started.elapsed().as_secs_f64()
            );
            Ok(value)
        }
        Err(error) => {
            eprintln!(
                "[{}] === FAILED: {label} after {:.2}s: {error} ===",
                utc_timestamp(),
                started.elapsed().as_secs_f64()
            );
            Err(error)
        }
    }
}

fn utc_timestamp() -> String {
    let seconds = match SystemTime::now().duration_since(UNIX_EPOCH) {
        Ok(duration) => duration.as_secs() as i64,
        Err(error) => -(error.duration().as_secs() as i64),
    };
    format_utc_timestamp(seconds)
}

fn format_utc_timestamp(seconds: i64) -> String {
    let days = seconds.div_euclid(86_400);
    let day_seconds = seconds.rem_euclid(86_400);

    // Convert days since the Unix epoch to a Gregorian calendar date.
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 }.div_euclid(146_097);
    let day_of_era = z - era * 146_097;
    let year_of_era =
        (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let mut year = year_of_era + era * 400;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let month_prime = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * month_prime + 2) / 5 + 1;
    let month = month_prime + if month_prime < 10 { 3 } else { -9 };
    year += i64::from(month <= 2);
    let hour = day_seconds / 3_600;
    let minute = day_seconds % 3_600 / 60;
    let second = day_seconds % 60;

    format!("{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}Z")
}

/// Opens a FASTQ file for buffered reading, automatically decompressing if gzipped.
pub fn open_fastq_reader(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let is_gz = match file.read_exact(&mut magic) {
        Ok(()) => magic == [0x1f, 0x8b],
        Err(_) => false,
    };
    file.seek(SeekFrom::Start(0))?;

    if is_gz {
        let gz = MultiGzDecoder::new(BufReader::with_capacity(128 * 1024, file));
        Ok(Box::new(BufReader::with_capacity(128 * 1024, gz)))
    } else {
        Ok(Box::new(BufReader::with_capacity(128 * 1024, file)))
    }
}

pub fn read_fastq_record<R: BufRead + ?Sized>(
    reader: &mut R,
) -> std::io::Result<Option<(String, String, String, String)>> {
    let mut header = String::new();
    let mut sequence = String::new();
    let mut plus = String::new();
    let mut quality = String::new();
    let mut line = String::new();
    if !read_fastq_record_into(
        reader,
        &mut header,
        &mut sequence,
        &mut plus,
        &mut quality,
        &mut line,
    )? {
        return Ok(None);
    }

    Ok(Some((header, sequence, plus, quality)))
}

pub(crate) fn read_fastq_record_into<R: BufRead + ?Sized>(
    reader: &mut R,
    header: &mut String,
    sequence: &mut String,
    plus: &mut String,
    quality: &mut String,
    line: &mut String,
) -> std::io::Result<bool> {
    fn read_line<R: BufRead + ?Sized>(reader: &mut R, line: &mut String) -> std::io::Result<bool> {
        line.clear();
        if reader.read_line(line)? == 0 {
            return Ok(false);
        }
        while line.ends_with(['\n', '\r']) {
            line.pop();
        }
        Ok(true)
    }

    header.clear();
    if !read_line(reader, header)? {
        return Ok(false);
    }
    if !header.starts_with('@') {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("invalid FASTQ header: {header}"),
        ));
    }

    sequence.clear();
    loop {
        if !read_line(reader, line)? {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "FASTQ record is missing its plus line",
            ));
        }
        if line.starts_with('+') {
            plus.clear();
            plus.push_str(line);
            break;
        }
        sequence.push_str(line);
    }

    quality.clear();
    quality.reserve(sequence.len());
    while quality.len() < sequence.len() {
        if !read_line(reader, line)? {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "FASTQ record has fewer quality characters than sequence bases",
            ));
        }
        quality.push_str(line);
    }
    if quality.len() != sequence.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "FASTQ record has more quality characters than sequence bases",
        ));
    }

    Ok(true)
}

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

static SEED: AtomicU64 = AtomicU64::new(0);
static HAS_SEED: AtomicBool = AtomicBool::new(false);

/// Set or clear the global random seed for deterministic execution.
pub fn set_seed(seed: Option<u64>) {
    match seed {
        Some(s) => {
            SEED.store(s, Ordering::Relaxed);
            HAS_SEED.store(true, Ordering::Relaxed);
        }
        None => {
            HAS_SEED.store(false, Ordering::Relaxed);
        }
    }
}

/// Returns whether a seed was configured.
pub fn has_seed() -> bool {
    HAS_SEED.load(Ordering::Relaxed)
}

/// Returns the configured seed, if any.
pub fn get_seed() -> Option<u64> {
    if has_seed() {
        Some(SEED.load(Ordering::Relaxed))
    } else {
        None
    }
}

/// Deterministically shuffle a slice using a simple 64-bit PRNG.
pub fn deterministic_shuffle<T>(vec: &mut [T], seed: u64) {
    let mut rng = seed ^ 0x517cc1b727220a95;
    for i in (1..vec.len()).rev() {
        rng ^= rng >> 12;
        rng ^= rng << 25;
        rng ^= rng >> 27;
        let r = (rng.wrapping_mul(0x2545f4914f6cdd1d) % ((i + 1) as u64)) as usize;
        vec.swap(i, r);
    }
}

/// Returns ordered items if a seed is set, preserving natural iteration order if not.
pub fn order_keys<T: Ord + Clone>(iter: impl IntoIterator<Item = T>) -> Vec<T> {
    let mut vec: Vec<T> = iter.into_iter().collect();
    if let Some(seed) = get_seed() {
        vec.sort();
        if seed != 0 {
            deterministic_shuffle(&mut vec, seed);
        }
    }
    vec
}

/// Returns ordered (key, value) entries by key if a seed is set.
pub fn order_entries_by_key<K: Ord, V>(iter: impl IntoIterator<Item = (K, V)>) -> Vec<(K, V)> {
    let mut vec: Vec<(K, V)> = iter.into_iter().collect();
    if let Some(seed) = get_seed() {
        vec.sort_by(|a, b| a.0.cmp(&b.0));
        if seed != 0 {
            deterministic_shuffle(&mut vec, seed);
        }
    }
    vec
}

/// Get the reverse-complement of a node (flip trailing '+' <-> '-').
pub fn rc_node(id: &str) -> String {
    if let Some(last) = id.chars().last() {
        if last == '+' {
            let base = &id[..id.len() - 1];
            return format!("{}-", base);
        } else if last == '-' {
            let base = &id[..id.len() - 1];
            return format!("{}+", base);
        }
    }
    id.to_string()
}

/// Delete a set of nodes (both orientations) from the graph and remove associated edges
pub fn delete_nodes_and_edges(graph: &mut OverlapGraph, nodes_to_delete: &HashSet<String>) {
    // Initialize set of nodes to remove
    let mut oriented_nodes_to_delete: HashSet<String> = HashSet::new();

    for node in nodes_to_delete.iter() {
        // add to set of nodes to delete
        oriented_nodes_to_delete.insert(node.clone());
        // add rc counterpart
        oriented_nodes_to_delete.insert(rc_node(node));
    }

    // Delete nodes from graph.nodes
    for oriented_node in oriented_nodes_to_delete.iter() {
        graph.nodes.remove(oriented_node);
    }

    // Remove edges pointing to removed nodes
    for (_src, node) in graph.nodes.iter_mut() {
        node.edges
            .retain(|e| !oriented_nodes_to_delete.contains(&e.target_id));
    }
}

/// Check if the bigraph is synchronized:
/// 1. Every node has a reverse complement.
/// 2. Ingoing edges of every node correspond to outgoing edges of its reverse complement.
pub fn check_synchronization(g: &OverlapGraph) {
    for n in g.nodes.keys() {
        // compute reverse complement node
        let n_rc = rc_node(n);

        // check that the reverse complement exists
        assert!(
            g.nodes.contains_key(&n_rc),
            "Reverse complement not found for {}. The bigraph is not synchronized.",
            n
        );

        // check that every outgoing edge has a counterpart in the reverse complement node
        if let Some(node) = g.nodes.get(n) {
            for e in &node.edges {
                let t = &e.target_id;
                let t_rc = rc_node(t);

                // get the reverse complement node for t_rc
                if let Some(t_rc_node) = g.nodes.get(&t_rc) {
                    // Check if there's a matching edge from t_rc to n_rc
                    assert!(
                        t_rc_node.edges.iter().any(|e| e.target_id == n_rc),
                        "Corresponding edge not found for {}. The bigraph is not synchronized.",
                        n
                    );
                } else {
                    panic!(
                        "Reverse complement node {} missing for target {}. The bigraph is not synchronized.",
                        t_rc, t
                    );
                }
            }
        }
    }
}

pub fn rev_comp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'N' | 'n' => 'N',
            _ => c,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{format_utc_timestamp, read_fastq_record, read_fastq_record_into};
    use std::io::Cursor;

    #[test]
    fn reads_wrapped_sequence_and_quality_lines() {
        let input = b"@read1 metadata\nACGT\nTGCA\n+read1\nIIII\n@@@@\n";
        let mut reader = Cursor::new(input);

        let record = read_fastq_record(&mut reader).unwrap().unwrap();

        assert_eq!(record.0, "@read1 metadata");
        assert_eq!(record.1, "ACGTTGCA");
        assert_eq!(record.2, "+read1");
        assert_eq!(record.3, "IIII@@@@");
        assert!(read_fastq_record(&mut reader).unwrap().is_none());
    }

    #[test]
    fn reads_fastq_records_into_reused_buffers() {
        let input = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\n!!!!\n";
        let mut reader = Cursor::new(input);
        let mut header = String::new();
        let mut sequence = String::new();
        let mut plus = String::new();
        let mut quality = String::new();
        let mut line = String::new();

        assert!(
            read_fastq_record_into(
                &mut reader,
                &mut header,
                &mut sequence,
                &mut plus,
                &mut quality,
                &mut line,
            )
            .unwrap()
        );
        assert_eq!(
            (&*header, &*sequence, &*plus, &*quality),
            ("@read1", "ACGT", "+", "IIII")
        );

        assert!(
            read_fastq_record_into(
                &mut reader,
                &mut header,
                &mut sequence,
                &mut plus,
                &mut quality,
                &mut line,
            )
            .unwrap()
        );
        assert_eq!(
            (&*header, &*sequence, &*plus, &*quality),
            ("@read2", "TGCA", "+", "!!!!")
        );
        assert!(
            !read_fastq_record_into(
                &mut reader,
                &mut header,
                &mut sequence,
                &mut plus,
                &mut quality,
                &mut line,
            )
            .unwrap()
        );
    }

    #[test]
    fn formats_utc_timestamps() {
        assert_eq!(format_utc_timestamp(0), "1970-01-01T00:00:00Z");
        assert_eq!(format_utc_timestamp(1_600_000_000), "2020-09-13T12:26:40Z");
    }
}

/// Run minimap2 without forwarding its progress output to the assembler's console.
pub fn run_minimap2(
    args: &[&str],
    threads: usize,
    target: &Path,
    query: &Path,
    out: &Path,
) -> std::io::Result<()> {
    let status = Command::new("minimap2")
        .args(args)
        .arg("-t")
        .arg(threads.max(1).to_string())
        .arg(target)
        .arg(query)
        .arg("-o")
        .arg(out)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()?;
    if status.success() {
        Ok(())
    } else {
        Err(std::io::Error::other(format!(
            "minimap2 failed ({status}) on {}",
            query.display()
        )))
    }
}
