use crate::compress_graph::CompressedGraph;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::process::Command;

pub fn remove_contained_segments(
    assembly: &mut CompressedGraph,
    fasta_path: &Path,
    output_dir: &Path,
    threads: usize,
) -> io::Result<usize> {
    if assembly.unitigs.is_empty() {
        return Ok(0);
    }

    let paf_path = output_dir.join("assembly_containment.paf");
    let status = Command::new("minimap2")
        .arg("-x")
        .arg("asm10")
        .arg("-DP")
        .arg("-t")
        .arg(threads.to_string())
        .arg(fasta_path)
        .arg(fasta_path)
        .arg("-o")
        .arg(&paf_path)
        .status()?;
    if !status.success() {
        return Err(io::Error::new(
            io::ErrorKind::Other,
            "minimap2 failed during assembly segment containment filtering",
        ));
    }

    let contained = contained_unitig_ids(BufReader::new(File::open(paf_path)?))?;
    let before = assembly.unitigs.len();
    assembly
        .unitigs
        .retain(|unitig| !contained.contains(&unitig.id));
    let retained: HashSet<_> = assembly.unitigs.iter().map(|unitig| unitig.id).collect();
    assembly
        .edges
        .retain(|edge| retained.contains(&edge.from) && retained.contains(&edge.to));

    let mut fasta = File::create(fasta_path)?;
    for unitig in &assembly.unitigs {
        let sequence = unitig.fasta_seq.as_deref().unwrap_or("");
        writeln!(
            fasta,
            ">unitig_{} len={}bp topology={}",
            unitig.id,
            sequence.len(),
            unitig.topology
        )?;
        writeln!(fasta, "{sequence}")?;
    }

    Ok(before - assembly.unitigs.len())
}

fn contained_unitig_ids(reader: impl BufRead) -> io::Result<HashSet<usize>> {
    let mut alignments: HashMap<(String, String), QueryAlignments> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 12 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid PAF record during assembly segment containment filtering",
            ));
        }

        let query_length = parse_paf_number(fields[1])?;
        let query_start = parse_paf_number(fields[2])?;
        let query_end = parse_paf_number(fields[3])?;
        let target_length = parse_paf_number(fields[6])?;

        if fields[0] == fields[5]
            || query_length == 0
            || query_start > query_end
            || query_end > query_length
        {
            continue;
        }

        let key = (fields[0].to_string(), fields[5].to_string());
        let alignment = alignments.entry(key).or_insert_with(|| QueryAlignments {
            query_length,
            target_length,
            intervals: Vec::new(),
        });
        alignment.intervals.push((query_start, query_end));
    }

    let mut contained = HashSet::new();
    for ((query_name, target_name), mut alignment) in alignments {
        if alignment.target_length < alignment.query_length
            || (alignment.target_length == alignment.query_length && target_name >= query_name)
            || covered_query_fraction(&mut alignment.intervals, alignment.query_length) <= 0.95
        {
            continue;
        }
        if let Some(id) = query_name
            .strip_prefix("unitig_")
            .and_then(|name| name.parse().ok())
        {
            contained.insert(id);
        }
    }
    Ok(contained)
}

struct QueryAlignments {
    query_length: u64,
    target_length: u64,
    intervals: Vec<(u64, u64)>,
}

fn covered_query_fraction(intervals: &mut [(u64, u64)], query_length: u64) -> f64 {
    intervals.sort_unstable();
    let mut covered = 0;
    let mut current: Option<(u64, u64)> = None;

    for &(start, end) in intervals.iter() {
        match current {
            Some((current_start, current_end)) if start <= current_end => {
                current = Some((current_start, current_end.max(end)));
            }
            Some((current_start, current_end)) => {
                covered += current_end - current_start;
                current = Some((start, end));
            }
            None => current = Some((start, end)),
        }
    }
    if let Some((start, end)) = current {
        covered += end - start;
    }

    covered as f64 / query_length as f64
}

fn parse_paf_number(value: &str) -> io::Result<u64> {
    value.parse().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid numeric field in PAF record",
        )
    })
}

#[cfg(test)]
mod tests {
    use super::contained_unitig_ids;
    use std::io::Cursor;

    #[test]
    fn removes_only_shorter_segments_passing_query_coverage_threshold() {
        let paf = concat!(
            "unitig_1\t100\t0\t96\t+\tunitig_2\t200\t20\t116\t96\t96\t60\n",
            "unitig_3\t100\t0\t95\t+\tunitig_4\t200\t20\t115\t95\t95\t60\n",
            "unitig_5\t100\t0\t99\t+\tunitig_6\t200\t20\t119\t99\t100\t60\n",
            "unitig_7\t200\t0\t198\t+\tunitig_8\t100\t0\t98\t198\t198\t60\n",
            "unitig_9\t100\t0\t99\t+\tunitig_9\t100\t0\t99\t99\t99\t60\n",
            "unitig_10\t100\t0\t99\t+\tunitig_11\t100\t0\t99\t99\t99\t60\n",
            "unitig_11\t100\t0\t99\t+\tunitig_10\t100\t0\t99\t99\t99\t60\n",
        );

        let contained = contained_unitig_ids(Cursor::new(paf)).unwrap();

        assert_eq!(contained, [1, 5, 11].into_iter().collect());
    }

    #[test]
    fn combines_non_overlapping_query_alignments_but_not_different_targets() {
        let paf = concat!(
            "unitig_20\t100\t0\t50\t+\tunitig_21\t200\t0\t50\t50\t50\t60\n",
            "unitig_20\t100\t50\t100\t-\tunitig_21\t200\t50\t100\t50\t50\t60\n",
            "unitig_22\t100\t0\t50\t+\tunitig_23\t200\t0\t50\t50\t50\t60\n",
            "unitig_22\t100\t50\t100\t+\tunitig_24\t200\t0\t50\t50\t50\t60\n",
        );

        let contained = contained_unitig_ids(Cursor::new(paf)).unwrap();

        assert_eq!(contained, [20].into_iter().collect());
    }
}
