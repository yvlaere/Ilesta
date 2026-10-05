use crate::compress_graph::CompressedGraph;
use std::collections::HashSet;
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
    let mut contained = HashSet::new();
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
        let matching_bases = parse_paf_number(fields[9])?;
        let alignment_block_length = parse_paf_number(fields[10])?;

        if fields[0] == fields[5]
            || target_length < query_length
            || (target_length == query_length && fields[5] >= fields[0])
            || query_length == 0
            || alignment_block_length == 0
        {
            continue;
        }

        let query_coverage = (query_end - query_start) as f64 / query_length as f64;
        let identity = matching_bases as f64 / alignment_block_length as f64;
        if query_coverage > 0.95 && identity > 0.99 {
            if let Some(id) = fields[0]
                .strip_prefix("unitig_")
                .and_then(|name| name.parse().ok())
            {
                contained.insert(id);
            }
        }
    }
    Ok(contained)
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
    fn removes_only_shorter_segments_passing_both_strict_thresholds() {
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

        assert_eq!(contained, [1, 11].into_iter().collect());
    }
}
