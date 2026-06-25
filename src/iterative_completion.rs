use crate::create_overlap_graph::OverlapGraph;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::Command;

#[derive(Debug, PartialEq)]
pub struct BridgeSupport {
    pub from_node: String,
    pub to_node: String,
    pub support_reads: usize,
    pub edge_len: u32,
    pub overlap_len: u32,
    pub identity: f64,
}

#[derive(Debug, PartialEq)]
pub struct CompletionEdge {
    pub from_node: String,
    pub to_node: String,
    pub edge_len: u32,
    pub overlap_len: u32,
    pub identity: f64,
}

fn run_minimap2_against_reference(
    query_fastq: &Path,
    reference_fasta: &Path,
    output_paf: &Path,
    threads: usize,
) -> std::io::Result<()> {
    let status = Command::new("minimap2")
        .arg("-x")
        .arg("ava-ont")
        .arg("-t")
        .arg(threads.to_string())
        .arg(reference_fasta)
        .arg(query_fastq)
        .arg("-o")
        .arg(output_paf)
        .status()?;

    if !status.success() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("minimap2 failed for {} against {}", query_fastq.display(), reference_fasta.display()),
        ));
    }

    Ok(())
}

fn parse_paf_bridges(
    paf_path: &Path,
    min_alignment_len: u32,
    min_identity: f64,
) -> std::io::Result<Vec<BridgeSupport>> {
    let file = File::open(paf_path)?;
    let reader = BufReader::new(file);

    let mut alignments_by_query: HashMap<String, Vec<(String, u32, f64)>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue;
        }

        let qname = fields[0];
        let qstart: u32 = fields[2].parse().unwrap_or(0);
        let qend: u32 = fields[3].parse().unwrap_or(0);
        let tname = fields[5];
        let nmatch: u32 = fields[9].parse().unwrap_or(0);
        let alignment_len = qend.saturating_sub(qstart);
        if alignment_len < min_alignment_len {
            continue;
        }

        let identity = if alignment_len > 0 {
            nmatch as f64 / alignment_len as f64
        } else {
            0.0
        };
        if identity < min_identity {
            continue;
        }

        alignments_by_query
            .entry(qname.to_string())
            .or_default()
            .push((tname.to_string(), alignment_len, identity));
    }

    let mut bridge_support: HashMap<(String, String), BridgeSupport> = HashMap::new();

    for targets in alignments_by_query.into_values() {
        let mut pairs = Vec::new();
        let mut seen_pairs_for_query: HashSet<(String, String)> = HashSet::new();
        for (target, alignment_len, identity) in &targets {
            for (other_target, other_len, other_identity) in &targets {
                if other_target == target {
                    continue;
                }
                let mut pair = (target.clone(), other_target.clone());
                if pair.0 > pair.1 {
                    std::mem::swap(&mut pair.0, &mut pair.1);
                }
                if seen_pairs_for_query.insert(pair.clone()) {
                    pairs.push((pair, *alignment_len, *other_len, *identity, *other_identity));
                }
            }
        }

        for ((from_node, to_node), from_len, to_len, from_id, to_id) in pairs {
            let entry = bridge_support.entry((from_node.clone(), to_node.clone())).or_insert(BridgeSupport {
                from_node: from_node.clone(),
                to_node: to_node.clone(),
                support_reads: 0,
                edge_len: 0,
                overlap_len: 0,
                identity: 0.0,
            });
            entry.support_reads += 1;
            let candidate_len = from_len.min(to_len);
            entry.edge_len = entry.edge_len.max(candidate_len);
            entry.overlap_len = entry.overlap_len.max(candidate_len);
            entry.identity = entry.identity.max((from_id + to_id) / 2.0);
        }
    }

    let mut bridges = Vec::new();
    for support in bridge_support.into_values() {
        let strong_support = support.support_reads >= 2
            || support.edge_len >= min_alignment_len.saturating_mul(2)
            || (support.edge_len >= min_alignment_len && support.identity >= min_identity + 0.05);

        if strong_support {
            bridges.push(support);
        }
    }

    bridges.sort_by(|a, b| b.support_reads.cmp(&a.support_reads).then_with(|| b.edge_len.cmp(&a.edge_len)));
    Ok(bridges)
}

pub fn run_completion_round(
    graph: &OverlapGraph,
    reads_fastq: &Path,
    unitigs_fasta: &Path,
    output_paf: &Path,
    min_alignment_len: u32,
    min_identity: f64,
    threads: usize,
) -> std::io::Result<Vec<CompletionEdge>> {
    run_minimap2_against_reference(reads_fastq, unitigs_fasta, output_paf, threads)?;
    let bridges = parse_paf_bridges(output_paf, min_alignment_len, min_identity)?;
    let mut edges = Vec::new();

    for bridge in bridges {
        let from_exists = graph.nodes.contains_key(&bridge.from_node);
        let to_exists = graph.nodes.contains_key(&bridge.to_node);
        if !from_exists || !to_exists {
            continue;
        }
        edges.push(CompletionEdge {
            from_node: bridge.from_node,
            to_node: bridge.to_node,
            edge_len: bridge.edge_len,
            overlap_len: bridge.overlap_len,
            identity: bridge.identity,
        });
    }

    Ok(edges)
}

#[cfg(test)]
mod tests {
    use super::{parse_paf_bridges, BridgeSupport};
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn parses_read_bridge_support_between_two_targets() {
        let mut tmp = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        tmp.push("target");
        tmp.push("tmp_completion_test.paf");
        let paf = concat!(
            "read1\t1000\t0\t800\t+\tunitig_0\t1200\t0\t800\t700\t800\t60\n",
            "read1\t1000\t0\t750\t+\tunitig_1\t1100\t0\t750\t650\t750\t60\n"
        );
        fs::write(&tmp, paf).unwrap();

        let bridges = parse_paf_bridges(&tmp, 500, 0.8).unwrap();

        assert_eq!(bridges.len(), 1);
        assert_eq!(
            bridges[0],
            BridgeSupport {
                from_node: "unitig_0".to_string(),
                to_node: "unitig_1".to_string(),
                support_reads: 1,
                edge_len: 750,
                overlap_len: 750,
                identity: 0.8708333333333333,
            }
        );

        let _ = fs::remove_file(tmp);
    }

    #[test]
    fn aggregates_multiple_supporting_reads_for_the_same_bridge() {
        let mut tmp = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        tmp.push("target");
        tmp.push("tmp_completion_test_aggregated.paf");
        let paf = concat!(
            "read1\t1000\t0\t800\t+\tunitig_0\t1200\t0\t800\t700\t800\t60\n",
            "read1\t1000\t0\t750\t+\tunitig_1\t1100\t0\t750\t650\t750\t60\n",
            "read2\t1000\t0\t900\t+\tunitig_0\t1200\t0\t900\t800\t900\t60\n",
            "read2\t1000\t0\t850\t+\tunitig_1\t1100\t0\t850\t750\t850\t60\n"
        );
        fs::write(&tmp, paf).unwrap();

        let bridges = parse_paf_bridges(&tmp, 500, 0.8).unwrap();

        assert_eq!(bridges.len(), 1);
        assert_eq!(bridges[0].support_reads, 2);
        assert_eq!(bridges[0].from_node, "unitig_0");
        assert_eq!(bridges[0].to_node, "unitig_1");
        assert_eq!(bridges[0].edge_len, 850);

        let _ = fs::remove_file(tmp);
    }
}
