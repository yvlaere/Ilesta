use crate::compress_graph::CompressedGraph;
use crate::create_overlap_graph::OverlapGraph;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::Command;

#[derive(Debug, PartialEq)]
pub struct BridgeSupport {
    pub from_node: String,
    pub from_orientation: char,
    pub to_node: String,
    pub to_orientation: char,
    pub support_reads: usize,
    pub edge_len: u32,
    pub overlap_len: u32,
    pub identity: f64,
}

#[derive(Debug, PartialEq)]
pub struct CompletionEdge {
    pub from_node: String,
    pub from_orientation: char,
    pub to_node: String,
    pub to_orientation: char,
    pub support_reads: usize,
    pub edge_len: u32,
    pub overlap_len: u32,
    pub identity: f64,
}

fn run_minimap2_against_reference(
    query_fastq: &Path,
    reference_fasta: &Path,
    output_paf: &Path,
    threads: usize,
    read_type: crate::cli::ReadType,
) -> std::io::Result<()> {
    let mut cmd = Command::new("minimap2");
    for arg in read_type.mapping_args() {
        cmd.arg(arg);
    }
    let status = cmd
        .arg("-N")
        .arg("50")
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
            format!(
                "minimap2 failed for {} against {}",
                query_fastq.display(),
                reference_fasta.display()
            ),
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

    let mut alignments_by_query: HashMap<String, Vec<(String, char, u32, u32, f64)>> =
        HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        let qname = fields[0];
        let qstart: u32 = fields[2].parse().unwrap_or(0);
        let qend: u32 = fields[3].parse().unwrap_or(0);
        let orientation = fields[4].chars().next().unwrap_or('+');
        let tname = fields[5];
        let nmatch: u32 = fields[9].parse().unwrap_or(0);
        let block_len: u32 = fields[10].parse().unwrap_or(0);
        let alignment_len = qend.saturating_sub(qstart);
        if alignment_len < min_alignment_len || block_len == 0 {
            continue;
        }

        let identity = nmatch as f64 / block_len as f64;
        if identity < min_identity {
            continue;
        }

        alignments_by_query
            .entry(qname.to_string())
            .or_default()
            .push((tname.to_string(), orientation, qstart, qend, identity));
    }

    let mut bridge_support: HashMap<(String, char, String, char), BridgeSupport> = HashMap::new();

    for mut alignments in alignments_by_query.into_values() {
        alignments.sort_by(|a, b| {
            a.2.cmp(&b.2)
                .then_with(|| b.3.cmp(&a.3))
                .then_with(|| a.0.cmp(&b.0))
        });
        let mut ordered_hits = Vec::new();
        let mut previous_end = 0;
        for hit in alignments {
            if hit.2 < previous_end {
                continue;
            }
            previous_end = hit.3;
            ordered_hits.push(hit);
        }

        for pair in ordered_hits.windows(2) {
            let (from_node, from_orientation, from_start, from_end, from_identity) = &pair[0];
            let (to_node, to_orientation, to_start, to_end, to_identity) = &pair[1];
            if from_node == to_node || from_end > to_start {
                continue;
            }
            let entry = bridge_support
                .entry((
                    from_node.clone(),
                    *from_orientation,
                    to_node.clone(),
                    *to_orientation,
                ))
                .or_insert(BridgeSupport {
                    from_node: from_node.clone(),
                    from_orientation: *from_orientation,
                    to_node: to_node.clone(),
                    to_orientation: *to_orientation,
                    support_reads: 0,
                    edge_len: 0,
                    overlap_len: 0,
                    identity: 0.0,
                });
            entry.support_reads += 1;
            let candidate_len = from_end
                .saturating_sub(*from_start)
                .min(to_end.saturating_sub(*to_start));
            entry.edge_len = entry.edge_len.max(candidate_len);
            entry.overlap_len = entry.overlap_len.max(candidate_len);
            entry.identity = entry.identity.max((from_identity + to_identity) / 2.0);
        }
    }

    let supports: Vec<_> = bridge_support.into_values().collect();
    let mut best_outgoing: HashMap<(String, char), (usize, usize)> = HashMap::new();
    let mut best_incoming: HashMap<(String, char), (usize, usize)> = HashMap::new();
    for support in &supports {
        update_best(
            &mut best_outgoing,
            (support.from_node.clone(), support.from_orientation),
            support.support_reads,
        );
        update_best(
            &mut best_incoming,
            (support.to_node.clone(), support.to_orientation),
            support.support_reads,
        );
    }

    let mut bridges = Vec::new();
    for support in supports {
        let unique_outgoing = best_outgoing
            .get(&(support.from_node.clone(), support.from_orientation))
            .is_some_and(|(score, count)| *score == support.support_reads && *count == 1);
        let unique_incoming = best_incoming
            .get(&(support.to_node.clone(), support.to_orientation))
            .is_some_and(|(score, count)| *score == support.support_reads && *count == 1);
        if support.support_reads >= 2 && unique_outgoing && unique_incoming {
            bridges.push(support);
        }
    }

    bridges.sort_by(|a, b| {
        b.support_reads
            .cmp(&a.support_reads)
            .then_with(|| b.edge_len.cmp(&a.edge_len))
            .then_with(|| a.from_node.cmp(&b.from_node))
            .then_with(|| a.from_orientation.cmp(&b.from_orientation))
            .then_with(|| a.to_node.cmp(&b.to_node))
            .then_with(|| a.to_orientation.cmp(&b.to_orientation))
    });
    Ok(bridges)
}

fn update_best<K: std::hash::Hash + Eq>(
    best: &mut HashMap<K, (usize, usize)>,
    key: K,
    score: usize,
) {
    match best.get_mut(&key) {
        Some((best_score, count)) if score > *best_score => {
            *best_score = score;
            *count = 1;
        }
        Some((best_score, count)) if score == *best_score => *count += 1,
        None => {
            best.insert(key, (score, 1));
        }
        _ => {}
    }
}

pub fn run_completion_round(
    _graph: &OverlapGraph,
    reads_fastq: &Path,
    unitigs_fasta: &Path,
    output_paf: &Path,
    min_alignment_len: u32,
    min_identity: f64,
    threads: usize,
    read_type: crate::cli::ReadType,
) -> std::io::Result<Vec<CompletionEdge>> {
    run_minimap2_against_reference(reads_fastq, unitigs_fasta, output_paf, threads, read_type)?;
    let bridges = parse_paf_bridges(output_paf, min_alignment_len, min_identity)?;
    let mut edges = Vec::new();

    for bridge in bridges {
        edges.push(CompletionEdge {
            from_node: bridge.from_node,
            from_orientation: bridge.from_orientation,
            to_node: bridge.to_node,
            to_orientation: bridge.to_orientation,
            support_reads: bridge.support_reads,
            edge_len: bridge.edge_len,
            overlap_len: bridge.overlap_len,
            identity: bridge.identity,
        });
    }

    Ok(edges)
}

pub fn apply_completion_edges(
    graph: &mut OverlapGraph,
    compressed: &CompressedGraph,
    edges: &[CompletionEdge],
) -> usize {
    let mut joins = 0;
    for edge in edges {
        let Some(from_id) = edge
            .from_node
            .strip_prefix("unitig_")
            .and_then(|id| id.parse::<usize>().ok())
        else {
            continue;
        };
        let Some(to_id) = edge
            .to_node
            .strip_prefix("unitig_")
            .and_then(|id| id.parse::<usize>().ok())
        else {
            continue;
        };
        let Some(from_unitig) = compressed
            .unitigs
            .iter()
            .find(|unitig| unitig.id == from_id)
        else {
            continue;
        };
        let Some(to_unitig) = compressed.unitigs.iter().find(|unitig| unitig.id == to_id) else {
            continue;
        };
        let (Some(from_first), Some(from_last), Some(to_first), Some(to_last)) = (
            from_unitig.members.first(),
            from_unitig.members.last(),
            to_unitig.members.first(),
            to_unitig.members.last(),
        ) else {
            continue;
        };

        let (from_node, to_node) = match (edge.from_orientation, edge.to_orientation) {
            ('+', '+') => (from_last.node_id.clone(), to_first.node_id.clone()),
            ('+', '-') => (
                from_last.node_id.clone(),
                crate::utils::rc_node(&to_last.node_id),
            ),
            ('-', '+') => (
                crate::utils::rc_node(&from_first.node_id),
                to_first.node_id.clone(),
            ),
            ('-', '-') => (
                crate::utils::rc_node(&from_first.node_id),
                crate::utils::rc_node(&to_last.node_id),
            ),
            _ => continue,
        };
        let from_length = node_interval_length(&from_node);
        let to_length = node_interval_length(&to_node);
        let was_added = graph.add_bridge_edge(&from_node, &to_node, from_length, 0, edge.identity);
        let reverse_added = graph.add_bridge_edge(
            &crate::utils::rc_node(&to_node),
            &crate::utils::rc_node(&from_node),
            to_length,
            0,
            edge.identity,
        );
        if was_added || reverse_added {
            joins += 1;
        }
    }
    joins
}

fn node_interval_length(node_id: &str) -> u32 {
    let Some((_, interval)) = node_id.rsplit_once(':') else {
        return 0;
    };
    let Some((start, end)) = interval[..interval.len().saturating_sub(1)].split_once('-') else {
        return 0;
    };
    end.parse::<u32>()
        .unwrap_or(0)
        .saturating_sub(start.parse::<u32>().unwrap_or(0))
}

#[cfg(test)]
mod tests {
    use super::{BridgeSupport, CompletionEdge, apply_completion_edges, parse_paf_bridges};
    use crate::compress_graph::{CompressedGraph, Unitig, UnitigMember};
    use crate::create_overlap_graph::OverlapGraph;
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn parses_read_bridge_support_between_two_targets() {
        let mut tmp = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        tmp.push("target");
        tmp.push("tmp_completion_test.paf");
        let paf = concat!(
            "read1\t1600\t800\t1550\t-\tunitig_1\t1100\t0\t750\t650\t750\t60\n",
            "read1\t1600\t0\t800\t+\tunitig_0\t1200\t0\t800\t700\t800\t60\n",
            "read2\t1600\t800\t1550\t-\tunitig_1\t1100\t0\t750\t650\t750\t60\n",
            "read2\t1600\t0\t800\t+\tunitig_0\t1200\t0\t800\t700\t800\t60\n"
        );
        fs::write(&tmp, paf).unwrap();

        let bridges = parse_paf_bridges(&tmp, 500, 0.8).unwrap();

        assert_eq!(bridges.len(), 1);
        assert_eq!(
            bridges[0],
            BridgeSupport {
                from_node: "unitig_0".to_string(),
                from_orientation: '+',
                to_node: "unitig_1".to_string(),
                to_orientation: '-',
                support_reads: 2,
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
            "read1\t1800\t0\t900\t+\tunitig_0\t1200\t0\t900\t800\t900\t60\n",
            "read1\t1800\t900\t1750\t+\tunitig_1\t1100\t0\t850\t750\t850\t60\n",
            "read2\t1800\t0\t900\t+\tunitig_0\t1200\t0\t900\t800\t900\t60\n",
            "read2\t1800\t900\t1750\t+\tunitig_1\t1100\t0\t850\t750\t850\t60\n"
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

    #[test]
    fn orders_equally_supported_bridges_by_endpoints() {
        let mut tmp = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        tmp.push("target");
        tmp.push("tmp_completion_test_ties.paf");
        let paf = concat!(
            "read1\t1200\t0\t600\t+\tunitig_b\t800\t0\t600\t580\t600\t60\n",
            "read1\t1200\t600\t1200\t+\tunitig_c\t800\t0\t600\t580\t600\t60\n",
            "read2\t1200\t0\t600\t+\tunitig_b\t800\t0\t600\t580\t600\t60\n",
            "read2\t1200\t600\t1200\t+\tunitig_c\t800\t0\t600\t580\t600\t60\n",
            "read3\t1200\t0\t600\t+\tunitig_a\t800\t0\t600\t580\t600\t60\n",
            "read3\t1200\t600\t1200\t+\tunitig_d\t800\t0\t600\t580\t600\t60\n",
            "read4\t1200\t0\t600\t+\tunitig_a\t800\t0\t600\t580\t600\t60\n",
            "read4\t1200\t600\t1200\t+\tunitig_d\t800\t0\t600\t580\t600\t60\n"
        );
        fs::write(&tmp, paf).unwrap();

        let bridges = parse_paf_bridges(&tmp, 500, 0.8).unwrap();

        assert_eq!(bridges.len(), 2);
        assert_eq!(bridges[0].from_node, "unitig_a");
        assert_eq!(bridges[1].from_node, "unitig_b");
        let _ = fs::remove_file(tmp);
    }

    #[test]
    fn drops_tied_outgoing_joins_as_ambiguous() {
        let mut tmp = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        tmp.push("target");
        tmp.push("tmp_completion_test_ambiguous.paf");
        let paf = concat!(
            "read1\t1200\t0\t600\t+\tunitig_0\t800\t0\t600\t580\t600\t60\n",
            "read1\t1200\t600\t1200\t+\tunitig_1\t800\t0\t600\t580\t600\t60\n",
            "read2\t1200\t0\t600\t+\tunitig_0\t800\t0\t600\t580\t600\t60\n",
            "read2\t1200\t600\t1200\t+\tunitig_1\t800\t0\t600\t580\t600\t60\n",
            "read3\t1200\t0\t600\t+\tunitig_0\t800\t0\t600\t580\t600\t60\n",
            "read3\t1200\t600\t1200\t+\tunitig_2\t800\t0\t600\t580\t600\t60\n",
            "read4\t1200\t0\t600\t+\tunitig_0\t800\t0\t600\t580\t600\t60\n",
            "read4\t1200\t600\t1200\t+\tunitig_2\t800\t0\t600\t580\t600\t60\n"
        );
        fs::write(&tmp, paf).unwrap();

        let bridges = parse_paf_bridges(&tmp, 500, 0.8).unwrap();
        assert!(bridges.is_empty());

        let _ = fs::remove_file(tmp);
    }

    #[test]
    fn applies_oriented_unitig_join_to_graph_endpoints() {
        let compressed = CompressedGraph {
            unitigs: vec![
                Unitig {
                    id: 0,
                    members: vec![UnitigMember {
                        node_id: "read_a:0-10+".to_string(),
                        edge: (String::new(), 0),
                    }],
                    fasta_seq: Some("AAAAAAAAAA".to_string()),
                    topology: 'l',
                },
                Unitig {
                    id: 1,
                    members: vec![UnitigMember {
                        node_id: "read_b:3-13+".to_string(),
                        edge: (String::new(), 0),
                    }],
                    fasta_seq: Some("CCCCCCCCCC".to_string()),
                    topology: 'l',
                },
            ],
            edges: Vec::new(),
        };
        let edge = CompletionEdge {
            from_node: "unitig_0".to_string(),
            from_orientation: '+',
            to_node: "unitig_1".to_string(),
            to_orientation: '-',
            support_reads: 2,
            edge_len: 10,
            overlap_len: 0,
            identity: 0.99,
        };
        let mut graph = OverlapGraph::new();

        assert_eq!(apply_completion_edges(&mut graph, &compressed, &[edge]), 1);
        assert_eq!(
            graph.nodes["read_a:0-10+"].edges[0].target_id,
            "read_b:3-13-"
        );
        assert_eq!(
            graph.nodes["read_b:3-13+"].edges[0].target_id,
            "read_a:0-10-"
        );
    }
}
