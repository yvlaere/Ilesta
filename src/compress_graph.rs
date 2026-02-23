use crate::utils;
use crate::create_overlap_graph::OverlapGraph;
use crate::alignment_filtering::Overlap;
/// graph compression module
/// creates a compressed graph of unitigs from an overlap graph
/// 1. get the indegree and outdegree of each node
/// 2. get non-circular unitigs (start at nodes with indegree != 1 or outdegree != 1)
/// 3. get circular unitigs (remaining unvisited nodes)
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

pub struct UnitigMember {
    pub node_id: String,
    // id of target node and edge length to that node
    pub edge: (String, u32),
}

pub struct Unitig {
    pub id: usize,
    pub members: Vec<UnitigMember>,
    pub fasta_seq: Option<String>,
}

pub struct UnitigEdge {
    pub from: usize,
    pub to: usize,
    pub edge_len: u32,
    pub overlap_len: u32,
    pub identity: f64,
}

pub struct CompressedGraph {
    pub unitigs: Vec<Unitig>,
    pub edges: Vec<UnitigEdge>,
}

impl CompressedGraph {
    pub fn write_gfa(&self, path: &str, overlaps: &HashMap<(usize, usize), Overlap>) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = std::fs::File::create(path)?;
        use std::io::Write;
        // header
        writeln!(file, "H\tVN:Z:1.0")?;

        // segments
        for u in &self.unitigs {
            let sid = format!("unitig_{}", u.id);
            let seq = u.fasta_seq.as_ref().map(|s| s.as_str()).unwrap_or("*");
            writeln!(file, "S\t{}\t{}", sid, seq)?;
        }

        // links (including circularising links)
        let mut edge_set = std::collections::HashSet::new();
        for e in &self.edges {
            let from = format!("unitig_{}", e.from);
            let to = format!("unitig_{}", e.to);
            let cigar = format!("{}M", e.overlap_len);
            writeln!(file, "L\t{}\t+\t{}\t+\t{}", from, to, cigar)?;
            edge_set.insert((e.from, e.to));
        }

        // Add circularising links for unitigs that are circular (first member's node_id == last member's edge.0)
        for u in &self.unitigs {
            if u.members.len() > 1 {
                let first = &u.members[0].node_id;
                let last_edge = &u.members[u.members.len() - 1].edge;
                let last = &last_edge.0;
                // If last edge points back to first node, and not already in edge_set
                if !last.is_empty() && last == first {
                    let from = format!("unitig_{}", u.id);
                    let to = format!("unitig_{}", u.id);
                    let cigar = format!("{}M", last_edge.1);
                    if !edge_set.contains(&(u.id, u.id)) {
                        writeln!(file, "L\t{}\t+\t{}\t+\t{}", from, to, cigar)?;
                    }
                }

                // if there is an overlap between the first and last unitig members in overlaps, add a link between them
                for ((_query_id, _target_id), o) in overlaps.iter() {
                    if (o.source_name == *first && o.sink_name == *last)  || 
                    (o.source_name == *last && o.sink_name == *first) || 
                    (o.source_name == utils::rc_node(last) && o.sink_name == utils::rc_node(first)) ||
                    (o.source_name == utils::rc_node(first) && o.sink_name == utils::rc_node(last)) {
                        let from = format!("unitig_{}", u.id);
                        let to = format!("unitig_{}", u.id);
                        let cigar = format!("{}M", o.overlap_len);
                        writeln!(file, "L\t{}\t+\t{}\t+\t{}", from, to, cigar)?;
                    }
                }
            }
        }

        // 'a' lines: annotate which reads were used to create each unitig
        for u in &self.unitigs {
            let sid = format!("unitig_{}", u.id);
            let mut reads = std::collections::HashSet::new();
            for m in &u.members {
                // node_id is like "read_id+" or "read_id-"
                let node_id = &m.node_id;
                if node_id.len() > 1 {
                    reads.insert(node_id[..node_id.len()-1].to_string());
                }
            }
            let mut reads_vec: Vec<_> = reads.into_iter().collect();
            reads_vec.sort();
            let reads_str = reads_vec.join(",");
            writeln!(file, "a\t{}\treads:{}", sid, reads_str)?;
        }

        Ok(())
    }
}

/// Main function: compress maximal non-branching paths into unitigs.
/// Preserves member lists and the overlap lengths between them.
pub fn compress_unitigs(graph: &OverlapGraph, fastq_path: &str, fasta_path: &str,
) -> CompressedGraph {
    // 1) create a map of indegrees
    let mut indegree: HashMap<String, usize> = HashMap::new();
    for id in graph.nodes.keys() {
        indegree.insert(id.clone(), 0);
    }
    for (_source_id, node) in &graph.nodes {
        for e in &node.edges {
            *indegree.entry(e.target_id.clone()).or_default() += 1;
        }
    }

    let mut visited: HashSet<String> = HashSet::new();
    let mut unitigs: Vec<Unitig> = Vec::new();

    // Helper to extract the single outgoing neighbor if outdeg == 1
    let out_single =
        |g: &OverlapGraph, cur: &str| -> Option<(String, u32)> {
            g.nodes.get(cur).and_then(|n| {
                if n.edges.len() == 1 {
                    let e = &n.edges[0];
                    Some((e.target_id.clone(), e.edge_len))
                } else {
                    None
                }
            })
        };

    // 2) non-circular unitigs, start unitigs at nodes where indegree != 1 || outdeg != 1
    for (id, node) in &graph.nodes {
        let indegree_i = *indegree.get(id).unwrap_or(&0);
        let outdeg_i = node.edges.len();

        // skip if already visited
        if visited.contains(id) {
            continue;
        }

        // start a unitig if indegree != 1 or outdegree != 1 (i.e., not a simple internal node)
        if indegree_i != 1 || outdeg_i != 1 {
            // create a unitig for each outgoing edge
            // if outdeg == 0, there will be no unitig, each outgoing edge will be zero-indexed
            for out_edge_i in 0..outdeg_i {
                // start a new unitig from id
                let mut members: Vec<UnitigMember> = Vec::new();
                let mut cur = id.clone();
                //members.push(UnitigMember { node_id: cur.clone(), overlap_from_prev: 0 });
                visited.insert(cur.clone());

                // check the next outgoing edge
                let (second, edge_len) = {
                    let e = &node.edges[out_edge_i];
                    (e.target_id.clone(), e.edge_len)
                };
                // push the first node into the unitig members
                members.push(UnitigMember {
                    node_id: cur.clone(),
                    edge: (second.clone(), edge_len),
                });

                let second_indegree = *indegree.get(&second).unwrap_or(&0);
                // check if the node breaks the chain
                if second_indegree != 1 {
                    continue;
                }
                // stop if second is already visited
                if visited.contains(&second) {
                    continue;
                }

                visited.insert(second.clone());
                cur = second.clone();

                // extend forward from second untill the end
                while let Some((next, edge_len)) = out_single(graph, &cur) {
                    let next_indegree = *indegree.get(&next).unwrap_or(&0);
                    // don't add the node that breaks the chain
                    if next_indegree != 1 {
                        break;
                    }
                    // stop if next is already visited
                    if visited.contains(&next) {
                        break;
                    }
                    // push cur to the unitig members
                    members.push(UnitigMember {
                        node_id: cur.clone(),
                        edge: (next.clone(), edge_len),
                    });
                    visited.insert(next.clone());
                    cur = next;
                }

                // add the final node
                members.push(UnitigMember {
                    node_id: cur.clone(),
                    edge: (String::new(), 0),
                });

                // create the unitig
                let uid = unitigs.len();
                unitigs.push(Unitig {
                    id: uid,
                    members,
                    fasta_seq: None,
                });
            }
        }
    }

    // 3) circular unitigs, handle remaining nodes that are still unvisited
    for id in graph.nodes.keys() {
        if visited.contains(id) {
            continue;
        }

        // start a circular unitig
        let mut cur = id.clone();
        let mut members: Vec<UnitigMember> = Vec::new();

        loop {
            visited.insert(cur.clone());
            // to follow, get the unique outgoing edge
            let next_edge = out_single(graph, &cur);
            if next_edge.is_none() {
                // shouldn't happen in pure cycle, break defensively
                break;
            }
            let (next, edge_len) = next_edge.unwrap();

            // push cur to the unitig members
            members.push(UnitigMember {
                node_id: cur.clone(),
                edge: (next.clone(), edge_len),
            });

            let next_indegree = *indegree.get(&next).unwrap_or(&0);

            // stop if the next node breaks the cycle structure (indegree != 1)
            if next_indegree != 1 {
                break;
            }
            // stop if next is already visited
            if visited.contains(&next) {
                break;
            }

            cur = next;
        }

        // register unitig
        let uid = unitigs.len();
        unitigs.push(Unitig {
            id: uid,
            members,
            fasta_seq: None,
        });
    }

    // Build edges between unitigs based on original overlap graph
    let mut node_to_unitig: HashMap<String, usize> = HashMap::new();
    for u in &unitigs {
        for m in &u.members {
            node_to_unitig.insert(m.node_id.clone(), u.id);
        }
    }

    let mut unitig_edge_map: HashMap<(usize, usize), UnitigEdge> = HashMap::new();
    for (source_id, node) in &graph.nodes {
        if let Some(&from_uid) = node_to_unitig.get(source_id) {
            for e in &node.edges {
                if let Some(&to_uid) = node_to_unitig.get(&e.target_id) {
                    if from_uid == to_uid {
                        continue;
                    }
                    let key = (from_uid, to_uid);
                    let entry = unitig_edge_map.entry(key).or_insert(UnitigEdge {
                        from: from_uid,
                        to: to_uid,
                        edge_len: e.edge_len,
                        overlap_len: e.overlap_len,
                        identity: e.identity,
                    });
                    // choose smallest edge_len and best identity
                    if e.edge_len < entry.edge_len {
                        entry.edge_len = e.edge_len;
                    }
                    if e.identity > entry.identity {
                        entry.identity = e.identity;
                    }
                    if e.overlap_len > entry.overlap_len {
                        entry.overlap_len = e.overlap_len;
                    }
                }
            }
        }
    }

    let edges: Vec<UnitigEdge> = unitig_edge_map.into_values().collect();

    // load fastq sequences
    println!("Loading FASTQ sequences from {}...", fastq_path);
    let fastq_seqs = load_fastq_sequences(fastq_path).unwrap();

    // generate fasta sequences for unitigs and write to fasta_path
    println!("Generating unitig sequences and writing to FASTA...");
    for unitig in unitigs.iter_mut() {
        let seq = unitig_sequence(unitig, graph, fastq_seqs.clone()).unwrap();
        unitig.fasta_seq = Some(seq);
    }
    // write to fasta file
    {
        let mut fasta_file = std::fs::File::create(fasta_path).unwrap();
        for unitig in unitigs.iter() {
            let seq = unitig.fasta_seq.as_ref().unwrap();
            let header = format!(">unitig_{} len={}bp\n", unitig.id, seq.len());
            use std::io::Write;
            fasta_file.write_all(header.as_bytes()).unwrap();
            fasta_file.write_all(seq.as_bytes()).unwrap();
            fasta_file.write_all(b"\n").unwrap();
        }
    }

    CompressedGraph { unitigs, edges }
}

fn load_fastq_sequences(fastq_path: &str) -> Result<HashMap<String, String>, String> {
    let mut seq_map: HashMap<String, String> = HashMap::new();

    let reader = match std::fs::File::open(fastq_path) {
        Ok(f) => f,
        Err(e) => return Err(format!("failed to open FASTQ file '{}': {}", fastq_path, e)),
    };
    let buf_reader = std::io::BufReader::new(reader);
    let mut lines = buf_reader.lines();

    while let Some(Ok(header)) = lines.next() {
        if !header.starts_with('@') {
            return Err(format!(
                "invalid FASTQ format: expected header line starting with '@', got '{}'",
                header
            ));
        }
        let seq = match lines.next() {
            Some(Ok(s)) => s,
            _ => return Err("invalid FASTQ format: missing sequence line".to_string()),
        };
        // skip plus line
        lines.next();
        // skip quality line
        lines.next();

        let id = header[1..].split_whitespace().next().unwrap().to_string();
        seq_map.insert(id, seq);
    }

    Ok(seq_map)
}

/// Build nucleotide sequence for `unitig` by concatenating node sequences and
/// removing overlaps recorded in UnitigMember.edge.(target, edge_len).
pub fn unitig_sequence(
    unitig: &Unitig,
    graph: &OverlapGraph,
    fastq_seqs: HashMap<String, String>,
) -> Result<String, String> {
    if unitig.members.is_empty() {
        println!("unitig has no members; cannot infer sequence");
    }

    // Helper to get sequence for a node id
    let get_seq = |node_id: &str| -> Result<String, String> {
        // Node ID must end with '+' or '-'
        let orientation = node_id
            .chars()
            .last()
            .ok_or_else(|| format!("node id '{}' has no orientation suffix (+/-)", node_id))?;
        let read_id = &node_id[..node_id.len() - 1];
        let _node = graph
            .nodes
            .get(node_id)
            .ok_or_else(|| format!("node_id '{}' not found in overlap graph", node_id))?;
        let seq = fastq_seqs.get(read_id).ok_or_else(|| {
            format!(
                "sequence for read_id '{}' not found in FASTQ sequences",
                read_id
            )
        })?;
        match orientation {
            '+' => Ok(seq.clone()),
            '-' => Ok(utils::rev_comp(seq)),
            _ => Err(format!(
                "invalid orientation '{}' in node id '{}'",
                orientation, node_id
            )),
        }
    };

    let mut out = String::new();

    // For each member, append the prefix untill the target (skip overlap bases)
    for member in &unitig.members {
        // an edge describes the target node and the edge length of the original node to reach that target
        let (target_id, edge_len) = &member.edge;
        let seq = get_seq(&member.node_id)?;
        let edge_len_usize = *edge_len as usize;

        if target_id.is_empty() {
            // This is the last node, append the entire sequence
            out.push_str(&seq);
        } else {
            if edge_len_usize > seq.len() {
                return Err(format!(
                    "edge length ({}) larger than seq length ({}) for node {} to target {}",
                    edge_len_usize,
                    seq.len(),
                    member.node_id,
                    target_id
                ));
            }

            // Append the non-overlapping prefix
            out.push_str(&seq[..edge_len_usize]);
        }
    }

    Ok(out)
}
