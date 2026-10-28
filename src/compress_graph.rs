use crate::alignment_filtering::Overlap;
use crate::create_overlap_graph::OverlapGraph;
use crate::utils;
/// graph compression module
/// creates a compressed graph of unitigs from an overlap graph
/// 1. get the indegree and outdegree of each node
/// 2. get non-circular unitigs (start at nodes with indegree != 1 or outdegree != 1)
/// 3. get circular unitigs (remaining unvisited nodes)
use std::collections::{HashMap, HashSet};

pub struct UnitigMember {
    pub node_id: String,
    // id of target node and edge length to that node
    pub edge: (String, u32),
}

pub struct Unitig {
    pub id: usize,
    pub members: Vec<UnitigMember>,
    pub fasta_seq: Option<String>,
    pub topology: char,
}

pub struct UnitigEdge {
    pub from: usize,
    pub to: usize,
    pub from_ori: char,
    pub to_ori: char,
    pub edge_len: u32,
    pub overlap_len: u32,
    pub identity: f64,
}

pub struct CompressedGraph {
    pub unitigs: Vec<Unitig>,
    pub edges: Vec<UnitigEdge>,
}

impl CompressedGraph {
    /// Write all unitig sequences as FASTA.
    pub fn write_fasta(&self, path: &std::path::Path) -> std::io::Result<()> {
        use std::io::Write;
        let mut out = std::io::BufWriter::new(std::fs::File::create(path)?);
        for u in &self.unitigs {
            let seq = u.fasta_seq.as_deref().unwrap_or("");
            writeln!(out, ">unitig_{} len={}bp topology={}\n{seq}", u.id, seq.len(), u.topology)?;
        }
        out.flush()
    }

    /// Read contigs from a FASTA file (no edges); `topology=c` in a header marks a circular contig.
    pub fn from_fasta(path: &std::path::Path) -> std::io::Result<CompressedGraph> {
        let mut unitigs: Vec<Unitig> = Vec::new();
        for line in std::io::BufRead::lines(std::io::BufReader::new(std::fs::File::open(path)?)) {
            let line = line?;
            if let Some(header) = line.strip_prefix('>') {
                let topology = if header.contains("topology=c") { 'c' } else { 'l' };
                unitigs.push(Unitig { id: unitigs.len(), members: Vec::new(), fasta_seq: Some(String::new()), topology });
            } else if let Some(seq) = unitigs.last_mut().and_then(|u| u.fasta_seq.as_mut()) {
                seq.push_str(line.trim());
            }
        }
        Ok(CompressedGraph { unitigs, edges: Vec::new() })
    }

    pub fn write_gfa(&self, path: &str, read_lines: bool) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = std::fs::File::create(path)?;
        use std::io::Write;
        // header
        writeln!(file, "H\tVN:Z:1.0")?;

        // segments
        for u in &self.unitigs {
            let sid = format!("unitig_{}_{}", u.id, u.topology);
            let seq = u.fasta_seq.as_deref().unwrap_or("*");
            writeln!(file, "S\t{}\t{}", sid, seq)?;
        }

        // Add circularising links for unitigs that are circular (first member's node_id == last member's edge.0)
        for u in &self.unitigs {
            if u.topology == 'c' {
                let from = format!("unitig_{}_{}", u.id, u.topology);
                let to = from.clone();
                let cigar = format!("{}M", 0);
                writeln!(file, "L\t{}\t+\t{}\t+\t{}", from, to, cigar)?;
                writeln!(file, "L\t{}\t-\t{}\t-\t{}", from, to, cigar)?;
            }
        }

        // other links
        let mut edge_set = std::collections::HashSet::new();
        for e in &self.edges {
            // get from unitig
            let from_utg = self
                .unitigs
                .iter()
                .find(|u| u.id == e.from)
                .ok_or("invalid edge: from unitig not found")?;
            let to_utg = self
                .unitigs
                .iter()
                .find(|u| u.id == e.to)
                .ok_or("invalid edge: to unitig not found")?;
            let from = format!("unitig_{}_{}", e.from, from_utg.topology);
            let to = format!("unitig_{}_{}", e.to, to_utg.topology);
            let cigar = format!("{}M", e.overlap_len);
            writeln!(
                file,
                "L\t{}\t{}\t{}\t{}\t{}",
                from, e.from_ori, to, e.to_ori, cigar
            )?;
            edge_set.insert((e.from, e.to));
        }

        // 'a' lines: annotate which reads were used to create each unitig
        for u in self.unitigs.iter().filter(|_| read_lines) {
            let utg_name = format!("unitig_{}_{}", u.id, u.topology);
            let mut utg_pos: u32 = 0;

            for m in &u.members {
                // Each member contributes the part of its read up to the next member, like
                // `unitig_sequence`; only the last member of a linear unitig (no target) contributes
                // its whole slice. The last member of a circular unitig points back to the first one.
                let contrib_len = if m.edge.0.is_empty() {
                    let (_read_name, read_start, read_end, _ori) =
                        parse_node_id(&m.node_id).map_err(|e| format!("a-line error: {}", e))?;
                    read_end - read_start
                } else {
                    m.edge.1
                };

                let (read_name, read_start, read_end, ori) =
                    parse_node_id(&m.node_id).map_err(|e| format!("a-line error: {}", e))?;

                let utg_start = utg_pos;
                let utg_end = utg_pos + contrib_len;

                writeln!(
                    file,
                    "a\t{}\t{}\t{}:{}-{}\t{}\t{}",
                    utg_name, utg_start, read_name, read_start, read_end, ori, utg_end
                )?;

                utg_pos = utg_end;
            }
        }

        Ok(())
    }
}

fn set_unitig_topologies(unitigs: &mut [Unitig], overlaps: &HashMap<(usize, usize), Overlap>) {
    for unitig in unitigs {
        if unitig.members.len() <= 1 {
            continue;
        }

        let first = &unitig.members[0].node_id;
        let last = &unitig.members[unitig.members.len() - 1].edge.0;
        if !last.is_empty() && last == first {
            unitig.topology = 'c';
        }

        for overlap in overlaps.values() {
            if (overlap.source_name == *first && overlap.sink_name == *last)
                || (overlap.source_name == *last && overlap.sink_name == *first)
                || (overlap.rc_source_name == *first && overlap.rc_sink_name == *last)
                || (overlap.rc_source_name == *last && overlap.rc_sink_name == *first)
            {
                unitig.topology = 'c';
                break;
            }
        }
    }
}

// helper function to parse the node_id format "read_name:start-end+/-" and extract the read_name and orientation
fn parse_node_id(node_id: &str) -> Result<(&str, u32, u32, char), String> {
    let ori = node_id
        .chars()
        .last()
        .ok_or_else(|| format!("node_id '{}' missing orientation", node_id))?;

    if ori != '+' && ori != '-' {
        return Err(format!("invalid orientation '{}' in '{}'", ori, node_id));
    }

    let core = &node_id[..node_id.len() - 1];

    let (read_name, range) = core
        .split_once(':')
        .ok_or_else(|| format!("node_id '{}' missing ':'", node_id))?;

    let (start, end) = range
        .split_once('-')
        .ok_or_else(|| format!("node_id '{}' missing '-'", node_id))?;

    let start: u32 = start
        .parse()
        .map_err(|_| format!("invalid start '{}' in '{}'", start, node_id))?;

    let end: u32 = end
        .parse()
        .map_err(|_| format!("invalid end '{}' in '{}'", end, node_id))?;

    if start >= end {
        return Err(format!(
            "invalid interval {}-{} in '{}'",
            start, end, node_id
        ));
    }

    Ok((read_name, start, end, ori))
}

/// Main function: compress maximal non-branching paths into unitigs.
/// Preserves member lists and the overlap lengths between them.
pub fn compress_unitigs(
    graph: &OverlapGraph,
    overlaps: &HashMap<(usize, usize), Overlap>,
    fastq_path: &std::path::Path,
) -> CompressedGraph {
    // 1) create a map of indegrees
    let mut indegree: HashMap<String, usize> = HashMap::new();
    for id in graph.nodes.keys() {
        indegree.insert(id.clone(), 0);
    }
    for node in graph.nodes.values() {
        for e in &node.edges {
            *indegree.entry(e.target_id.clone()).or_default() += 1;
        }
    }

    let mut visited: HashSet<String> = HashSet::new();
    let mut unitigs: Vec<Unitig> = Vec::new();

    // Helper to extract the single outgoing neighbor if outdeg == 1
    let out_single = |g: &OverlapGraph, cur: &str| -> Option<(String, u32)> {
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
    let node_ids = crate::utils::order_keys(graph.nodes.keys().cloned());
    for id in &node_ids {
        let node = &graph.nodes[id];
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
                    topology: 'l',
                });
            }
        }
    }

    // 3) circular unitigs, handle remaining nodes that are still unvisited
    let circular_ids = crate::utils::order_keys(graph.nodes.keys().cloned());
    for id in &circular_ids {
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
            topology: 'l',
        });
    }

    // 4) build edges between unitigs based on original overlap graph

    // remove empty unitigs (can happen if there are isolated nodes with no edges)
    unitigs.retain(|u| !u.members.is_empty());
    let mut unitig_ends: HashMap<usize, (String, String)> = HashMap::new();
    // (start_node_id, end_node_id)


    let mut node_to_unitig: HashMap<String, usize> = HashMap::new();
    for u in &unitigs {
        let start = u.members.first().unwrap().node_id.clone();
        let end = u.members.last().unwrap().node_id.clone();
        unitig_ends.insert(u.id, (start, end));

        for m in &u.members {
            node_to_unitig.insert(m.node_id.clone(), u.id);
        }
    }

    let mut unitig_edge_map: HashMap<(usize, usize, char, char), UnitigEdge> = HashMap::new();
    for (source_id, node) in &graph.nodes {
        if let Some(&from_uid) = node_to_unitig.get(source_id) {
            let (from_start, from_end) = &unitig_ends[&from_uid];
            let from_is_start = source_id == from_start;
            let _from_is_end = source_id == from_end;
            let from_orient = if from_is_start { '+' } else { '-' };

            for e in &node.edges {
                if let Some(&to_uid) = node_to_unitig.get(&e.target_id) {
                    if from_uid == to_uid {
                        continue;
                    }

                    let (to_start, to_end) = &unitig_ends[&to_uid];
                    let to_is_start = e.target_id == *to_start;
                    let _to_is_end = e.target_id == *to_end;
                    let to_orient = if to_is_start { '+' } else { '-' };

                    let key = (from_uid, to_uid, from_orient, to_orient);
                    let entry = unitig_edge_map.entry(key).or_insert(UnitigEdge {
                        from: from_uid,
                        to: to_uid,
                        from_ori: from_orient,
                        to_ori: to_orient,
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

    let mut edges: Vec<UnitigEdge> = unitig_edge_map.into_values().collect();

    // remove one orientation of every unitig pair

    // helper function to strip the orientation from a node_id
    fn stripped_node_id(node_id: &str) -> String {
        let (read, start, end, _) = parse_node_id(node_id).unwrap();
        format!("{}:{}-{}", read, start, end)
    }

    // helper function to get the signature (vector of member names) of a unitig
    fn unitig_signature(u: &Unitig) -> Vec<String> {
        let mut sig: Vec<String> = u
            .members
            .iter()
            .map(|m| stripped_node_id(&m.node_id))
            .collect();

        // sort the signature to make it comparable
        sig.sort_unstable();
        sig
    }

    // vector of all reads that have been used in unitigs
    let mut seen: HashSet<String> = HashSet::new();
    // unitigs that should be kept
    let mut keep: HashSet<usize> = HashSet::new();

    unitigs.sort_unstable_by_key(|u| std::cmp::Reverse(u.members.len()));

    for u in unitigs.iter() {
        let sig = unitig_signature(u);

        for n in sig.iter() {
            if !seen.contains(n) {
                for s in sig.iter() {
                    seen.insert(s.clone());
                }
                keep.insert(u.id);
                break;
            }
        }
    }

    // remove unitigs
    unitigs.retain(|u| keep.contains(&u.id));

    // remove edges touching dropped unitigs
    edges.retain(|e| keep.contains(&e.from) && keep.contains(&e.to));

    if utils::has_seed() {
        edges.sort_unstable_by_key(|edge| (edge.from, edge.to, edge.from_ori, edge.to_ori));
    }

    // load the sequences of the reads used by the unitigs
    let needed: HashSet<&str> = unitigs
        .iter()
        .flat_map(|u| u.members.iter())
        .filter_map(|m| parse_node_id(&m.node_id).ok().map(|(name, ..)| name))
        .collect();
    let fastq_seqs = load_fastq_sequences(fastq_path, &needed).unwrap();

    for unitig in unitigs.iter_mut() {
        unitig.fasta_seq = Some(unitig_sequence(unitig, &fastq_seqs).unwrap());
    }
    set_unitig_topologies(&mut unitigs, overlaps);

    CompressedGraph { unitigs, edges }
}

fn load_fastq_sequences(
    fastq_path: &std::path::Path,
    needed: &HashSet<&str>,
) -> std::io::Result<HashMap<String, String>> {
    let mut seq_map = HashMap::new();
    let mut reader = utils::open_fastq_reader(fastq_path)?;
    while let Some((header, seq, _, _)) = utils::read_fastq_record(reader.as_mut())? {
        let id = header[1..].split_whitespace().next().unwrap_or("");
        if needed.contains(id) {
            seq_map.insert(id.to_string(), seq);
        }
    }
    Ok(seq_map)
}

/// Build the nucleotide sequence of `unitig` by concatenating node sequences and
/// removing the overlaps recorded in UnitigMember.edge.(target, edge_len).
pub fn unitig_sequence(
    unitig: &Unitig,
    fastq_seqs: &HashMap<String, String>,
) -> Result<String, String> {
    let get_seq = |node_id: &str| -> Result<String, String> {
        let (read_name, start, end, ori) = parse_node_id(node_id)?;
        let full_seq = fastq_seqs
            .get(read_name)
            .ok_or_else(|| format!("sequence for read '{read_name}' not found in FASTQ"))?;
        let slice = &full_seq[start as usize..end as usize];
        Ok(if ori == '+' { slice.to_string() } else { utils::rev_comp(slice) })
    };

    let mut out = String::new();
    for member in &unitig.members {
        // the last member has no target and contributes its full sequence
        let (target_id, edge_len) = &member.edge;
        let seq = get_seq(&member.node_id)?;
        if target_id.is_empty() {
            out.push_str(&seq);
        } else {
            let n = *edge_len as usize;
            if n > seq.len() {
                return Err(format!(
                    "edge length ({n}) larger than seq length ({}) for node {} to target {target_id}",
                    seq.len(),
                    member.node_id
                ));
            }
            out.push_str(&seq[..n]);
        }
    }
    Ok(out)
}