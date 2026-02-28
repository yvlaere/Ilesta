/// Transitive edge reduction module
/// transitive edges are redundant edges that don't add any information to the graph
/// Say read 1 overlaps with read 2 and read 2 overlaps with read 3 and read 1 also overlaps with read 3, then this last overlap is redundant, represented by a transitive edge
/// Algorithm based on https://doi.org/10.1093/bioinformatics/bti1114
use crate::create_overlap_graph::OverlapGraph;
use crate::utils;

use std::collections::{HashMap, HashSet};

/// Enum for node marking (Vacant, In-play, Eliminated)
#[derive(Copy, Clone, PartialEq)]
enum Mark {
    Vacant,
    InPlay,
    Eliminated,
}

/// Reduce transitive edges
pub fn reduce_transitive_edges(g: &mut OverlapGraph, fuzz: u32) {
    // prepare node list to iterate deterministically and avoid borrow conflicts
    let node_keys: Vec<String> = g.nodes.keys().cloned().collect();

    // mark: per-node status (Vacant/InPlay/Eliminated)
    let mut mark: HashMap<String, Mark> = HashMap::with_capacity(g.nodes.len());

    // reduced set: node pairs (from, to) that should be removed
    let mut reduced: HashSet<(String, String)> = HashSet::new();

    // initialize all nodes as vacant and ensure edges are sorted ascending by length
    for n in &node_keys {
        mark.insert(n.clone(), Mark::Vacant);
        if let Some(node) = g.nodes.get_mut(n) {
            node.sort_edges();
        }
    }

    // main loop: For every node compare nodes encountered two steps into the future with those encountered one step into the future
    for n1 in &node_keys {
        // skip if node not present (may have been removed) or no outgoing edges
        let out_edges = match g.nodes.get(n1) {
            Some(node) => &node.edges,
            None => continue,
        };
        if out_edges.is_empty() {
            continue;
        }

        // 1) mark all direct neighbors of n1 as InPlay
        for e in out_edges.iter() {
            mark.insert(e.target_id.clone(), Mark::InPlay);
        }

        // 2) compute longest outgoing edge length from n1 (last after sort) + fuzz
        let longest = {
            let last_len = out_edges.last().unwrap().edge_len as u64;
            last_len + fuzz as u64
        };

        // 3) For each n2 (outgoing from n1), check n2->n3 edges
        for e_n2 in out_edges.iter() {
            let n2 = &e_n2.target_id;
            let len_n1n2 = e_n2.edge_len;

            // skip n2 if not InPlay
            if mark.get(n2).copied() != Some(Mark::InPlay) {
                continue;
            }

            // get node n2
            if let Some(node2) = g.nodes.get(n2) {
                for e_n3 in node2.edges.iter() {
                    let n3 = &e_n3.target_id;
                    let len_n2n3 = e_n3.edge_len;
                    // if path length n1->n2->n3 <= longest then candidate for elimination
                    let path_len = len_n2n3 as u64 + len_n1n2 as u64;
                    if path_len <= longest
                        && mark.get(n3).copied() == Some(Mark::InPlay) {
                            mark.insert(n3.clone(), Mark::Eliminated);
                        }
                }
            }
        }

        // 4) Additional rule: if n2->n3 is very small (< fuzz) or is the smallest outgoing edge of n2,
        // then eliminate n3 if it is InPlay.
        for e in out_edges.iter() {
            let n2 = &e.target_id;
            if let Some(node2) = g.nodes.get(n2) {
                // find min outgoing length for n2, if any
                let min_len_opt = node2.edges.iter().map(|e| e.edge_len).min();

                for e_n3 in node2.edges.iter() {
                    let n3 = &e_n3.target_id;
                    let len_n2n3 = e_n3.edge_len;
                    let do_eliminate = if len_n2n3 < fuzz {
                        true
                    } else if let Some(min_len) = min_len_opt {
                        len_n2n3 == min_len
                    } else {
                        false
                    };

                    if do_eliminate && mark.get(n3).copied() == Some(Mark::InPlay) {
                        mark.insert(n3.clone(), Mark::Eliminated);
                    }
                }
            }
        }

        // 5) Mark edges from n1 to eliminated nodes for removal, then reset marks of the direct neighbors to vacant
        for e_n2 in out_edges.iter() {
            let n2 = &e_n2.target_id;
            if mark.get(n2).copied() == Some(Mark::Eliminated) {
                reduced.insert((n1.clone(), n2.clone()));
            }
            // reset mark back to vacant for next iteration
            mark.insert(n2.clone(), Mark::Vacant);
        }
    } // end for n1

    // 6) Actually remove reduced edges from the graph (and remove reverse-complement counterparts if present)

    // Collect all edges to remove first
    let mut edges_to_remove: Vec<(String, String)> = Vec::new();

    for n1 in g.nodes.keys() {
        if let Some(node) = g.nodes.get(n1) {
            for e in &node.edges {
                let n2 = &e.target_id;
                if reduced.contains(&(n1.clone(), n2.clone())) {
                    // Add both the edge and its reverse complement
                    edges_to_remove.push((n1.clone(), n2.clone()));
                    edges_to_remove.push((utils::rc_node(n2), utils::rc_node(n1)));
                }
            }
        }
    }

    // Now remove all edges in a separate pass
    for (from, to) in edges_to_remove {
        if let Some(node) = g.nodes.get_mut(&from) {
            node.remove_edge(&to);
        }
    }
}
