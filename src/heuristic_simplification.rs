use crate::create_overlap_graph::OverlapGraph;
use crate::utils;

/// Remove short edges (edges with overlap length below threshold)
/// For each node with multiple outgoing edges, keep only those with overlap length >= drop_ratio * best_overlap_len
pub fn remove_short_edges(graph: &mut OverlapGraph, drop_ratio: f64) -> usize {
    let mut n_short = 0;

    // iterate over a snapshot of current node keys (no mutation while iterating)
    let keys: Vec<String> = crate::utils::order_keys(graph.nodes.keys().cloned());

    for node_id in keys {
        // Get the outgoing edges for this node
        let edges_to_remove: Vec<String> = if let Some(node) = graph.nodes.get(&node_id) {
            // Skip if less than 2 outgoing edges
            if node.edges.len() < 2 {
                continue;
            }

            // Find the maximum overlap length
            let max_overlap = node.edges.iter().map(|e| e.overlap_len).max().unwrap_or(0);

            // Calculate threshold
            let threshold = (max_overlap as f64 * drop_ratio + 0.499) as u32;

            // Collect edges that are below the threshold
            node.edges
                .iter()
                .filter(|e| e.overlap_len < threshold)
                .map(|e| e.target_id.clone())
                .collect()
        } else {
            continue;
        };

        // Remove the short edges
        for target_id in edges_to_remove {
            if let Some(node) = graph.nodes.get_mut(&node_id) {
                node.remove_edge(&target_id);
                n_short += 1;
                // Remove the reverse complement edge
                let rc_target_id = utils::rc_node(&target_id);
                let rc_node_id = utils::rc_node(&node_id);
                if let Some(rc_target) = graph.nodes.get_mut(&rc_target_id) {
                    rc_target.remove_edge(&rc_node_id);
                    n_short += 1;
                }
            }
        }
    }
    n_short
    //eprintln!("remove short edges: removed {} short edges", n_short);
}

/// Delete multi-arcs: when a node has multiple arcs to the same target,
/// keep the best (largest overlap_len, then highest identity) and remove the others.
pub fn remove_multi_edges(graph: &mut OverlapGraph) -> usize {
    use std::collections::HashMap;

    let mut n_multi = 0usize;
    let keys: Vec<String> = crate::utils::order_keys(graph.nodes.keys().cloned());

    for src in keys {
        let edges_snapshot = match graph.nodes.get(&src) {
            Some(n) => n.edges.clone(),
            None => continue,
        };
        if edges_snapshot.len() < 2 {
            continue;
        }

        let mut counts: HashMap<String, usize> = HashMap::new();
        let mut best_idx: HashMap<String, usize> = HashMap::new();
        for (i, e) in edges_snapshot.iter().enumerate() {
            *counts.entry(e.target_id.clone()).or_insert(0) += 1;
            match best_idx.get(&e.target_id) {
                Some(&bi) => {
                    let be = &edges_snapshot[bi];
                    if e.overlap_len > be.overlap_len
                        || (e.overlap_len == be.overlap_len && e.identity > be.identity)
                    {
                        best_idx.insert(e.target_id.clone(), i);
                    }
                }
                None => {
                    best_idx.insert(e.target_id.clone(), i);
                }
            }
        }

        if edges_snapshot.len() == best_idx.len() {
            continue;
        }

        let mut chosen: Vec<(String, crate::create_overlap_graph::EdgeInfo)> = best_idx
            .into_iter()
            .map(|(t, i)| (t, edges_snapshot[i].clone()))
            .collect();
        chosen.sort_by(|a, b| a.0.cmp(&b.0));
        let new_edges: Vec<crate::create_overlap_graph::EdgeInfo> =
            chosen.into_iter().map(|(_t, e)| e).collect();

        let mut removed_targets: Vec<String> = Vec::new();
        for (t, cnt) in counts.into_iter() {
            if cnt > 1 {
                removed_targets.push(t);
                n_multi += cnt - 1;
            }
        }

        if let Some(node_mut) = graph.nodes.get_mut(&src) {
            node_mut.edges = new_edges;
        }

        for tgt in removed_targets {
            if let Some(tnode) = graph.nodes.get_mut(&tgt) {
                tnode.remove_edge(&src);
            }
            let src_rc = utils::rc_node(&src);
            let tgt_rc = utils::rc_node(&tgt);
            if let Some(tnode_rc) = graph.nodes.get_mut(&tgt_rc) {
                tnode_rc.remove_edge(&src_rc);
            }
            if let Some(src_rc_node) = graph.nodes.get_mut(&src_rc) {
                src_rc_node.remove_edge(&tgt_rc);
            }
        }
    }

    n_multi
}