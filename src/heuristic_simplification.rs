use crate::create_overlap_graph::OverlapGraph;
use crate::utils;

/// Remove short edges (edges with overlap length below threshold)
/// For each node with multiple outgoing edges, keep only those with overlap length >= drop_ratio * best_overlap_len
pub fn remove_short_edges(graph: &mut OverlapGraph, drop_ratio: f64) -> usize {
    let mut n_short = 0;

    // iterate over a snapshot of current node keys (no mutation while iterating)
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();

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
                //if let Some(target_node) = graph.nodes.get_mut(&target_id) {
                //    target_node.remove_edge(&node_id);
                //}
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
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();

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

/// Cut small bi-loops: patterns where v->...->x and w->v, w->x exist
/// If overlap(w->v) > overlap(w->x), remove the w->x edge (keep the longer path)
pub fn cut_biloop(graph: &mut OverlapGraph, max_ext: usize) {
    let mut cnt = 0;

    // iterate over a snapshot of current node keys
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();

    for v in keys {
        // Skip if node doesn't have multiple outgoing edges
        if let Some(v_node) = graph.nodes.get(&v) {
            if v_node.edges.len() < 2 {
                continue;
            }
        } else {
            continue;
        }

        // Try to extend from this node up to max_ext steps
        let extended_path = extend_path(graph, &v, max_ext);
        if extended_path.is_empty() || extended_path.len() < 2 {
            continue;
        }

        // Get the last node in the extended path (x)
        let x = &extended_path[extended_path.len() - 1];

        // Find incoming edges to v (w->v edges)
        let incoming_to_v: Vec<(String, u32)> = graph
            .nodes
            .iter()
            .filter_map(|(node_id, node)| {
                node.edges
                    .iter()
                    .find(|e| &e.target_id == &v)
                    .map(|e| (node_id.clone(), e.overlap_len))
            })
            .collect();

        // For each incoming node w, check if it also has an edge to x
        for (w, ov) in incoming_to_v {
            if let Some(w_node) = graph.nodes.get(&w) {
                // Find overlap from w to x
                if let Some(edge_to_x) = w_node.edges.iter().find(|e| &e.target_id == x) {
                    let ox = edge_to_x.overlap_len;

                    // If overlap(w->v) > overlap(w->x), remove w->x edge
                    if ov > ox {
                        if let Some(w_mut) = graph.nodes.get_mut(&w) {
                            w_mut.remove_edge(x);
                        }
                        // Remove reverse edge
                        if let Some(x_mut) = graph.nodes.get_mut(x) {
                            x_mut.remove_edge(&w);
                        }
                        cnt += 1;
                    }
                }
            }
        }
    }

    if cnt > 0 {
        //eprintln!("[heuristic_simplification::cut_biloop] cut {} small bi-loops", cnt);
    }
}

/// Helper function to extend a path from a starting node up to max_ext edges
/// Returns the sequence of nodes visited (including start node)
fn extend_path(graph: &OverlapGraph, start: &str, max_ext: usize) -> Vec<String> {
    let mut path = vec![start.to_string()];
    let mut current = start.to_string();

    for _ in 0..max_ext {
        // Get the single outgoing edge (if it exists and is unique)
        if let Some(node) = graph.nodes.get(&current) {
            // Only extend if there's exactly one outgoing edge
            if node.edges.len() == 1 {
                let next = &node.edges[0].target_id;
                // Avoid cycles
                if !path.contains(next) {
                    path.push(next.clone());
                    current = next.clone();
                } else {
                    break;
                }
            } else {
                // Stop if we hit a branching point
                break;
            }
        } else {
            break;
        }
    }

    path
}

/// Cut internal sequences: when there's a simple linear chain between two branching
/// nodes, remove the internal sequences (reads) found along that chain up to
/// `max_ext` steps.
pub fn cut_internal(graph: &mut OverlapGraph, max_ext: usize) {
    // build indegree map
    let mut indegree: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for id in graph.nodes.keys() {
        indegree.insert(id.clone(), 0);
    }
    for (_id, node) in &graph.nodes {
        for e in &node.edges {
            *indegree.entry(e.target_id.clone()).or_default() += 1;
        }
    }

    let mut removed_reads = 0usize;

    // snapshot of keys to avoid borrowing while mutating
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();

    for v in keys {
        // snapshot outgoing targets for v to avoid borrowing graph while mutating
        let outgoing_targets: Vec<String> = match graph.nodes.get(&v) {
            Some(n) => n.edges.iter().map(|e| e.target_id.clone()).collect(),
            None => continue,
        };
        if outgoing_targets.len() < 2 {
            continue;
        }

        // for each outgoing target, try to follow a simple chain
        for target in outgoing_targets {
            let mut path: Vec<String> = Vec::new();
            let mut cur = target.clone();
            let mut steps = 0usize;

            // follow while nodes are simple (indegree==1 && outdegree==1)
            while steps < max_ext {
                // ensure node exists
                let cur_node = match graph.nodes.get(&cur) {
                    Some(n) => n,
                    None => break,
                };
                let in_deg = *indegree.get(&cur).unwrap_or(&0);
                let out_deg = cur_node.edges.len();

                // stop extension if this node is not a simple internal node
                if in_deg != 1 || out_deg != 1 {
                    break;
                }

                // record this internal node and advance
                path.push(cur.clone());
                let next = cur_node.edges[0].target_id.clone();
                // avoid cycles
                if path.contains(&next) {
                    break;
                }
                cur = next;
                steps += 1;
            }

            // check if chain ended at a branching node (both ends branching)
            if path.is_empty() {
                continue;
            }
            if let Some(end_node) = graph.nodes.get(&cur) {
                if end_node.edges.len() < 2 {
                    continue;
                }
            } else {
                continue;
            }

            // delete the internal reads (both orientations)
            for internal in path {
                if internal.len() < 1 {
                    continue;
                }
                let read = &internal[..internal.len() - 1];
                let plus = format!("{}+", read);
                let minus = format!("{}-", read);
                let removed_plus = graph.nodes.remove(&plus);
                let removed_minus = graph.nodes.remove(&minus);
                if removed_plus.is_some() || removed_minus.is_some() {
                    removed_reads += 1;
                }
            }
        }
    }

    // cleanup: remove edges that point to missing nodes
    let existing: std::collections::HashSet<String> = graph.nodes.keys().cloned().collect();
    for node in graph.nodes.values_mut() {
        node.edges.retain(|e| existing.contains(&e.target_id));
    }

    if removed_reads > 0 {
        //eprintln!("[heuristic_simplification::cut_internal] cut {} internal sequences", removed_reads);
    }
}

/// Remove low identity from nodes with multiple outgoing edges
pub fn remove_weak(graph: &mut OverlapGraph) {
    // iterate over a snapshot of current node keys (no mutation while iterating)
    let keys: Vec<String> = graph.nodes.keys().cloned().collect();
    for n in keys.into_iter() {
        // check the amount of outgoing edges
        let outgoing = match graph.nodes.get(&n) {
            Some(n) => n
                .edges
                .iter()
                .map(|e| (e.target_id.clone(), e.identity))
                .collect::<Vec<_>>(),
            None => continue,
        };

        // skip if less than 2 outgoing edges
        if outgoing.len() < 2 {
            continue;
        }

        // only keep edge with highest identity, remove others
        let mut max_identity: f64 = -1.0;
        let mut best_target: Option<String> = None;
        for (target_id, identity) in outgoing.iter() {
            if *identity > max_identity {
                max_identity = *identity;
                best_target = Some(target_id.clone());
            }
        }

        // remove all edges except the best one
        for (target_id, _) in outgoing.iter() {
            if Some(target_id) != best_target.as_ref() {
                if let Some(node) = graph.nodes.get_mut(&n) {
                    node.remove_edge(target_id);
                    // remove the reverse edge as well
                    if let Some(target_node) = graph.nodes.get_mut(target_id) {
                        target_node.remove_edge(&n);
                    }
                }
            }
        }
    }
}
