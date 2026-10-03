use crate::create_overlap_graph::OverlapGraph;
use crate::utils;
/// Bubble removal module
/// using a "tour bus" style (BFS) search
/// 1. for each node `u` with out-degree >= 2, consider every pair of distinct
///   outgoing neighbors (v, w).
/// 2. perform a bounded BFS from v and from w (depth limited by `max_bubble_len`) to
///   discover meeting nodes m where the two searches converge.
/// 3. for the best meeting node (smallest combined depth), reconstruct the two
///   u->...->m paths, compute a simple path score (number of nodes, percentage identity) and pick
///   the higher-scoring path to keep. If the higher score is at least
///   `min_support_ratio * lower_score` (e.g. 1.1 to require 10% stronger), then the
///   lower-scoring path is removed (internal nodes removed, excluding u and m).
/// 4. removals are RC-aware: the reverse-complement node for each removed node is
///   also removed, and all incoming edges to removed nodes are purged.
use std::collections::{HashMap, HashSet, VecDeque};

/// Path metrics for scoring
#[derive(Clone, Default)]
struct PathMetrics {
    read_count: u32,
    total_overlap_len: u32,
    avg_identity: f64,
}

/// Bounded BFS from a start node
///   - parent map (node -> parent)
///   - depth map (node -> depth from start)
///   - path metrics map (node -> PathMetrics)
fn bfs_limited(
    graph: &OverlapGraph,
    start: &str,
    max_depth: usize,
) -> (
    HashMap<String, Option<String>>,
    HashMap<String, usize>,
    HashMap<String, PathMetrics>,
) {
    // initialize
    let mut parent: HashMap<String, Option<String>> = HashMap::new(); // map a node to its previous node in the path
    let mut depth: HashMap<String, usize> = HashMap::new(); // map a node to its depth from start
    let mut metrics: HashMap<String, PathMetrics> = HashMap::new(); // map a node to its path metrics

    // vecdeque is a double ended queue, allows efficient popping from front
    let mut q: VecDeque<String> = VecDeque::new();
    parent.insert(start.to_string(), None);
    depth.insert(start.to_string(), 0);
    metrics.insert(
        start.to_string(),
        PathMetrics {
            read_count: 1,
            total_overlap_len: 0,
            avg_identity: 1.0,
        },
    );
    q.push_back(start.to_string());

    // BFS loop (uses a fifo queue, push to the back, pop from the front)
    while let Some(cur) = q.pop_front() {
        // if depth limit is reached, stop adding nodes
        let cur_depth = *depth.get(&cur).unwrap_or(&0);
        if cur_depth >= max_depth {
            continue;
        }

        // expand neighbours
        if let Some(node) = graph.nodes.get(&cur) {
            for edge in &node.edges {
                // if the edge target is unseen, add it to parent/depth/metrics
                if !depth.contains_key(&edge.target_id) {
                    parent.insert(edge.target_id.clone(), Some(cur.clone()));
                    depth.insert(edge.target_id.clone(), cur_depth + 1);

                    // update path metrics, get metrics from current node and update them for the target node
                    let mut new_metrics = metrics.get(&cur).cloned().unwrap_or_default();
                    new_metrics.read_count += 1;
                    new_metrics.total_overlap_len += edge.overlap_len;
                    // update running average of identity
                    let old_weight = new_metrics.total_overlap_len as f64 - edge.overlap_len as f64;
                    let new_identity = if old_weight > 0.0 {
                        (new_metrics.avg_identity * old_weight
                            + edge.identity * edge.overlap_len as f64)
                            / new_metrics.total_overlap_len as f64
                    } else {
                        edge.identity
                    };
                    new_metrics.avg_identity = new_identity;

                    metrics.insert(edge.target_id.clone(), new_metrics);
                    q.push_back(edge.target_id.clone());
                }
            }
        }
    }

    (parent, depth, metrics)
}

/// Reconstruct path from source to sink using the parent map returned by bfs_limited
/// If sink is not reachable, returns an empty Vec
fn reconstruct_path(
    parent: &HashMap<String, Option<String>>,
    source: &str,
    sink: &str,
) -> Vec<String> {
    // the path will be reconstructed in reverse
    let mut path_rev: Vec<String> = Vec::new();
    let mut cur = sink.to_string();
    path_rev.push(cur.clone());

    while let Some(Some(p)) = parent.get(&cur) {
        cur = p.clone();
        path_rev.push(cur.clone());
    }
    // We expect the last item to be source
    if path_rev.last().map(|s| s.as_str()) != Some(source) {
        return Vec::new();
    }

    path_rev.reverse();
    path_rev
}

/// Remove simple bubbles in the overlap graph.
pub fn remove_bubbles(graph: &mut OverlapGraph, max_bubble_len: usize, min_support_ratio: f64) {
    if max_bubble_len == 0 {
        return;
    }

    // snapshot of nodes to iterate safely
    let node_keys: Vec<String> = crate::utils::order_keys(graph.nodes.keys().cloned());

    for n in node_keys.iter() {
        // get outgoing neighbors (clone so we don't borrow across mutation)
        let mut outgoing = match graph.nodes.get(n) {
            Some(n) => n
                .edges
                .iter()
                .map(|e| (e.target_id.clone(), e.edge_len))
                .collect::<Vec<_>>(),
            None => continue,
        };
        if crate::utils::has_seed() {
            outgoing.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
        }

        // skip if less than 2 outgoing edges
        if outgoing.len() < 2 {
            continue;
        }

        // consider every unordered pair of outgoing neighbors
        for i in 0..outgoing.len() {
            for j in (i + 1)..outgoing.len() {
                let start_a = &outgoing[i].0;
                let start_b = &outgoing[j].0;

                // skip identical starts, shouldn't happen though
                if start_a == start_b {
                    continue;
                }

                let (parent_a, depth_a, score_a) = bfs_limited(graph, start_a, max_bubble_len);
                let (parent_b, depth_b, score_b) = bfs_limited(graph, start_b, max_bubble_len);

                // find meeting nodes
                let reached_a: HashSet<&String> = depth_a.keys().collect();
                let reached_b: HashSet<&String> = depth_b.keys().collect();

                // check intersection of reached nodes
                let mut meetings: Vec<(&String, usize)> = Vec::new(); // (node, combined_depth)
                for node in reached_a.intersection(&reached_b) {
                    let d = depth_a.get(*node).unwrap_or(&usize::MAX)
                        + depth_b.get(*node).unwrap_or(&usize::MAX);
                    meetings.push((node.clone(), d));
                }

                // skip if no common node was reached
                if meetings.is_empty() {
                    continue;
                }

                // pick best sink node: minimal combined depth
                if crate::utils::has_seed() {
                    meetings.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(b.0)));
                } else {
                    meetings.sort_unstable_by_key(|k| k.1);
                }
                let (meet_node, _meet_depth) = meetings[0];

                // reconstruct paths start_a -> meet_node and start_b -> meet_node
                let path_a = reconstruct_path(&parent_a, start_a, meet_node);
                let path_b = reconstruct_path(&parent_b, start_b, meet_node);
                if path_a.is_empty() || path_b.is_empty() {
                    continue;
                }

                // get metrics for both paths
                let metrics_a = score_a.get(meet_node).cloned().unwrap_or_default();
                let metrics_b = score_b.get(meet_node).cloned().unwrap_or_default();
                let depth_a = *depth_a.get(meet_node).unwrap_or(&usize::MAX);
                let depth_b = *depth_b.get(meet_node).unwrap_or(&usize::MAX);

                // calculate composite scores
                // weight factors can be adjusted based on importance
                const OVERLAP_WEIGHT: f64 = 1.0;
                const IDENTITY_WEIGHT: f64 = 2.0;
                const READ_COUNT_WEIGHT: f64 = 1.5;

                let score_a = (metrics_a.total_overlap_len as f64 * OVERLAP_WEIGHT)
                    + (metrics_a.avg_identity * IDENTITY_WEIGHT * 100.0)
                    + (metrics_a.read_count as f64 * READ_COUNT_WEIGHT);

                let score_b = (metrics_b.total_overlap_len as f64 * OVERLAP_WEIGHT)
                    + (metrics_b.avg_identity * IDENTITY_WEIGHT * 100.0)
                    + (metrics_b.read_count as f64 * READ_COUNT_WEIGHT);

                // if both paths have no score (unexpected), skip
                if score_a == 0.0 && score_b == 0.0 {
                    continue;
                }

                // compare paths: higher score wins
                // if scores equal, shorter path (less depth) wins
                let (loser_path, winner_score) =
                    if score_a > score_b || (score_a == score_b && depth_a < depth_b) {
                        (path_b.clone(), score_a)
                    } else if score_b > score_a || (score_a == score_b && depth_b < depth_a) {
                        (path_a.clone(), score_b)
                    } else {
                        // Exactly equal - skip this bubble
                        continue;
                    };

                // require that the winner has enough support (based on score difference and min_support_ratio)
                let loser_score = if score_a > score_b { score_b } else { score_a };
                if loser_score * min_support_ratio > winner_score {
                    continue;
                }

                // nodes to remove: all nodes on loser_path excluding the sink node
                // also exclude the source node n, typically the path starts at the neighbor of n
                let mut to_remove: HashSet<String> = HashSet::new();
                for node in loser_path.into_iter() {
                    if node == *meet_node {
                        break;
                    }
                    // defensive: don't remove n
                    if node == *n {
                        continue;
                    }
                    to_remove.insert(node);
                }

                if to_remove.is_empty() {
                    continue;
                }

                // perform RC-aware removal
                utils::delete_nodes_and_edges(graph, &to_remove);

                break;
            }
        }
    }
}
