use crate::create_overlap_graph::OverlapGraph;
use crate::utils;
use std::collections::{HashMap, HashSet};

/// Check if the bigraph is synchronized:
/// 1. Every node has a reverse complement.
/// 2. Ingoing edges of every node correspond to outgoing edges of its reverse complement.
pub fn check_synchronization(g: &OverlapGraph) {
    for n in g.nodes.keys() {
        // compute reverse complement node
        let n_rc = utils::rc_node(n);

        // check that the reverse complement exists
        assert!(
            g.nodes.contains_key(&n_rc),
            "Reverse complement not found for {}. The bigraph is not synchronized.",
            n
        );

        // check that every outgoing edge has a counterpart in the reverse complement node
        if let Some(node) = g.nodes.get(n) {
            for e in &node.edges {
                let t = &e.target_id;
                let t_rc = utils::rc_node(t);

                // get the reverse complement node for t_rc
                if let Some(t_rc_node) = g.nodes.get(&t_rc) {
                    // Check if there's a matching edge from t_rc to n_rc
                    assert!(
                        t_rc_node.edges.iter().any(|e| e.target_id == n_rc),
                        "Corresponding edge not found for {}. The bigraph is not synchronized.",
                        n
                    );
                } else {
                    panic!(
                        "Reverse complement node {} missing for target {}. The bigraph is not synchronized.",
                        t_rc, t
                    );
                }
            }
        }
    }
}

/// Find weakly connected components ("clustered reads") of the graph.
pub fn weakly_connected_components(graph: &OverlapGraph) -> Vec<Vec<String>> {
    // Build an undirected adjacency list
    // Because it is undirected, we can move through both incoming and outgoing edges, meaning we can reach all nodes in a component.
    let mut adjacency_list: HashMap<String, Vec<String>> = HashMap::new();

    // ensure every node in graph.nodes has an entry
    for node_id in graph.nodes.keys() {
        adjacency_list.entry(node_id.clone()).or_default();
    }

    // Populate adjacency using outgoing edges (and add reverse edges to make undirected)
    for (source_id, node) in &graph.nodes {
        // Each edge is stored as EdgeInfo
        for e in &node.edges {
            let target_id = &e.target_id;
            adjacency_list
                .entry(source_id.clone())
                .or_default()
                .push(target_id.clone());
            adjacency_list
                .entry(target_id.clone())
                .or_default()
                .push(source_id.clone());
        }
    }

    // Find components with DFS/stack
    let mut visited: HashSet<String> = HashSet::new();
    let mut components: Vec<Vec<String>> = Vec::new();

    for start in adjacency_list.keys() {
        // check if already visited
        if visited.contains(start) {
            continue;
        }

        // new component
        let mut component: Vec<String> = Vec::new();
        let mut stack: Vec<String> = vec![start.clone()];
        visited.insert(start.clone());

        while let Some(current) = stack.pop() {
            component.push(current.clone());
            if let Some(neighbors) = adjacency_list.get(&current) {
                for neighbor in neighbors {
                    if !visited.contains(neighbor) {
                        visited.insert(neighbor.clone());
                        stack.push(neighbor.clone());
                    }
                }
            }
        }

        components.push(component);
    }

    components
}

/// Convenience: return component sizes sorted descending
pub fn component_sizes_sorted(graph: &OverlapGraph) -> Vec<usize> {
    let mut sizes: Vec<usize> = weakly_connected_components(graph)
        .into_iter()
        .map(|c| c.len())
        .collect();
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    sizes
}

/// Analyze node degrees to understand graph connectivity and compressibility
pub fn analyze_degrees(graph: &OverlapGraph) -> (HashMap<usize, usize>, HashMap<usize, usize>) {
    let mut indegree_dist: HashMap<usize, usize> = HashMap::new();
    let mut outdegree_dist: HashMap<usize, usize> = HashMap::new();

    // Count indegrees first
    let mut indegrees: HashMap<String, usize> = HashMap::new();
    for node in graph.nodes.values() {
        for e in &node.edges {
            *indegrees.entry(e.target_id.clone()).or_default() += 1;
        }
    }

    // Now collect distributions
    for node_id in graph.nodes.keys() {
        let in_deg = indegrees.get(node_id).copied().unwrap_or(0);
        let out_deg = graph.nodes.get(node_id).map(|n| n.edges.len()).unwrap_or(0);

        *indegree_dist.entry(in_deg).or_default() += 1;
        *outdegree_dist.entry(out_deg).or_default() += 1;
    }

    (indegree_dist, outdegree_dist)
}

/// Fraction of nodes that are compressible (in==1 && out==1) at the oriented-node level.
pub fn compressible_node_stats(graph: &OverlapGraph) -> (usize, usize, f64) {
    // compute indegrees
    let mut indegrees: HashMap<String, usize> = HashMap::new();
    for (src, node) in &graph.nodes {
        for e in &node.edges {
            *indegrees.entry(e.target_id.clone()).or_default() += 1;
        }
        // ensure src exists in map
        indegrees.entry(src.clone()).or_default();
    }

    let mut compressible = 0usize;
    let mut total = 0usize;
    for (node_id, node) in &graph.nodes {
        let in_deg = *indegrees.get(node_id).unwrap_or(&0);
        let out_deg = node.edges.len();
        total += 1;
        if in_deg == 1 && out_deg == 1 {
            compressible += 1;
        }
    }
    let frac = if total == 0 {
        0.0
    } else {
        compressible as f64 / total as f64
    };
    (compressible, total, frac)
}

/// Find tips and measure tip-lengths (walk forward from nodes with indeg==0)
/// max_walk limits how far we follow a chain (safety).
pub fn tip_length_distribution(graph: &OverlapGraph, max_walk: usize) -> Vec<usize> {
    // build indegrees first
    let mut indegrees: HashMap<String, usize> = HashMap::new();
    for (src, node) in &graph.nodes {
        for e in &node.edges {
            *indegrees.entry(e.target_id.clone()).or_default() += 1;
        }
        indegrees.entry(src.clone()).or_default();
    }

    let mut lengths: Vec<usize> = Vec::new();
    for start in graph.nodes.keys() {
        let in_deg = *indegrees.get(start).unwrap_or(&0);
        if in_deg != 0 {
            continue; // not a tip start
        }

        // follow forward while nodes are linear (in==1 && out==1)
        let mut cur = start.clone();
        let mut len = 0usize;
        let mut steps = 0usize;
        let mut visited_local: HashSet<String> = HashSet::new();
        while steps < max_walk {
            if visited_local.contains(&cur) {
                break;
            } // cycle safety
            visited_local.insert(cur.clone());

            let node = match graph.nodes.get(&cur) {
                Some(n) => n,
                None => break,
            };
            if node.edges.is_empty() {
                break;
            }
            // if more than one outgoing edge we stop counting the linear extension
            if node.edges.len() != 1 {
                break;
            }

            // move to next node
            let next = node.edges[0].target_id.clone();
            len += 1;
            steps += 1;

            let next_in = *indegrees.get(&next).unwrap_or(&0);
            if next_in != 1 {
                break;
            }
            cur = next;
        }
        lengths.push(len);
    }

    lengths
}

/// Simple branching summary: return top-k nodes by (in_deg + out_deg)
pub fn branching_summary(graph: &OverlapGraph, top_k: usize) -> Vec<(String, usize, usize)> {
    // compute indegrees
    let mut indegrees: HashMap<String, usize> = HashMap::new();
    for (src, node) in &graph.nodes {
        for e in &node.edges {
            *indegrees.entry(e.target_id.clone()).or_default() += 1;
        }
        indegrees.entry(src.clone()).or_default();
    }

    let mut v: Vec<(String, usize, usize, usize)> = Vec::new(); // id, in, out, sum
    for (id, node) in &graph.nodes {
        let in_deg = *indegrees.get(id).unwrap_or(&0);
        let out_deg = node.edges.len();
        let sum = in_deg + out_deg;
        if in_deg > 1 || out_deg > 1 {
            v.push((id.clone(), in_deg, out_deg, sum));
        }
    }
    // sort by sum desc
    v.sort_unstable_by(|a, b| b.3.cmp(&a.3));
    v.into_iter()
        .take(top_k)
        .map(|(id, in_deg, out_deg, _)| (id, in_deg, out_deg))
        .collect()
}
