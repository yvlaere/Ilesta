use crate::alignment_filtering::Overlap;
/// Overlap graph creation module
/// read overlaps from alignment filtering module and build the overlap graph
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{BufWriter, Write};

/// Edge info containing all the metrics we track
#[derive(Clone)]
pub struct EdgeInfo {
    pub target_id: String,
    pub edge_len: u32,
    pub overlap_len: u32,
    pub identity: f64,
}

/// A node in the overlap graph. Earch read is represented by two nodes: "<read_name>+" and "<read_name>-"
/// One for the origininal orientation and one for the reverse complement
/// Each node has directed edges to other nodes with associated edge lengths
pub struct Node {
    node_id: String,
    pub edges: Vec<EdgeInfo>,
}

impl Node {
    /// Create a new node with given id and no edges
    fn new(node_id: String) -> Self {
        Self {
            node_id,
            edges: Vec::new(),
        }
    }

    /// Add a directed edge to a node, if an edge to the target already exists, we ignore (avoid duplicates)
    fn add_edge(&mut self, target_node: &str, edge_len: u32, overlap_len: u32, identity: f64) {
        if self.edges.iter().any(|e| e.target_id == target_node) {
            // multiple edges to the same target
            // currently impossible due to the overlap filtering
            // handling might change in the future to handle this case
            // silently ignore duplicate edges but log for debugging
            eprintln!(
                "Warning: duplicate edge {} -> {} ignored",
                self.node_id, target_node
            );
            return;
        }
        self.edges.push(EdgeInfo {
            target_id: target_node.to_owned(),
            edge_len,
            overlap_len,
            identity,
        });
    }

    /// Remove a directed edge to the node with target_node id
    pub fn remove_edge(&mut self, target_node: &str) {
        if let Some(pos) = self.edges.iter().position(|e| e.target_id == target_node) {
            self.edges.swap_remove(pos);
        }
    }

    /// Sort edges by length (ascending).
    pub fn sort_edges(&mut self) {
        self.edges.sort_unstable_by_key(|e| e.edge_len);
    }
}

/// Overlap graph containing nodes keyed by their node id
pub struct OverlapGraph {
    pub nodes: HashMap<String, Node>,
}

impl OverlapGraph {
    /// Create a new empty overlap graph
    fn new() -> Self {
        Self {
            nodes: HashMap::new(),
        }
    }

    /// Add a node to the graph if it does not already exist, if it already exists do nothing
    fn add_node(&mut self, node_id: String) {
        self.nodes
            .entry(node_id.clone())
            .or_insert_with(|| Node::new(node_id));
    }

    /// Add a directed edge from from_id to to_id with given edge length and metrics
    fn add_edge(
        &mut self,
        from_id: &str,
        to_id: &str,
        edge_len: u32,
        overlap_len: u32,
        identity: f64,
    ) {
        // ensure nodes exist
        if !self.nodes.contains_key(from_id) || !self.nodes.contains_key(to_id) {
            panic!("add_edge: nodes must exist before adding edge");
        }
        if let Some(node) = self.nodes.get_mut(from_id) {
            node.add_edge(to_id, edge_len, overlap_len, identity);
        }
    }

    /// Write the overlap graph to a DOT file for visualization
    pub fn write_dot<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        // helper functions
        fn escape_dot(s: &str) -> String {
            s.replace('\\', "\\\\").replace('"', "\\\"")
        }

        fn degree_color(degree: usize) -> &'static str {
            match degree {
                0 => "gray",
                1 => "black",
                2 => "blue",
                3..=5 => "orange",
                _ => "red",
            }
        }

        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        // gather stats
        let mut outdegrees: HashMap<String, usize> = HashMap::new();
        let mut indegrees: HashMap<String, usize> = HashMap::new();
        for node in self.nodes.values() {
            *outdegrees.entry(node.node_id.clone()).or_insert(0) += node.edges.len();
            for e in &node.edges {
                *indegrees.entry(e.target_id.clone()).or_insert(0) += 1;
            }
        }

        writeln!(w, "digraph OverlapGraph {{")?;
        writeln!(w, "  rankdir=LR;")?;
        writeln!(w, "  node [shape=box fontname=\"Helvetica\"];")?;
        writeln!(w, "  edge [fontname=\"Helvetica\"];")?;
        writeln!(w)?;

        // Emit nodes explicitly (optional but useful for styling later)
        for node_id in self.nodes.keys() {
            writeln!(
                w,
                "  \"{}\" [style=filled fillcolor={} ];",
                escape_dot(node_id),
                degree_color(*outdegrees.get(node_id).unwrap_or(&0))
            )?;
        }

        writeln!(w)?;

        // Emit edges
        for node in self.nodes.values() {
            let from = escape_dot(&node.node_id);
            for e in &node.edges {
                let to = escape_dot(&e.target_id);
                writeln!(
                    w,
                    "  \"{}\" -> \"{}\";",
                    from,
                    to //"  \"{}\" -> \"{}\" [label=\"len={} ovl={} id={:.3}\"];",
                       //from, to, e.edge_len, e.overlap_len, e.identity
                )?;
            }
        }

        writeln!(w, "}}")?;
        Ok(())
    }
}

pub struct AlignmentFilteringOutput {
    pub overlaps: HashMap<(usize, usize), Overlap>,
}

impl AlignmentFilteringOutput {
    // serialize the overlaps
    pub fn serialize_overlaps(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::create(path)?;
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, &self.overlaps)?;

        Ok(())
    }
}

/// Build overlap graph from overlaps
pub fn run_create_overlap_graph(
    overlaps: &HashMap<(usize, usize), Overlap>,
) -> Result<OverlapGraph, io::Error> {
    println!("=== OVERLAP GRAPH CREATION ===");
    let mut g = OverlapGraph::new();

    for ((_query_id, _target_id), o) in overlaps.iter() {
        // add overlap to the graph

        // original orientation
        g.add_node(o.source_name.clone());
        g.add_node(o.sink_name.clone());
        g.add_edge(
            &o.source_name,
            &o.sink_name,
            o.edge_len,
            o.overlap_len,
            o.identity,
        );

        // reverse complement counterpart:
        g.add_node(o.rc_source_name.clone());
        g.add_node(o.rc_sink_name.clone());
        g.add_edge(
            &o.rc_source_name,
            &o.rc_sink_name,
            o.rc_edge_len,
            o.overlap_len,
            o.identity,
        );
    }

    // graph stats
    let edge_count: usize = g.nodes.values().map(|n| n.edges.len()).sum();
    let node_count = g.nodes.len();
    let node_to_edge_ratio = node_count as f64 / edge_count as f64;
    println!("Graph nodes: {}", node_count);
    println!("Graph edges: {}", edge_count);
    println!("Node to edge ratio: {:.4}", node_to_edge_ratio);
    println!("=== OVERLAP GRAPH CREATION FINISHED ===");
    return Ok(g);
}
