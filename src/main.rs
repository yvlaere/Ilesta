mod align_reads;
mod alignment_filtering;
mod bubble_removal;
mod cli;
mod compress_graph;
mod consensus;
mod create_overlap_graph;
mod heuristic_simplification;
mod plasmid_rescue;
mod scaffold;
mod segment_dedup;
mod tip_trimming;
mod transitive_edge_reduction;
mod utils;

use alignment_filtering::Overlap;
use clap::Parser;
use cli::{AssembleArgs, Cli, Commands};
use create_overlap_graph::OverlapGraph;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

fn main() -> Result<(), Box<dyn Error>> {
    match Cli::parse().command {
        Commands::Align(args) => {
            let out_dir = Path::new(&args.output_dir);
            std::fs::create_dir_all(out_dir)?;
            check_exists(&args.reads.reads_fq)?;
            align_reads::align_reads(&args.reads, out_dir)?;
        }
        Commands::AlignmentFiltering(args) => {
            utils::set_seed(args.seed);
            let out_dir = Path::new(&args.output_dir);
            std::fs::create_dir_all(out_dir)?;
            let out = alignment_filtering::run_alignment_filtering(&out_dir.join(&args.paf), &args.overlap)?;
            out.serialize_overlaps(&args.output_overlaps)?;
            println!("Wrote overlaps to {}", args.output_overlaps);
        }
        Commands::Scaffold(args) => {
            let out = Path::new(&args.output);
            let work_dir = out.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
            let reads = Path::new(&args.reads_fq);
            let mut graph = compress_graph::CompressedGraph::from_fasta(Path::new(&args.contigs))?;
            let stats = scaffold::scaffold_graph(&mut graph, reads, work_dir, &args.scaffold, args.threads, args.read_type)?;
            println!("Applied {} bridges ({} circularised)", stats.accepted, stats.circularised);
            if args.polish_rounds > 0 {
                let params = consensus::ConsensusParams {
                    rounds: args.polish_rounds,
                    window_len: 500,
                    min_depth: 3,
                    max_depth: 25,
                    threads: args.threads,
                    read_type: args.read_type,
                };
                consensus::polish_graph(&mut graph, reads, work_dir, &params)?;
            }
            graph.write_fasta(out)?;
        }
        Commands::Assemble(args) => assemble(&args)?,
    }
    Ok(())
}

fn check_exists(path: &str) -> Result<(), Box<dyn Error>> {
    if Path::new(path).exists() {
        Ok(())
    } else {
        Err(format!("Reads file {path} does not exist").into())
    }
}

fn edge_count(graph: &OverlapGraph) -> usize {
    graph.nodes.values().map(|n| n.edges.len()).sum()
}

/// Iterative graph cleanup: transitive reduction, multi-edges, short edges, bubbles,
/// isolated nodes and tips.
fn cleanup_graph(graph: &mut OverlapGraph, args: &AssembleArgs) {
    println!("\n=== STARTING GRAPH CLEANUP ===");
    utils::check_synchronization(graph);
    for iteration in 1..=args.cleanup_iterations {
        println!("\n=== Cleanup Iteration {iteration} ===");
        let before = edge_count(graph);
        transitive_edge_reduction::reduce_transitive_edges(graph, args.fuzz);
        println!("Removed {} edges with transitive edge reduction", before - edge_count(graph));
        println!("Removed {} multi-edges", heuristic_simplification::remove_multi_edges(graph));
        println!("Removed {} short edges", heuristic_simplification::remove_short_edges(graph, args.short_edge_ratio));
        let before = graph.nodes.len();
        bubble_removal::remove_bubbles(graph, args.max_bubble_length as usize, args.min_support_ratio);
        println!("Removed {} bubble nodes (including RCs)", before - graph.nodes.len());

        // nodes without any edge
        let targets: HashSet<&str> =
            graph.nodes.values().flat_map(|n| n.edges.iter().map(|e| e.target_id.as_str())).collect();
        let isolated: HashSet<String> = graph
            .nodes
            .iter()
            .filter(|(id, n)| n.edges.is_empty() && !targets.contains(id.as_str()))
            .map(|(id, _)| id.clone())
            .collect();
        utils::delete_nodes_and_edges(graph, &isolated);
        println!("Removed {} isolated nodes", isolated.len());

        let before = graph.nodes.len();
        tip_trimming::trim_tips(graph, args.max_tip_len as usize);
        println!("Removed {} nodes by tip trimming", before - graph.nodes.len());
    }
    utils::check_synchronization(graph);
    println!("\n=== GRAPH CLEANUP COMPLETE ===");
    println!("Final graph has {} nodes and {} edges", graph.nodes.len(), edge_count(graph));
}

fn assemble(args: &AssembleArgs) -> Result<(), Box<dyn Error>> {
    utils::set_seed(args.seed);
    let out_dir = Path::new(&args.output_dir);
    std::fs::create_dir_all(out_dir)?;
    check_exists(&args.reads.reads_fq)?;
    let (threads, read_type) = (args.reads.threads, args.reads.read_type);

    println!("=== READ FILTERING AND ALIGNMENT ===");
    let reads = align_reads::align_reads(&args.reads, out_dir)?;
    println!("Filtered read count: {}", reads.filtered_read_count);
    println!("Reads used in primary assembly: {}", reads.primary_read_count);

    // overlaps: provided, or filtered from the alignments
    let overlaps_path: PathBuf = match &args.overlaps {
        Some(path) => PathBuf::from(path),
        None => {
            let path = out_dir.join(format!("{}.overlaps.bin", args.output_prefix));
            let out = alignment_filtering::run_alignment_filtering(&out_dir.join(&args.reads.paf), &args.overlap)?;
            out.serialize_overlaps(path.to_str().ok_or("invalid output path")?)?;
            println!("Wrote overlaps to {}", path.display());
            path
        }
    };
    let mut overlaps: HashMap<(usize, usize), Overlap> = bincode::deserialize_from(BufReader::new(File::open(&overlaps_path)?))?;
    let mut graph = create_overlap_graph::run_create_overlap_graph(&overlaps)?;
    cleanup_graph(&mut graph, args);
    if args.write_dot {
        graph.write_dot(out_dir.join("graph.dot"))?;
    }

    println!("\n=== COMPRESSING UNITIGS ===");
    let out_path = out_dir.join(format!("{}.fa", args.output_prefix));
    let mut assembly_reads = reads.primary_reads.clone();
    let mut compressed = compress_graph::compress_unitigs(&graph, &overlaps, &assembly_reads);
    println!("Primary assembly produced {} unitigs", compressed.unitigs.len());

    let mut rescued_nodes = HashSet::new();
    if args.rescue_enabled() {
        println!("\n=== PLASMID RESCUE ===");
        compressed.write_fasta(&out_path)?;
        let rescue = plasmid_rescue::run_rescue_stage(
            &reads.filtered_reads, &out_path, out_dir, &args.reads, &args.overlap, &mut graph, &mut overlaps,
        )?;
        if let Some(rescue) = rescue {
            let combined = out_dir.join(format!("{}.combined.fq", args.output_prefix));
            let mut file = File::create(&combined)?;
            std::io::copy(&mut File::open(&reads.primary_reads)?, &mut file)?;
            std::io::copy(&mut File::open(&rescue.reads)?, &mut file)?;
            assembly_reads = combined;
            rescued_nodes = rescue.node_ids;
            compressed = compress_graph::compress_unitigs(&graph, &overlaps, &assembly_reads);
        }
    }

    compressed.write_fasta(&out_path)?;
    let removed = segment_dedup::remove_contained_segments(&mut compressed, &out_path, out_dir, threads)?;
    println!("Removed {removed} contained assembly segments");

    if args.scaffold.rounds > 0 {
        println!("\n=== CONTIG BRIDGING ===");
        let stats = scaffold::scaffold_graph(&mut compressed, &assembly_reads, out_dir, &args.scaffold, threads, read_type)?;
        println!("Bridging: {} bridges applied, {} circularised", stats.accepted, stats.circularised);
    }
    if args.consensus_rounds > 0 {
        println!("\n=== CONSENSUS POLISHING ===");
        let params = consensus::ConsensusParams {
            rounds: args.consensus_rounds,
            window_len: args.consensus_window_len,
            min_depth: args.consensus_min_depth,
            max_depth: args.consensus_max_depth,
            threads,
            read_type,
        };
        let stats = consensus::polish_graph(&mut compressed, &assembly_reads, out_dir, &params)?;
        println!("Consensus: {} rounds, {} -> {} bases", stats.rounds_run, stats.bases_before, stats.bases_after);
    }
    if assembly_reads != reads.primary_reads {
        std::fs::remove_file(&assembly_reads)?;
    }

    compressed.write_fasta(&out_path)?;
    let gfa_path = out_dir.join(format!("{}.gfa", args.output_prefix));
    compressed.write_gfa(gfa_path.to_str().ok_or("invalid output path")?)?;
    println!("Assembly produced {} unitigs (written to {} and {})", compressed.unitigs.len(), out_path.display(), gfa_path.display());
    if !rescued_nodes.is_empty() {
        let rescued: Vec<_> = compressed
            .unitigs
            .iter()
            .filter(|u| u.members.iter().any(|m| rescued_nodes.contains(&m.node_id)))
            .collect();
        let bases: usize = rescued.iter().map(|u| u.fasta_seq.as_ref().map_or(0, String::len)).sum();
        println!("Rescued unitigs: {} ({bases} bases)", rescued.len());
    }
    println!("\n=== ASSEMBLY COMPLETE ===");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    fn fixture_overlaps(reverse: bool) -> HashMap<(usize, usize), Overlap> {
        let mut pairs = vec![(0, 2), (1, 2), (2, 5), (0, 5), (1, 5), (3, 5), (4, 5)];
        if reverse {
            pairs.reverse();
        }
        pairs
            .into_iter()
            .map(|(source, target)| {
                let overlap = Overlap {
                    source_name: format!("read{}:0-8+", source),
                    sink_name: format!("read{}:0-8+", target),
                    rc_source_name: format!("read{}:0-8-", target),
                    rc_sink_name: format!("read{}:0-8-", source),
                    edge_len: if source < 2 && target == 5 { 2 } else { 1 },
                    rc_edge_len: 1,
                    overlap_len: 6,
                    identity: 99.0,
                };
                ((source, target), overlap)
            })
            .collect()
    }

    #[test]
    fn seeded_outputs_are_reproducible() {
        let directory = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("target")
            .join(format!("seeded_reproducibility_{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        let fastq_path = directory.join("reads.fq");
        let fasta_path = directory.join("unitigs.fa");
        let gfa_path = directory.join("unitigs.gfa");
        let dot_path = directory.join("graph.dot");
        let overlaps_path = directory.join("overlaps.bin");
        let fastq: String = (0..6)
            .map(|read_id| format!("@read{read_id}\nACGTACGT\n+\nIIIIIIII\n"))
            .collect();
        std::fs::write(&fastq_path, fastq).unwrap();

        for seed in [0, 42] {
            utils::set_seed(Some(seed));
            let mut expected = None;
            for attempt in 0..32 {
                let output = alignment_filtering::AlignmentFilteringOutput {
                    overlaps: fixture_overlaps(attempt % 2 != 0),
                };
                output
                    .serialize_overlaps(overlaps_path.to_str().unwrap())
                    .unwrap();
                let overlap_bytes = std::fs::read(&overlaps_path).unwrap();
                let overlaps: HashMap<(usize, usize), Overlap> =
                    bincode::deserialize(&overlap_bytes).unwrap();
                let original: BTreeMap<_, _> = output.overlaps.iter().collect();
                let decoded: BTreeMap<_, _> = overlaps.iter().collect();
                assert_eq!(
                    bincode::serialize(&original).unwrap(),
                    bincode::serialize(&decoded).unwrap()
                );
                let mut graph = create_overlap_graph::run_create_overlap_graph(&overlaps).unwrap();
                transitive_edge_reduction::reduce_transitive_edges(&mut graph, 0);
                let remaining = &graph.nodes["read5:0-8-"].edges;
                assert_eq!(remaining.len(), 3);
                assert!(remaining.iter().all(|edge| {
                    edge.target_id != "read0:0-8-" && edge.target_id != "read1:0-8-"
                }));
                let snapshot: BTreeMap<_, _> = graph
                    .nodes
                    .iter()
                    .map(|(node_id, node)| {
                        let edges: Vec<_> = node
                            .edges
                            .iter()
                            .map(|edge| (edge.target_id.clone(), edge.edge_len))
                            .collect();
                        (node_id.clone(), edges)
                    })
                    .collect();
                heuristic_simplification::remove_multi_edges(&mut graph);
                heuristic_simplification::remove_short_edges(&mut graph, 0.8);
                bubble_removal::remove_bubbles(&mut graph, 100, 1.1);
                graph.write_dot(&dot_path).unwrap();
                let compressed = compress_graph::compress_unitigs(&graph, &overlaps, &fastq_path);
                compressed.write_fasta(&fasta_path).unwrap();
                assert!(compressed.edges.len() >= 2);
                compressed.write_gfa(gfa_path.to_str().unwrap()).unwrap();
                let result = (
                    snapshot,
                    overlap_bytes,
                    std::fs::read(&fasta_path).unwrap(),
                    std::fs::read(&gfa_path).unwrap(),
                    std::fs::read(&dot_path).unwrap(),
                );
                if let Some(expected) = &expected {
                    assert_eq!(&result, expected, "seed {seed}, attempt {attempt}");
                } else {
                    expected = Some(result);
                }
            }
        }
        utils::set_seed(None);
        std::fs::remove_dir_all(directory).unwrap();
    }
}