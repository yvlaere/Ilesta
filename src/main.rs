mod align_reads;
mod alignment_filtering;
mod bubble_removal;
mod circularise;
mod cli;
mod compress_graph;
mod consensus;
mod create_overlap_graph;
mod heuristic_simplification;
mod rescue;
mod scaffold;
mod segment_dedup;
mod tip_trimming;
mod transitive_edge_reduction;
mod ui;
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
            print_header("Align", &[("reads", args.reads.reads_fq.clone()), ("output", args.output_dir.clone())]);
            ui::plan(1);
            ui::stage("Read processing and alignment");
            align_reads::align_reads(
                &args.reads,
                out_dir,
                align_reads::ReadSelection::TargetCoverage,
            )?;
            print_footer();
        }
        Commands::AlignmentFiltering(args) => {
            utils::set_seed(Some(args.seed));
            let out_dir = Path::new(&args.output_dir);
            std::fs::create_dir_all(out_dir)?;
            print_header("Alignment Filtering", &[("alignments", args.paf.clone()), ("output", args.output_overlaps.clone())]);
            ui::plan(1);
            ui::stage("Overlap filtering");
            utils::timed_step("Filter alignments and write overlaps", || {
                let out = alignment_filtering::run_alignment_filtering(
                    &out_dir.join(&args.paf),
                    &args.overlap,
                    true,
                )?;
                out.serialize_overlaps(&args.output_overlaps)?;
                ui::detail(format!("wrote {} overlaps to {}", ui::count(out.overlaps.len() as u64), args.output_overlaps));
                Ok::<_, Box<dyn Error>>(())
            })?;
            print_footer();
        }
        Commands::Scaffold(args) => {
            let out = Path::new(&args.output);
            let work_dir = out
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            let reads = Path::new(&args.reads_fq);
            print_header("Scaffold", &[("contigs", args.contigs.clone()), ("reads", args.reads_fq.clone()), ("output", args.output.clone())]);
            ui::plan(1 + usize::from(args.polish_rounds > 0));
            let mut graph = compress_graph::CompressedGraph::from_fasta(Path::new(&args.contigs))?;
            ui::stage("Contig bridging");
            let stats = scaffold::scaffold_graph(&mut graph, reads, work_dir, &args.scaffold, args.threads, args.read_type)?;
            ui::detail(format!("bridges applied: {} · contigs circularised: {}", stats.accepted, stats.circularised));
            if args.polish_rounds > 0 {
                let params = consensus::ConsensusParams {
                    rounds: args.polish_rounds,
                    window_len: 500,
                    min_depth: 3,
                    max_depth: 25,
                    threads: args.threads,
                    read_type: args.read_type,
                };
                ui::stage("Consensus polishing");
                let stats = consensus::polish_graph(&mut graph, reads, work_dir, &params)?;
                ui::detail(format!(
                    "{} rounds · {} → {} · {} windows changed",
                    stats.rounds_run,
                    ui::bp(stats.bases_before as u64),
                    ui::bp(stats.bases_after as u64),
                    ui::count(stats.windows_changed as u64)
                ));
            }
            graph.write_fasta(out)?;
            print_footer();
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

/// The command line as typed, quoted so that it can be pasted back into a shell.
fn command_line() -> String {
    let quote = |arg: String| {
        let plain = |c: char| c.is_ascii_alphanumeric() || "_-./:=,+@%".contains(c);
        if !arg.is_empty() && arg.chars().all(plain) { arg } else { format!("'{}'", arg.replace('\'', "'\\''")) }
    };
    let mut args = std::env::args();
    let program = args.next().map(|p| Path::new(&p).file_name().map_or(p.clone(), |n| n.to_string_lossy().into_owned()));
    program.into_iter().chain(args).map(quote).collect::<Vec<_>>().join(" ")
}

/// Run header: version, start time, the exact command line and the main inputs.
fn print_header(command: &str, items: &[(&str, String)]) {
    ui::banner(&format!("Ilesta {} · {command}", env!("CARGO_PKG_VERSION")));
    ui::detail(format!("{:<10}{}", "command", command_line()));
    ui::detail(format!("{:<10}{}", "started", utils::utc_timestamp()));
    for (name, value) in items {
        ui::detail(format!("{name:<10}{value}"));
    }
}

fn print_footer() {
    ui::banner("Done");
    let times = ui::stage_times();
    let width = times.iter().map(|t| t.0.chars().count()).max().unwrap_or(0) + 3;
    for (name, seconds) in &times {
        ui::detail(format!("{name:<width$}{}", ui::duration(*seconds)));
    }
    if times.len() > 1 {
        let total = ui::duration(times.iter().map(|t| t.1).sum());
        ui::detail(format!("{:<width$}{}", "total", ui::bold(&total)));
    }
}

fn edge_count(graph: &OverlapGraph) -> usize {
    graph.nodes.values().map(|n| n.edges.len()).sum()
}

/// Iterative graph cleanup: transitive reduction, multi-edges, short edges, bubbles,
/// isolated nodes and tips.
fn cleanup_graph(graph: &mut OverlapGraph, args: &AssembleArgs) {
    utils::check_synchronization(graph);
    for iteration in 1..=args.cleanup_iterations {
        let before = edge_count(graph);
        transitive_edge_reduction::reduce_transitive_edges(graph, args.fuzz);
        let transitive = before - edge_count(graph);
        let multi_edges = heuristic_simplification::remove_multi_edges(graph);
        let short_edges =
            heuristic_simplification::remove_short_edges(graph, args.short_edge_ratio);
        let before = graph.nodes.len();
        bubble_removal::remove_bubbles(
            graph,
            args.max_bubble_length as usize,
            args.min_support_ratio,
        );
        let bubbles = before - graph.nodes.len();

        // nodes without any edge
        let targets: HashSet<&str> = graph
            .nodes
            .values()
            .flat_map(|n| n.edges.iter().map(|e| e.target_id.as_str()))
            .collect();
        let isolated: HashSet<String> = graph
            .nodes
            .iter()
            .filter(|(id, n)| n.edges.is_empty() && !targets.contains(id.as_str()))
            .map(|(id, _)| id.clone())
            .collect();
        utils::delete_nodes_and_edges(graph, &isolated);

        let before = graph.nodes.len();
        tip_trimming::trim_tips(graph, args.max_tip_len as usize);
        ui::detail(format!(
            "iteration {iteration}: removed {} transitive, {} multi-, {} short edges; {} bubble, {} isolated, {} tip nodes",
            ui::count(transitive as u64),
            ui::count(multi_edges as u64),
            ui::count(short_edges as u64),
            ui::count(bubbles as u64),
            ui::count(isolated.len() as u64),
            ui::count((before - graph.nodes.len()) as u64)
        ));
    }
    utils::check_synchronization(graph);
    ui::detail(format!(
        "cleaned graph: {} nodes · {} directed edges",
        ui::count(graph.nodes.len() as u64),
        ui::count(edge_count(graph) as u64)
    ));
}

/// Run the assembly. Intermediate files are in `<output-dir>/<prefix>.work`, which is removed
/// afterwards unless `--keep-intermediate` is given (then they stay in the output directory).
/// After a failure the directory is kept so the run can be inspected.
fn assemble(args: &AssembleArgs) -> Result<(), Box<dyn Error>> {
    utils::set_seed(Some(args.seed));
    let out_dir = Path::new(&args.output_dir);
    std::fs::create_dir_all(out_dir)?;
    check_exists(&args.reads.reads_fq)?;
    let work_dir = if args.keep_intermediate {
        out_dir.to_path_buf()
    } else {
        out_dir.join(format!("{}.work", args.output_prefix))
    };
    std::fs::create_dir_all(&work_dir)?;

    let result = run_assembly(args, out_dir, &work_dir);
    match &result {
        Ok(()) if !args.keep_intermediate => {
            std::fs::remove_dir_all(&work_dir)?;
        }
        Err(_) if !args.keep_intermediate => {
            eprintln!("Intermediate files were kept in {} for inspection", work_dir.display());
        }
        _ => {}
    }
    result
}

fn run_assembly(args: &AssembleArgs, out_dir: &Path, work_dir: &Path) -> Result<(), Box<dyn Error>> {
    let (threads, read_type) = (args.reads.threads, args.reads.read_type);
    let polishing = args.consensus_rounds > 0;
    let fasta_path = out_dir.join(format!("{}.fa", args.output_prefix));
    let gfa_path = out_dir.join(format!("{}.gfa", args.output_prefix));
    // sequence of the current contigs, used as mapping input between the stages
    let draft_path = work_dir.join(format!("{}.draft.fa", args.output_prefix));

    print_header(
        "assemble",
        &[
            ("reads", args.reads.reads_fq.clone()),
            ("output", format!("{} ({}, {})", out_dir.display(), fasta_path.file_name().unwrap().to_string_lossy(), gfa_path.file_name().unwrap().to_string_lossy())),
            (
                "settings",
                format!("{threads} threads · {} reads · seed {} · target coverage {}x", format!("{read_type:?}").to_lowercase(), args.seed, args.reads.target_coverage),
            ),
        ],
    );
    ui::plan(5 + usize::from(args.rescue_enabled()) + usize::from(polishing));

    // ---- 1. reads ---------------------------------------------------------
    ui::stage("Read processing and alignment");
    let reads = align_reads::align_reads(
        &args.reads,
        work_dir,
        align_reads::ReadSelection::TargetCoverage,
    )?;

    // ---- 2. overlaps: provided, or filtered from the alignments ------------
    ui::stage("Overlap filtering");
    let overlaps_path = match &args.overlaps {
        Some(path) => {
            ui::detail(format!("using precomputed overlaps from {path}"));
            PathBuf::from(path)
        }
        None => {
            let path = work_dir.join(format!("{}.overlaps.bin", args.output_prefix));
            utils::timed_step("Filter read overlaps", || {
                let out = alignment_filtering::run_alignment_filtering(
                    &work_dir.join(&args.reads.paf),
                    &args.overlap,
                    true,
                )?;
                out.serialize_overlaps(path.to_str().ok_or("invalid output path")?)?;
                Ok::<_, Box<dyn Error>>(())
            })?;
            path
        }
    };

    // ---- 3. graph ---------------------------------------------------------
    ui::stage("Graph construction");
    let (mut overlaps, mut graph): (HashMap<(usize, usize), Overlap>, OverlapGraph) =
        utils::timed_step("Create overlap graph", || {
            let overlaps: HashMap<(usize, usize), Overlap> =
                bincode::deserialize_from(BufReader::new(File::open(&overlaps_path)?))?;
            let graph = create_overlap_graph::run_create_overlap_graph(&overlaps)?;
            Ok::<_, Box<dyn Error>>((overlaps, graph))
        })?;
    utils::timed_step("Clean overlap graph", || {
        cleanup_graph(&mut graph, args);
        Ok::<_, Box<dyn Error>>(())
    })?;
    if args.write_dot {
        utils::timed_step("Write overlap graph visualization", || {
            graph.write_dot(out_dir.join("graph.dot"))
        })?;
    }

    let mut assembly_reads = reads.primary_reads.clone();
    let mut compressed = utils::timed_step("Compress graph into contigs", || {
        Ok::<_, Box<dyn Error>>(compress_graph::compress_contigs(
            &graph,
            &overlaps,
            &assembly_reads,
        ))
    })?;
    report_contigs("contigs", &compressed);

    // ---- 4. rescue --------------------------------------------------------
    let mut rescued_nodes = HashSet::new();
    if args.rescue_enabled() {
        ui::stage("Rescue of poorly represented reads");
        compressed.write_fasta(&draft_path)?;
        let rescue = rescue::run_rescue_stage(
            &reads.filtered_reads,
            &draft_path,
            work_dir,
            &args.reads,
            &cli::OverlapOpts { min_overlap_length: args.rescue_min_overlap_length, ..args.overlap.clone() },
            &mut graph,
            &mut overlaps,
        )?;
        if let Some(rescue) = rescue {
            let combined = work_dir.join(format!("{}.combined.fq", args.output_prefix));
            let mut file = File::create(&combined)?;
            std::io::copy(&mut File::open(&reads.primary_reads)?, &mut file)?;
            std::io::copy(&mut File::open(&rescue.reads)?, &mut file)?;
            assembly_reads = combined;
            rescued_nodes = rescue.node_ids;
            compressed = compress_graph::compress_contigs(&graph, &overlaps, &assembly_reads);
            report_contigs("contigs after rescue", &compressed);
        }
    }

    // ---- 5. refinement ----------------------------------------------------
    ui::stage("Contig refinement");
    compressed.write_fasta(&draft_path)?;
    let removed = utils::timed_step("Remove contained contigs", || {
        segment_dedup::remove_contained_segments(&mut compressed, &draft_path, work_dir, threads)
    })?;
    ui::detail(format!("contained contigs removed: {removed}"));
    report_contigs("contigs after dedup", &compressed);

    if args.scaffold.rounds > 0 {
        let stats = scaffold::scaffold_graph(&mut compressed, &assembly_reads, work_dir, &args.scaffold, threads, read_type)?;
        ui::detail(format!("bridges applied: {} · contigs circularised: {}", stats.accepted, stats.circularised));
        report_contigs("contigs after bridging", &compressed);
    }
    let circularised = utils::timed_step("Circularise contigs with overlapping ends", || {
        circularise::circularise_self_overlaps(&mut compressed, work_dir, threads)
    })?;
    ui::detail(format!("contigs circularised by overlapping ends: {circularised}"));

    // ---- 6. consensus -----------------------------------------------------
    if polishing {
        ui::stage("Consensus polishing");
        let params = consensus::ConsensusParams {
            rounds: args.consensus_rounds,
            window_len: args.consensus_window_len,
            min_depth: args.consensus_min_depth,
            max_depth: args.consensus_max_depth,
            threads,
            read_type,
        };
        let stats = consensus::polish_graph(&mut compressed, &assembly_reads, work_dir, &params)?;
        ui::detail(format!(
            "{} rounds · {} → {} · {} windows changed, {} lacked support",
            stats.rounds_run,
            ui::bp(stats.bases_before as u64),
            ui::bp(stats.bases_after as u64),
            ui::count(stats.windows_changed as u64),
            ui::count(stats.windows_unsupported as u64)
        ));

        // Second pass on polished sequences: copies that were too diverged to be recognised in
        // the raw contigs (and those created by bridging) align cleanly now.
        compressed.write_fasta(&draft_path)?;
        let removed = utils::timed_step("Remove contained contigs (polished)", || {
            segment_dedup::remove_contained_segments(&mut compressed, &draft_path, work_dir, threads)
        })?;
        ui::detail(format!("contained contigs removed after polishing: {removed}"));
    }

    // ---- 7. output --------------------------------------------------------
    ui::stage("Output");
    // The `a` lines of the GFA place reads on the contig sequence. Polishing changes that
    // sequence, so they are only written when it is off.
    compressed.write_fasta(&fasta_path)?;
    compressed.write_gfa(gfa_path.to_str().ok_or("invalid output path")?, !polishing)?;
    if !rescued_nodes.is_empty() {
        let rescued = compressed
            .contigs
            .iter()
            .filter(|u| u.members.iter().any(|m| rescued_nodes.contains(&m.node_id)))
            .count();
        ui::detail(format!("contigs containing rescued reads: {rescued}"));
    }
    ui::detail(format!("wrote {} and {}", fasta_path.display(), gfa_path.display()));
    if assembly_reads != reads.primary_reads {
        std::fs::remove_file(&assembly_reads)?;
    }
    print_assembly_summary(&compressed);
    print_footer();
    Ok(())
}

/// One-line contig overview for a stage: count, total length and the largest contigs.
fn report_contigs(stage: &str, graph: &compress_graph::CompressedGraph) {
    let mut lens: Vec<(usize, usize, char)> = graph
        .contigs
        .iter()
        .map(|u| (u.fasta_seq.as_ref().map_or(0, String::len), u.id, u.topology))
        .collect();
    lens.sort_unstable_by(|a, b| b.cmp(a));
    let total: u64 = lens.iter().map(|l| l.0 as u64).sum();
    let largest: Vec<String> = lens
        .iter()
        .take(3)
        .map(|&(len, _, topology)| format!("{}{}", ui::bp(len as u64), if topology == 'c' { " (circular)" } else { "" }))
        .collect();
    let largest = if largest.is_empty() { String::new() } else { format!(" · largest {}", largest.join(", ")) };
    ui::detail(format!("{stage}: {} · {}{largest}", ui::count(lens.len() as u64), ui::bp(total)));
}

/// Final overview: totals and one line per contig.
fn print_assembly_summary(graph: &compress_graph::CompressedGraph) {
    ui::banner("Assembly");
    let mut contigs: Vec<(usize, &compress_graph::contig)> = graph
        .contigs
        .iter()
        .map(|u| (u.fasta_seq.as_ref().map_or(0, String::len), u))
        .collect();
    contigs.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.id.cmp(&b.1.id)));
    let total: u64 = contigs.iter().map(|c| c.0 as u64).sum();
    let circular = contigs.iter().filter(|c| c.1.topology == 'c').count();
    let gc = |seq: &str| {
        let gc = seq.bytes().filter(|b| matches!(b.to_ascii_uppercase(), b'G' | b'C')).count();
        100.0 * gc as f64 / seq.len().max(1) as f64
    };
    ui::detail(format!(
        "{} contigs ({circular} circular) · {} · N50 {} · largest {}",
        contigs.len(),
        ui::bp(total),
        ui::bp(ui::n50(contigs.iter().map(|c| c.0 as u64))),
        ui::bp(contigs.first().map_or(0, |c| c.0 as u64))
    ));
    for (len, contig) in contigs.iter().take(10) {
        ui::detail(format!(
            "contig_{:<5}{:>12}  {:<9}GC {:.1}%",
            contig.id,
            ui::count(*len as u64),
            if contig.topology == 'c' { "circular" } else { "linear" },
            gc(contig.fasta_seq.as_deref().unwrap_or(""))
        ));
    }
    if contigs.len() > 10 {
        ui::detail(format!("… and {} more", contigs.len() - 10));
    }
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
        let fasta_path = directory.join("contigs.fa");
        let gfa_path = directory.join("contigs.gfa");
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
                let compressed = compress_graph::compress_contigs(&graph, &overlaps, &fastq_path);
                compressed.write_fasta(&fasta_path).unwrap();
                assert!(compressed.edges.len() >= 2);
                compressed.write_gfa(gfa_path.to_str().unwrap(), true).unwrap();
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
