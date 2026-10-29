use crate::compress_graph::CompressedGraph;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Minimum fraction of matching bases (matches / alignment block) of an alignment that counts as
/// containment. This is measured on base-level alignments (`-c`) of contigs that are still
/// unpolished, so two copies of the same sequence typically agree on ~90-95% of the bases.
const MIN_IDENTITY: f64 = 0.75;
/// Minimum fraction of the query that needs to be covered
const MIN_QUERY_COVERAGE: f64 = 0.90;
/// minimap2 options of the all-vs-all contig alignment. `-c` is essential: without it minimap2
/// reports the matches of the seed chain (a small fraction of the bases at ONT error rates)
const MINIMAP2_ARGS: &[&str] = &["-x", "asm20", "-DP", "-c"];
/// Removed contigs of at least this length are reported (replicons, not rescue fragments).
const REPORT_MIN_LEN: u64 = 20_000;

/// Remove contigs that are (>90%) contained in a longer one. `fasta` holds the current contigs
/// and is only used as mapping input. Returns the number of removed contigs.
pub fn remove_contained_segments(
    assembly: &mut CompressedGraph,
    fasta: &Path,
    work_dir: &Path,
    threads: usize,
) -> io::Result<usize> {
    if assembly.unitigs.is_empty() {
        return Ok(0);
    }
    let paf = work_dir.join("assembly_containment.paf");
    crate::utils::run_minimap2(MINIMAP2_ARGS, threads, fasta, fasta, &paf)?;
    let contained = contained_unitig_ids(BufReader::new(File::open(paf)?))?;
    let before = assembly.unitigs.len();
    assembly.unitigs.retain(|u| !contained.contains(&u.id));
    let retained: HashSet<_> = assembly.unitigs.iter().map(|u| u.id).collect();
    assembly.edges.retain(|e| retained.contains(&e.from) && retained.contains(&e.to));
    Ok(before - assembly.unitigs.len())
}

fn contained_unitig_ids(reader: impl BufRead) -> io::Result<HashSet<usize>> {
    let mut alignments: HashMap<(String, String), QueryAlignments> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 12 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid PAF record during assembly segment containment filtering",
            ));
        }

        // Only real alignments count: with `-P` minimap2 also reports long, loose chains between
        // repeat copies (a few percent matching bases) that would cover most of a plasmid.
        let (matches, block) = (parse_paf_number(fields[9])?, parse_paf_number(fields[10])?);
        if block == 0 || (matches as f64) < MIN_IDENTITY * block as f64 {
            continue;
        }

        let query_length = parse_paf_number(fields[1])?;
        let query_start = parse_paf_number(fields[2])?;
        let query_end = parse_paf_number(fields[3])?;
        let target_length = parse_paf_number(fields[6])?;

        if fields[0] == fields[5]
            || query_length == 0
            || query_start > query_end
            || query_end > query_length
        {
            continue;
        }

        let key = (fields[0].to_string(), fields[5].to_string());
        let alignment = alignments.entry(key).or_insert_with(|| QueryAlignments {
            query_length,
            target_length,
            intervals: Vec::new(),
        });
        alignment.intervals.push((query_start, query_end));
    }

    let mut contained = HashSet::new();
    for ((query_name, target_name), mut alignment) in alignments {
        if alignment.target_length < alignment.query_length
            || (alignment.target_length == alignment.query_length && target_name >= query_name)
            || covered_query_fraction(&mut alignment.intervals, alignment.query_length) <= MIN_QUERY_COVERAGE
        {
            continue;
        }
        if alignment.query_length >= REPORT_MIN_LEN {
            let covered = covered_query_fraction(&mut alignment.intervals, alignment.query_length);
            crate::ui::detail(format!(
                "removing {query_name} ({}): {:.1}% covered by {target_name} ({})",
                crate::ui::bp(alignment.query_length as u64),
                100.0 * covered,
                crate::ui::bp(alignment.target_length as u64)
            ));
        }

        if let Some(id) = query_name
            .strip_prefix("contig_")
            .and_then(|name| name.parse().ok())
        {
            contained.insert(id);
        }
    }
    Ok(contained)
}

struct QueryAlignments {
    query_length: u64,
    target_length: u64,
    intervals: Vec<(u64, u64)>,
}

fn covered_query_fraction(intervals: &mut [(u64, u64)], query_length: u64) -> f64 {
    intervals.sort_unstable();
    let mut covered = 0;
    let mut current: Option<(u64, u64)> = None;

    for &(start, end) in intervals.iter() {
        match current {
            Some((current_start, current_end)) if start <= current_end => {
                current = Some((current_start, current_end.max(end)));
            }
            Some((current_start, current_end)) => {
                covered += current_end - current_start;
                current = Some((start, end));
            }
            None => current = Some((start, end)),
        }
    }
    if let Some((start, end)) = current {
        covered += end - start;
    }

    covered as f64 / query_length as f64
}

fn parse_paf_number(value: &str) -> io::Result<u64> {
    value.parse().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid numeric field in PAF record",
        )
    })
}

#[cfg(test)]
mod tests {
    use super::contained_unitig_ids;
    use std::io::Cursor;

    #[test]
    fn removes_only_shorter_segments_passing_query_coverage_threshold() {
        let paf = concat!(
            "contig_1\t100\t0\t96\t+\tcontig_2\t200\t20\t116\t96\t96\t60\n",
            "contig_3\t100\t0\t95\t+\tcontig_4\t200\t20\t115\t95\t95\t60\n",
            "contig_5\t100\t0\t99\t+\tcontig_6\t200\t20\t119\t99\t100\t60\n",
            "contig_7\t200\t0\t198\t+\tcontig_8\t100\t0\t98\t198\t198\t60\n",
            "contig_9\t100\t0\t99\t+\tcontig_9\t100\t0\t99\t99\t99\t60\n",
            "contig_10\t100\t0\t99\t+\tcontig_11\t100\t0\t99\t99\t99\t60\n",
            "contig_11\t100\t0\t99\t+\tcontig_10\t100\t0\t99\t99\t99\t60\n",
            // loose chain between repeat copies: 99% of the query is spanned, 2% of it matches
            "contig_12\t100\t0\t99\t+\tcontig_13\t200\t20\t119\t2\t99\t0\n",
        );

        let contained = contained_unitig_ids(Cursor::new(paf)).unwrap();

        assert_eq!(contained, [1, 5, 11].into_iter().collect());
    }

    #[test]
    fn combines_non_overlapping_query_alignments_but_not_different_targets() {
        let paf = concat!(
            "contig_20\t100\t0\t50\t+\tcontig_21\t200\t0\t50\t50\t50\t60\n",
            "contig_20\t100\t50\t100\t-\tcontig_21\t200\t50\t100\t50\t50\t60\n",
            "contig_22\t100\t0\t50\t+\tcontig_23\t200\t0\t50\t50\t50\t60\n",
            "contig_22\t100\t50\t100\t+\tcontig_24\t200\t0\t50\t50\t50\t60\n",
        );

        let contained = contained_unitig_ids(Cursor::new(paf)).unwrap();

        assert_eq!(contained, [20].into_iter().collect());
    }
}

    /// A 6 kb copy of a longer contig with ~4% substitutions and indels (what two unpolished
    /// contigs of the same region look like) must be recognised as contained; an unrelated
    /// contig must not.
    #[test]
    #[ignore = "requires minimap2 on PATH"]
    fn removes_a_diverged_copy_but_not_an_unrelated_contig() {
        let mut state = 12345u64;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        let random = |n: usize, next: &mut dyn FnMut() -> usize| -> Vec<u8> {
            (0..n).map(|_| b"ACGT"[next() % 4]).collect()
        };
        let big = random(30_000, &mut next);
        let mut copy = Vec::new();
        for &b in &big[10_000..16_000] {
            match next() % 100 {
                0..=1 => {}                                   // deletion
                2 => { copy.push(b); copy.push(b"ACGT"[next() % 4]); } // insertion
                3 => copy.push(b"ACGT"[next() % 4]),          // substitution
                _ => copy.push(b),
            }
        }
        let other = random(8_000, &mut next);

        let dir = std::env::temp_dir().join(format!("ilesta-dedup-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("contigs.fa");
        let mut text = String::new();
        for (i, seq) in [&big, &copy, &other].iter().enumerate() {
            text += &format!(">contig_{i} topology=l\n{}\n", String::from_utf8_lossy(seq));
        }
        std::fs::write(&fasta, text).unwrap();
        let mut graph = crate::compress_graph::CompressedGraph::from_fasta(&fasta).unwrap();
        let removed = remove_contained_segments(&mut graph, &fasta, &dir, 1).unwrap();
        let _ = std::fs::remove_dir_all(&dir);

        let kept: Vec<usize> = graph.unitigs.iter().map(|u| u.id).collect();
        assert_eq!((removed, kept), (1, vec![0, 2]));
    }