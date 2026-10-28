//! Circularisation of linear unitigs whose two ends overlap.
//!
//! When the overlap graph does not close a cycle (low coverage at the junction, a contained
//! read removed at the wrong place), the unitig of a circular replicon is a linear walk that
//! runs past its own start: the last tens of kb repeat the first tens of kb. Such a unitig is
//! recognised by aligning its head to its tail; if the head aligns to the tail end-to-end
//! (same strand, both ends reached), the duplicated copy is cut off and the unitig is marked
//! circular

use crate::compress_graph::CompressedGraph;
use std::fs;
use std::io;
use std::path::Path;

/// Smallest overlap that is taken as evidence for a ring. Shorter shared ends are typically
/// repeats (IS elements, rRNA operons) at the boundary of a contig, not a closed replicon.
const MIN_OVERLAP: usize = 10_000;
/// Largest piece of each end that is compared.
const MAX_WINDOW: usize = 300_000;
/// How far the alignment may stop short of the outermost bases of the unitig.
const END_SLACK: usize = 1_000;
/// Minimum matches / alignment block of the end-to-end alignment.
const MIN_IDENTITY: f64 = 0.8;

/// From the PAF of `head` (query) against `tail` (target), where `head` is the first
/// `window` bases and `tail` the last `window` bases of a unitig of `len` bases, return the
/// range `[start, end)` of the unitig that forms the ring, if its ends overlap.
fn ring_range(paf: &str, len: usize, window: usize) -> Option<(usize, usize)> {
    let mut best: Option<(usize, usize, usize)> = None; // (overlap, ring start, ring end)
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 || f[4] != "+" {
            continue;
        }
        let num = |i: usize| f[i].parse::<usize>().ok();
        let (qs, qe, tlen, ts, te, matches, block) =
            (num(2)?, num(3)?, num(6)?, num(7)?, num(8)?, num(9)?, num(10)?);
        if block == 0 || (matches as f64) < MIN_IDENTITY * block as f64 {
            continue;
        }
        let overlap = qe.saturating_sub(qs);
        // the head must start at the start of the unitig, the tail must end at its end
        if qs > END_SLACK || te + END_SLACK < tlen || overlap < MIN_OVERLAP || overlap * 2 > len {
            continue;
        }
        let ring_end = len - window + ts;
        if ring_end <= qs || ring_end - qs < len / 2 {
            continue;
        }
        if best.is_none_or(|b| overlap > b.0) {
            best = Some((overlap, qs, ring_end));
        }
    }
    best.map(|(_, s, e)| (s, e))
}

/// Circularise the linear unitigs of `graph` whose ends overlap. Returns how many were changed.
pub fn circularise_self_overlaps(
    graph: &mut CompressedGraph,
    work_dir: &Path,
    threads: usize,
) -> io::Result<usize> {
    let mut changed = 0;
    for u in graph.unitigs.iter_mut() {
        let Some(seq) = u.fasta_seq.as_ref() else { continue };
        let len = seq.len();
        if u.topology == 'c' || len < 2 * MIN_OVERLAP {
            continue;
        }
        let window = (len / 2).min(MAX_WINDOW);
        let head = work_dir.join("circ_head.fa");
        let tail = work_dir.join("circ_tail.fa");
        let paf = work_dir.join("circ.paf");
        fs::write(&head, format!(">head\n{}\n", &seq[..window]))?;
        fs::write(&tail, format!(">tail\n{}\n", &seq[len - window..]))?;
        crate::utils::run_minimap2(
            &["-x", "asm20", "-c", "--secondary=no"],
            threads,
            &tail,
            &head,
            &paf,
        )?;
        let text = fs::read_to_string(&paf)?;
        for p in [&head, &tail, &paf] {
            let _ = fs::remove_file(p);
        }
        if let Some((start, end)) = ring_range(&text, len, window) {
            crate::ui::detail(format!(
                "circularising unitig_{} ({}): its ends overlap by {}",
                u.id,
                crate::ui::bp(len as u64),
                crate::ui::bp((len - (end - start)) as u64)
            ));
            u.fasta_seq = Some(seq[start..end].to_string());
            u.topology = 'c';
            // the read tiling no longer describes the trimmed sequence
            u.members = Vec::new();
            changed += 1;
        }
    }
    Ok(changed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_the_ring_from_an_end_to_end_alignment() {
        // 100 kb unitig, window 50 kb: head[0..30000] == tail[20000..50000] (tail starts at 50000)
        let paf = "head\t50000\t0\t30000\t+\ttail\t50000\t20000\t50000\t29000\t30000\t60\n";
        assert_eq!(ring_range(paf, 100_000, 50_000), Some((0, 70_000)));
        // an alignment that does not reach the end of the tail is a repeat, not a ring
        let paf = "head\t50000\t0\t30000\t+\ttail\t50000\t10000\t40000\t29000\t30000\t60\n";
        assert_eq!(ring_range(paf, 100_000, 50_000), None);
        // opposite strand, short overlap and low identity are rejected
        let rev = "head\t50000\t0\t30000\t-\ttail\t50000\t20000\t50000\t29000\t30000\t60\n";
        let short = "head\t50000\t0\t5000\t+\ttail\t50000\t45000\t50000\t4900\t5000\t60\n";
        let loose = "head\t50000\t0\t30000\t+\ttail\t50000\t20000\t50000\t9000\t30000\t60\n";
        for paf in [rev, short, loose] {
            assert_eq!(ring_range(paf, 100_000, 50_000), None);
        }
    }

    /// A 60 kb ring that was walked 20 kb too far (the last 20 kb repeat the first 20 kb with
    /// ~3% noise) is cut back to 60 kb and marked circular; a unitig without overlapping ends
    /// is left alone.
    #[test]
    #[ignore = "requires minimap2 on PATH"]
    fn trims_and_marks_an_overrun_unitig() {
        let mut state = 99u64;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        let ring: Vec<u8> = (0..60_000).map(|_| b"ACGT"[next() % 4]).collect();
        let mut dup = Vec::new();
        for &b in &ring[..20_000] {
            match next() % 100 {
                0 => {}
                1 => {
                    dup.push(b);
                    dup.push(b"ACGT"[next() % 4]);
                }
                2 => dup.push(b"ACGT"[next() % 4]),
                _ => dup.push(b),
            }
        }
        let mut overrun = ring.clone();
        overrun.extend_from_slice(&dup);
        let unrelated: Vec<u8> = (0..50_000).map(|_| b"ACGT"[next() % 4]).collect();

        let dir = std::env::temp_dir().join(format!("ilesta-circ-test-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("c.fa");
        let text = format!(
            ">unitig_0 topology=l\n{}\n>unitig_1 topology=l\n{}\n",
            String::from_utf8_lossy(&overrun),
            String::from_utf8_lossy(&unrelated)
        );
        fs::write(&fasta, text).unwrap();
        let mut g = CompressedGraph::from_fasta(&fasta).unwrap();
        let n = circularise_self_overlaps(&mut g, &dir, 1).unwrap();
        let _ = fs::remove_dir_all(&dir);
        assert_eq!(n, 1);
        assert_eq!(g.unitigs[0].topology, 'c');
        let got = g.unitigs[0].fasta_seq.as_ref().unwrap().len();
        assert!(got.abs_diff(60_000) <= 300, "ring length {got}");
        assert_eq!(g.unitigs[1].topology, 'l');
        assert_eq!(g.unitigs[1].fasta_seq.as_ref().unwrap().len(), 50_000);
    }
}