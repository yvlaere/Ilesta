//! Consensus (self-polishing) step
//!
//! Method:
//! 1. Map the assembly reads back to the contigs with `minimap2 -c`.
//! 2. Cut every contig into windows (default ~500 bp) and use each CIGAR to
//!      split every read alignment at the window borders.
//!   3. Per window, build a partial-order alignment (POA) graph: the draft
//!      window is the backbone, each read segment is aligned to the graph
//!      and merged in.
//!   4. The heaviest path through the POA graph is the new window sequence.
//!   5. The windows are concatenated; the whole procedure can be repeated
//!
//! Circular contigs (`topology == 'c'`) are mapped with a flank of the
//! opposite end on each side so that reads crossing the origin contribute to
//! the windows on both sides of it, and the windows are tiled around the ring.
//!

use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::cli::ReadType;
use crate::compress_graph::CompressedGraph;

// ---------------------------------------------------------------------------
// Parameters / statistics
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct ConsensusParams {
    /// Number of polishing rounds (0 disables the step).
    pub rounds: u32,
    /// Target window length in bp.
    pub window_len: usize,
    /// Minimum number of read segments for a window to be re-called;
    /// windows with less support keep the draft sequence.
    pub min_depth: usize,
    /// Maximum number of read segments per window (best identity first).
    pub max_depth: usize,
    pub threads: usize,
    pub read_type: ReadType,
}

#[derive(Debug, Default, Clone)]
pub struct ConsensusStats {
    pub rounds_run: u32,
    pub windows_total: usize,
    pub windows_changed: usize,
    pub windows_unsupported: usize,
    pub bases_before: usize,
    pub bases_after: usize,
}

// POA scoring (linear gaps).
const MATCH: i32 = 5;
const MISMATCH: i32 = -4;
const GAP: i32 = -8;
/// Slack (in backbone positions) around the window span of a read segment.
const SLACK: u32 = 24;
/// Minimum half-width of the DP band around the interpolated read position. The band is widened
/// by the length difference between the read segment and its backbone span, so that indels of
/// up to a window length in the draft can be repaired.
const BAND: usize = 24;
/// Upper limit for the half-width of the DP band.
const MAX_BAND: usize = 1500;
/// Minimum identity of a whole read alignment to be used.
const MIN_IDENTITY: f64 = 0.6;
/// A path node is kept only if at least this fraction of the local median
/// read depth supports it (majority rule).
const MIN_SUPPORT: f64 = 0.5;
/// Flank (bp) added to each side of circular contigs while mapping.
const CIRC_FLANK: usize = 2000;
/// Context added on both sides of a window. Reads are cut at the extended borders and only the
/// consensus of the core is kept: the bases next to a cut are the least reliable (a read that
/// ends there has no context to resolve homopolymers/indels), so without a margin the number
/// of avoidable errors grows with the number of window borders.
const WINDOW_MARGIN: usize = 50;
/// Contigs shorter than this are left untouched.
const MIN_CONTIG_LEN: usize = 100;

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Polish all contigs of `graph` in place.
///
/// `reads_path` is the FASTQ(.gz) that was used for assembly (it is only read,
/// never modified).  Intermediate files are written to `work_dir` and removed.
pub fn polish_graph(
    graph: &mut CompressedGraph,
    reads_path: &Path,
    work_dir: &Path,
    params: &ConsensusParams,
) -> Result<ConsensusStats, Box<dyn Error>> {
    let mut stats = ConsensusStats::default();
    if params.rounds == 0 || graph.unitigs.is_empty() {
        return Ok(stats);
    }

    let mut contigs: Vec<Contig> = graph
        .unitigs
        .iter()
        .map(|u| Contig {
            seq: u
                .fasta_seq
                .as_deref()
                .unwrap_or("")
                .bytes()
                .map(|b| b.to_ascii_uppercase())
                .collect(),
            circular: u.topology == 'c',
        })
        .collect();
    stats.bases_before = contigs.iter().map(|c| c.seq.len()).sum();

    for round in 1..=params.rounds {
        let round_stats = crate::utils::timed_step(
            &format!("Consensus round {round}/{}", params.rounds),
            || polish_round(&mut contigs, reads_path, work_dir, params, round),
        )?;
        stats.rounds_run += 1;
        stats.windows_total = round_stats.windows_total;
        stats.windows_changed += round_stats.windows_changed;
        stats.windows_unsupported = round_stats.windows_unsupported;
        crate::ui::detail(format!(
            "round {round}: {} of {} windows changed, {} without enough read support",
            crate::ui::count(round_stats.windows_changed as u64),
            crate::ui::count(round_stats.windows_total as u64),
            crate::ui::count(round_stats.windows_unsupported as u64)
        ));
        if round_stats.windows_changed == 0 {
            break; // converged
        }
    }
    stats.bases_after = contigs.iter().map(|c| c.seq.len()).sum();

    for (unitig, contig) in graph.unitigs.iter_mut().zip(contigs.into_iter()) {
        unitig.fasta_seq = Some(String::from_utf8(contig.seq)?);
    }

    Ok(stats)
}

// ---------------------------------------------------------------------------
// One polishing round
// ---------------------------------------------------------------------------

struct Contig {
    seq: Vec<u8>,
    circular: bool,
}

/// Window layout of one contig.
struct Layout {
    len: usize,
    flank: usize,
    circular: bool,
    /// Core borders in contig coordinates: bounds[0]=0, bounds[K]=len.
    bounds: Vec<usize>,
    /// Extended interval of every window (contig coordinates).
    xs: Vec<usize>,
    xe: Vec<usize>,
    /// Border positions (all `xs`, `xe`) in coordinates of the (flanked) mapping reference.
    ext_bounds: Vec<usize>,
}

impl Layout {
    fn new(len: usize, circular: bool, window_len: usize, margin: usize) -> Layout {
        let flank = if circular { CIRC_FLANK.min(len / 2) } else { 0 };
        let k = ((len as f64 / window_len as f64).round() as usize).max(1);
        let bounds: Vec<usize> = (0..=k).map(|i| i * len / k).collect();
        // an extended interval may overlap its neighbours but never reach past them
        let min_core = bounds.windows(2).map(|w| w[1] - w[0]).min().unwrap_or(len);
        let margin = margin.min(min_core.saturating_sub(2) / 2);
        let xs: Vec<usize> = (0..k).map(|i| bounds[i].saturating_sub(margin)).collect();
        let xe: Vec<usize> = (0..k).map(|i| (bounds[i + 1] + margin).min(len)).collect();
        let ext_len = len + 2 * flank;
        let mut ext_bounds = Vec::new();
        let shifts: &[isize] = if circular { &[-1, 0, 1] } else { &[0] };
        for &t in shifts {
            for &b in xs.iter().chain(&xe).chain(&[0, len]) {
                let x = flank as isize + b as isize + t * len as isize;
                if x > 0 && (x as usize) < ext_len {
                    ext_bounds.push(x as usize);
                }
            }
        }
        ext_bounds.sort_unstable();
        ext_bounds.dedup();
        Layout {
            len,
            flank,
            circular,
            bounds,
            xs,
            xe,
            ext_bounds,
        }
    }

    fn windows(&self) -> usize {
        self.bounds.len() - 1
    }

    /// Contig coordinates of an ext-coordinate interval (not crossing a border).
    fn to_contig(&self, ts: usize, te: usize) -> Option<(usize, usize)> {
        let (mut ts, mut te) = (ts, te);
        if self.circular {
            if ts < self.flank {
                ts += self.len;
                te += self.len;
            } else if ts >= self.flank + self.len {
                ts -= self.len;
                te -= self.len;
            }
        }
        if ts < self.flank {
            return None;
        }
        let y = ts - self.flank;
        let y_end = te - self.flank;
        if y >= self.len || y_end > self.len {
            return None;
        }
        Some((y, y_end))
    }

    /// Windows (as a half-open index range) whose extended interval contains `[y, y_end)`.
    fn covering(&self, y: usize, y_end: usize) -> std::ops::Range<usize> {
        let last = self.xs.partition_point(|&x| x <= y);
        let first = self.xe.partition_point(|&x| x < y_end);
        first..last.max(first)
    }
}

#[derive(Clone)]
struct Seg {
    read: u32,
    rev: bool,
    /// Query interval in the orientation of the alignment.
    qs: u32,
    qe: u32,
    qlen: u32,
    /// Offsets inside the window.
    b: u32,
    e: u32,
    ident: f32,
}

fn polish_round(
    contigs: &mut [Contig],
    reads_path: &Path,
    work_dir: &Path,
    p: &ConsensusParams,
    round: u32,
) -> Result<ConsensusStats, Box<dyn Error>> {
    let threads = p.threads.max(1);

    // ---- 1. write the mapping reference ---------------------------------
    let ref_path = work_dir.join(format!("consensus_ref_{round}.fa"));
    let paf_path = work_dir.join(format!("consensus_aln_{round}.paf"));
    let mut layouts: Vec<Option<Layout>> = Vec::new();
    {
        let mut f = std::io::BufWriter::new(std::fs::File::create(&ref_path)?);
        for (i, c) in contigs.iter().enumerate() {
            if c.seq.len() < MIN_CONTIG_LEN {
                layouts.push(None);
                continue;
            }
            let lay = Layout::new(c.seq.len(), c.circular, p.window_len.max(50), WINDOW_MARGIN);
            writeln!(f, ">c{i}")?;
            if lay.flank > 0 {
                f.write_all(&c.seq[c.seq.len() - lay.flank..])?;
                f.write_all(&c.seq)?;
                f.write_all(&c.seq[..lay.flank])?;
            } else {
                f.write_all(&c.seq)?;
            }
            f.write_all(b"\n")?;
            layouts.push(Some(lay));
        }
        f.flush()?;
    }

    // ---- 2. map reads -----------------------------------------------------
    let mut args: Vec<&str> = p.read_type.mapping_args().to_vec();
    args.extend(["-c", "--secondary=no", "-v", "1"]);
    crate::utils::run_minimap2(&args, threads, &ref_path, reads_path, &paf_path)?;

    // ---- 3. split alignments into window segments ------------------------
    let mut windows: Vec<Vec<Vec<Seg>>> = layouts
        .iter()
        .map(|l| match l {
            Some(l) => vec![Vec::new(); l.windows()],
            None => Vec::new(),
        })
        .collect();
    let mut names: HashMap<String, u32> = HashMap::new();
    {
        let reader = BufReader::new(std::fs::File::open(&paf_path)?);
        for line in reader.lines() {
            let line = line?;
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 12 {
                continue;
            }
            let (Ok(qlen), Ok(qstart), Ok(qend), Ok(tstart), Ok(matches), Ok(block)) = (
                f[1].parse::<usize>(),
                f[2].parse::<usize>(),
                f[3].parse::<usize>(),
                f[7].parse::<usize>(),
                f[9].parse::<usize>(),
                f[10].parse::<usize>(),
            ) else {
                continue;
            };
            let Some(ci) = f[5].strip_prefix('c').and_then(|s| s.parse::<usize>().ok()) else {
                continue;
            };
            let Some(Some(lay)) = layouts.get(ci) else {
                continue;
            };
            let mut cigar = None;
            let mut primary = true;
            for tag in &f[12..] {
                if let Some(c) = tag.strip_prefix("cg:Z:") {
                    cigar = Some(c);
                } else if let Some(t) = tag.strip_prefix("tp:A:") {
                    primary = t == "P";
                }
            }
            let (Some(cigar), true) = (cigar, primary) else {
                continue;
            };
            if block == 0 || (matches as f64) / (block as f64) < MIN_IDENTITY {
                continue;
            }
            let rev = f[4] == "-";
            let q0 = if rev { qlen - qend } else { qstart };
            let ident = matches as f32 / block as f32;
            let mut read_id: Option<u32> = None;
            let wins = &mut windows[ci];
            // pieces between two borders; a window collects the pieces of its extended interval
            let mut acc: Vec<(usize, usize, usize, usize, usize)> = Vec::new(); // (k, qs, qe, y, y_end)
            split_cigar(cigar, tstart, q0, &lay.ext_bounds, |ts, te, qs, qe| {
                let Some((y, y_end)) = lay.to_contig(ts, te) else {
                    return;
                };
                for k in lay.covering(y, y_end) {
                    match acc.iter_mut().rev().find(|a| a.0 == k && a.4 == y) {
                        Some(a) => {
                            a.2 = qe;
                            a.4 = y_end;
                        }
                        None => acc.push((k, qs, qe, y, y_end)),
                    }
                }
            });
            for (k, qs, qe, y, y_end) in acc {
                let win_len = lay.xe[k] - lay.xs[k];
                let (b, e) = (y - lay.xs[k], y_end - lay.xs[k]);
                if qe <= qs || (e - b) * 2 < win_len {
                    continue;
                }
                let id = *read_id.get_or_insert_with(|| {
                    if let Some(&id) = names.get(f[0]) {
                        id
                    } else {
                        let id = names.len() as u32;
                        names.insert(f[0].to_owned(), id);
                        id
                    }
                });
                wins[k].push(Seg {
                    read: id,
                    rev,
                    qs: qs as u32,
                    qe: qe as u32,
                    qlen: qlen as u32,
                    b: b as u32,
                    e: e as u32,
                    ident,
                });
            }
        }
    }
    let _ = std::fs::remove_file(&paf_path);
    let _ = std::fs::remove_file(&ref_path);

    // ---- 4. keep the best `max_depth` segments per window ----------------
    let mut needed = vec![false; names.len()];
    for (ci, per_contig) in windows.iter_mut().enumerate() {
        let Some(lay) = layouts[ci].as_ref() else {
            continue;
        };
        for (k, segs) in per_contig.iter_mut().enumerate() {
            let wl = (lay.xe[k] - lay.xs[k]) as u32;
            if segs.len() > p.max_depth {
                segs.sort_by(|a, b| {
                    let fa = a.b <= SLACK && a.e + SLACK >= wl;
                    let fb = b.b <= SLACK && b.e + SLACK >= wl;
                    fb.cmp(&fa)
                        .then(b.ident.partial_cmp(&a.ident).unwrap())
                        .then(a.read.cmp(&b.read))
                });
                segs.truncate(p.max_depth);
            }
            for s in segs.iter() {
                needed[s.read as usize] = true;
            }
        }
    }

    // ---- 5. load the sequences of the selected reads ---------------------
    let mut reads: Vec<Vec<u8>> = vec![Vec::new(); names.len()];
    {
        let mut reader = crate::utils::open_fastq_reader(reads_path)?;
        while let Some((header, seq, _, _)) = crate::utils::read_fastq_record(&mut reader)? {
            let name = header[1..].split_whitespace().next().unwrap_or("");
            if let Some(&id) = names.get(name) {
                if needed[id as usize] && reads[id as usize].is_empty() {
                    reads[id as usize] = seq.bytes().map(encode).collect();
                }
            }
        }
    }

    // ---- 6. per-window consensus (multi-threaded) ------------------------
    let mut tasks: Vec<(usize, usize)> = Vec::new();
    for (ci, per_contig) in windows.iter().enumerate() {
        for k in 0..per_contig.len() {
            tasks.push((ci, k));
        }
    }
    let next = AtomicUsize::new(0);
    let contigs_ro: &[Contig] = contigs;
    let mut results: Vec<(usize, Option<Vec<u8>>)> = Vec::with_capacity(tasks.len());
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..threads {
            handles.push(scope.spawn(|| {
                let mut out = Vec::new();
                let mut poa = Poa::new();
                let mut scratch = DpScratch::default();
                loop {
                    let t = next.fetch_add(1, Ordering::Relaxed);
                    if t >= tasks.len() {
                        break;
                    }
                    let (ci, k) = tasks[t];
                    let lay = layouts[ci].as_ref().unwrap();
                    let draft: Vec<u8> = contigs_ro[ci].seq[lay.xs[k]..lay.xe[k]]
                        .iter()
                        .map(|&b| encode(b))
                        .collect();
                    let segs = &windows[ci][k];
                    let res = if segs.len() < p.min_depth {
                        None
                    } else {
                        let layers: Vec<Layer> = segs
                            .iter()
                            .map(|s| make_layer(s, &reads[s.read as usize]))
                            .collect();
                        // keep only the bases whose backbone anchor lies in the core
                        let (c0, c1) = ((lay.bounds[k] - lay.xs[k]) as u32, (lay.bounds[k + 1] - lay.xs[k]) as u32);
                        poa_window(&mut poa, &mut scratch, &draft, &layers).map(|cons| {
                            cons.into_iter()
                                .filter(|&(_, a)| a >= c0 && a < c1)
                                .map(|x| x.0)
                                .collect::<Vec<u8>>()
                        })
                    };
                    out.push((t, res));
                }
                out
            }));
        }
        for h in handles {
            results.extend(h.join().expect("consensus worker panicked"));
        }
    });
    results.sort_by_key(|r| r.0);

    // ---- 7. stitch windows back together ---------------------------------
    let mut stats = ConsensusStats::default();
    let mut new_seqs: Vec<Vec<u8>> = contigs.iter().map(|_| Vec::new()).collect();
    for (t, res) in results {
        let (ci, k) = tasks[t];
        let lay = layouts[ci].as_ref().unwrap();
        let draft = &contigs[ci].seq[lay.bounds[k]..lay.bounds[k + 1]];
        stats.windows_total += 1;
        match res {
            Some(cons) => {
                let cons: Vec<u8> = cons.into_iter().map(decode).collect();
                if cons != draft {
                    stats.windows_changed += 1;
                }
                new_seqs[ci].extend_from_slice(&cons);
            }
            None => {
                stats.windows_unsupported += 1;
                new_seqs[ci].extend_from_slice(draft);
            }
        }
    }
    for (c, s) in contigs.iter_mut().zip(new_seqs.into_iter()) {
        if !s.is_empty() {
            c.seq = s;
        }
    }
    Ok(stats)
}

// ---------------------------------------------------------------------------
// CIGAR splitting
// ---------------------------------------------------------------------------

/// Walk a CIGAR and report, for every stretch of the target between two
/// consecutive `bounds` (plus the partial first/last stretch), the matching
/// target interval `[ts,te)` and query interval `[qs,qe)`.
/// Insertions that sit exactly on a border go to the window on the right.
fn split_cigar<F: FnMut(usize, usize, usize, usize)>(
    cigar: &str,
    tstart: usize,
    qstart: usize,
    bounds: &[usize],
    mut emit: F,
) {
    let mut t = tstart;
    let mut q = qstart;
    let mut seg_t = t;
    let mut seg_q = q;
    let mut bi = bounds.partition_point(|&b| b <= t);
    let mut num = 0usize;
    for ch in cigar.bytes() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch - b'0') as usize;
            continue;
        }
        let n = num;
        num = 0;
        match ch {
            b'M' | b'=' | b'X' | b'D' | b'N' => {
                let consumes_q = matches!(ch, b'M' | b'=' | b'X');
                let mut remaining = n;
                loop {
                    let room = bounds.get(bi).map(|&b| b - t);
                    match room {
                        Some(room) if remaining >= room => {
                            t += room;
                            if consumes_q {
                                q += room;
                            }
                            remaining -= room;
                            emit(seg_t, t, seg_q, q);
                            seg_t = t;
                            seg_q = q;
                            bi += 1;
                            if remaining == 0 {
                                break;
                            }
                        }
                        _ => {
                            t += remaining;
                            if consumes_q {
                                q += remaining;
                            }
                            break;
                        }
                    }
                }
            }
            b'I' | b'S' => q += n,
            _ => {}
        }
    }
    if t > seg_t {
        emit(seg_t, t, seg_q, q);
    }
}

// ---------------------------------------------------------------------------
// Sequence helpers (codes: A=0 C=1 G=2 T=3 N=4)
// ---------------------------------------------------------------------------

fn encode(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

fn decode(c: u8) -> u8 {
    b"ACGTN"[c.min(4) as usize]
}

struct Layer {
    seq: Vec<u8>,
    b: u32,
    e: u32,
}

fn make_layer(s: &Seg, read: &[u8]) -> Layer {
    let seq: Vec<u8> = if s.rev {
        let fs = (s.qlen - s.qe) as usize;
        let fe = (s.qlen - s.qs) as usize;
        read[fs..fe]
            .iter()
            .rev()
            .map(|&c| if c < 4 { 3 - c } else { 4 })
            .collect()
    } else {
        read[s.qs as usize..s.qe as usize].to_vec()
    };
    Layer {
        seq,
        b: s.b,
        e: s.e,
    }
}

// ---------------------------------------------------------------------------
// Partial-order alignment
// ---------------------------------------------------------------------------

const NONE: u32 = u32::MAX;

struct Poa {
    base: Vec<u8>,
    /// Backbone position the node belongs to (insertions inherit the anchor of
    /// the node before them).
    anchor: Vec<u32>,
    inp: Vec<Vec<(u32, u32)>>,
    outp: Vec<Vec<(u32, u32)>>,
    /// Nodes in the same alignment column (mismatch variants).
    sib: Vec<Vec<u32>>,
    /// Number of read segments (backbone excluded) that pass through the node.
    cov: Vec<u32>,
}

#[derive(Default)]
struct DpScratch {
    h: Vec<i32>,
    tb: Vec<u32>,
    rows: Vec<u32>,
    row_of: Vec<u32>,
    js: Vec<u32>,
    indeg: Vec<u32>,
    /// Substitution scores per base code (rows 0..5) and read position.
    prof: Vec<i32>,
    /// Running maxima of the insertion pass.
    pv: Vec<i32>,
}

impl Poa {
    fn new() -> Poa {
        Poa {
            base: Vec::new(),
            anchor: Vec::new(),
            inp: Vec::new(),
            outp: Vec::new(),
            sib: Vec::new(),
            cov: Vec::new(),
        }
    }

    fn reset(&mut self, backbone: &[u8]) {
        self.base.clear();
        self.anchor.clear();
        self.inp.clear();
        self.outp.clear();
        self.sib.clear();
        self.cov.clear();
        for (i, &b) in backbone.iter().enumerate() {
            let v = self.new_node(b, i as u32);
            if i > 0 {
                self.add_edge(v - 1, v);
            }
        }
    }

    fn new_node(&mut self, b: u8, anchor: u32) -> u32 {
        self.base.push(b);
        self.anchor.push(anchor);
        self.inp.push(Vec::new());
        self.outp.push(Vec::new());
        self.sib.push(Vec::new());
        self.cov.push(0);
        (self.base.len() - 1) as u32
    }

    fn add_edge(&mut self, a: u32, b: u32) {
        if let Some(e) = self.outp[a as usize].iter_mut().find(|e| e.0 == b) {
            e.1 += 1;
            if let Some(e2) = self.inp[b as usize].iter_mut().find(|e| e.0 == a) {
                e2.1 += 1;
            }
        } else {
            self.outp[a as usize].push((b, 1));
            self.inp[b as usize].push((a, 1));
        }
    }

    /// Kahn topological order; `None` if the graph contains a cycle.
    fn topo(&self, indeg: &mut Vec<u32>) -> Option<Vec<u32>> {
        let n = self.base.len();
        indeg.clear();
        indeg.extend(self.inp.iter().map(|v| v.len() as u32));
        let mut order: Vec<u32> = (0..n as u32).filter(|&v| indeg[v as usize] == 0).collect();
        let mut head = 0;
        while head < order.len() {
            let v = order[head] as usize;
            head += 1;
            for &(w, _) in &self.outp[v] {
                indeg[w as usize] -= 1;
                if indeg[w as usize] == 0 {
                    order.push(w);
                }
            }
        }
        if order.len() == n { Some(order) } else { None }
    }

    /// Align `seq` globally to the part of the graph whose anchors lie in
    /// `[lo, hi]` (global in the read, free start/end in the graph subject to
    /// the anchor limits).  The DP is banded: a node with backbone anchor `a`
    /// is only paired with read positions within `band` (`BAND` plus the length
    /// difference of read segment and backbone span) of the linear
    /// interpolation `(a - b0) * m / span`.  Returns (node or -1, read
    /// position or -1) pairs.
    fn align(
        &self,
        seq: &[u8],
        order: &[u32],
        (lo, hi, start_max_anchor, end_min_anchor): (u32, u32, u32, u32),
        (b0, span): (u32, u32),
        sc: &mut DpScratch,
    ) -> Option<Vec<(i32, i32)>> {
        const NEG: i32 = i32::MIN / 2;
        let m = seq.len();
        let band = (BAND + m.abs_diff(span as usize)).min(MAX_BAND);
        let ww = 2 * band + 2;
        sc.rows.clear();
        sc.row_of.clear();
        sc.row_of.resize(self.base.len(), 0);
        for &v in order {
            let a = self.anchor[v as usize];
            if a >= lo && a <= hi {
                sc.rows.push(v);
                sc.row_of[v as usize] = sc.rows.len() as u32; // 1-based
            }
        }
        let nr = sc.rows.len();
        if nr == 0 {
            return None;
        }
        // Row r holds the columns [js[r], js[r] + ww) of the DP; row 0 is the
        // virtual start (h = j * GAP).
        // Every active cell is written before it is read; get() returns NEG outside the band.
        if sc.h.len() < (nr + 1) * ww {
            sc.h.resize((nr + 1) * ww, NEG);
            sc.tb.resize((nr + 1) * ww, 0);
        }
        sc.js.clear();
        sc.js.push(0);
        for r in 1..=nr {
            let a = self.anchor[sc.rows[r - 1] as usize] as i64 - b0 as i64;
            let c = (a * m as i64 / span.max(1) as i64).clamp(0, m as i64);
            sc.js.push((c - band as i64).max(0) as u32);
        }
        let js = std::mem::take(&mut sc.js);
        let get = |h: &[i32], p: usize, j: usize| -> i32 {
            if p == 0 {
                return j as i32 * GAP;
            }
            let o = j.wrapping_sub(js[p] as usize);
            if o < ww { h[p * ww + o] } else { NEG }
        };
        sc.prof.clear();
        for b in 0..5u8 {
            sc.prof.extend(
                seq.iter()
                    .map(|&x| if x == b && b < 4 { MATCH } else { MISMATCH }),
            );
        }
        let mut preds: Vec<u32> = Vec::with_capacity(8);
        for r in 1..=nr {
            let v = sc.rows[r - 1] as usize;
            preds.clear();
            for &(p, _) in &self.inp[v] {
                let pr = sc.row_of[p as usize];
                if pr != 0 {
                    preds.push(pr);
                }
            }
            if preds.is_empty() || self.anchor[v] <= start_max_anchor {
                preds.push(0);
            }
            let bv = self.base[v] as usize;
            let j0 = js[r] as usize;
            let jend = (j0 + ww - 1).min(m);
            let len = jend - j0 + 1;
            let (done, rest) = sc.h.split_at_mut(r * ww);
            let (_, tb_rest) = sc.tb.split_at_mut(r * ww);
            let cur = &mut rest[..len];
            let cur_tb = &mut tb_rest[..len];
            cur.fill(NEG);
            cur_tb.fill(0);

            // Apply candidates in the original order; strict improvements preserve tie-breaking.
            let jd0 = j0.max(1);
            for &p in &preds {
                let p = p as usize;
                let code = (p as u32) << 2;
                if jd0 > jend {
                    continue;
                }
                let sub = &sc.prof[bv * m + jd0 - 1..bv * m + jend];
                if p == 0 {
                    let (hc, tc) = (&mut cur[jd0 - j0..], &mut cur_tb[jd0 - j0..]);
                    for (((h, t), &sb), c) in
                        hc.iter_mut().zip(tc.iter_mut()).zip(sub).zip(jd0 - 1..)
                    {
                        let d = c as i32 * GAP + sb;
                        *t = if d > *h { code } else { *t };
                        *h = (*h).max(d);
                    }
                } else {
                    let jp = js[p] as usize;
                    let lo = (jp + 1).clamp(jd0, jend + 1);
                    let hi = (jp + ww + 1).clamp(lo, jend + 1);
                    for j in (jd0..lo).chain(hi..=jend) {
                        let d = NEG + sub[j - jd0];
                        if d > cur[j - j0] {
                            cur[j - j0] = d;
                            cur_tb[j - j0] = code;
                        }
                    }
                    if lo < hi {
                        let src = &done[p * ww + lo - 1 - jp..p * ww + hi - 1 - jp];
                        let sub_in = &sub[lo - jd0..hi - jd0];
                        let (hin, tin) =
                            (&mut cur[lo - j0..hi - j0], &mut cur_tb[lo - j0..hi - j0]);
                        for (((h, t), &sv), &sb) in
                            hin.iter_mut().zip(tin.iter_mut()).zip(src).zip(sub_in)
                        {
                            let d = sv + sb;
                            *t = if d > *h { code } else { *t };
                            *h = (*h).max(d);
                        }
                    }
                }
            }

            for &p in &preds {
                let p = p as usize;
                let code = ((p as u32) << 2) | 1;
                if p == 0 {
                    for ((h, t), c) in cur.iter_mut().zip(cur_tb.iter_mut()).zip(j0..) {
                        let d = c as i32 * GAP + GAP;
                        *t = if d > *h { code } else { *t };
                        *h = (*h).max(d);
                    }
                } else {
                    let jp = js[p] as usize;
                    let lo = jp.clamp(j0, jend + 1);
                    let hi = (jp + ww).clamp(lo, jend + 1);
                    if lo < hi {
                        let src = &done[p * ww + lo - jp..p * ww + hi - jp];
                        let (hin, tin) =
                            (&mut cur[lo - j0..hi - j0], &mut cur_tb[lo - j0..hi - j0]);
                        for ((h, t), &sv) in hin.iter_mut().zip(tin.iter_mut()).zip(src) {
                            let d = sv + GAP;
                            *t = if d > *h { code } else { *t };
                            *h = (*h).max(d);
                        }
                    }
                }
            }

            if len > 1 {
                let ins_code = ((r as u32) << 2) | 2;
                sc.pv.clear();
                sc.pv.reserve(len);
                let mut run = cur[0];
                sc.pv.push(run);
                for (i, &h) in cur.iter().enumerate().skip(1) {
                    sc.pv.push(run);
                    run = run.max(h - i as i32 * GAP);
                }
                let ins = &sc.pv[1..];
                for (((h, t), &mx), i) in cur[1..]
                    .iter_mut()
                    .zip(cur_tb[1..].iter_mut())
                    .zip(ins)
                    .zip(1i32..)
                {
                    let d = mx + i * GAP;
                    *t = if d > *h { ins_code } else { *t };
                    *h = (*h).max(d);
                }
            }
        }
        // best end
        let mut best_r = 0usize;
        let mut best_s = NEG;
        for pass in 0..2 {
            for r in 1..=nr {
                let v = sc.rows[r - 1] as usize;
                if pass == 0 && self.anchor[v] < end_min_anchor {
                    continue;
                }
                let s = get(&sc.h, r, m);
                if s > best_s {
                    best_s = s;
                    best_r = r;
                }
            }
            if best_r != 0 {
                break;
            }
        }
        if best_r == 0 {
            sc.js = js;
            return None;
        }
        // traceback
        let mut out: Vec<(i32, i32)> = Vec::with_capacity(m + 16);
        let (mut r, mut j) = (best_r, m);
        while r != 0 || j != 0 {
            if r == 0 {
                out.push((-1, j as i32 - 1));
                j -= 1;
                continue;
            }
            let code = sc.tb[r * ww + j - js[r] as usize];
            let p = (code >> 2) as usize;
            match code & 3 {
                0 => {
                    out.push((sc.rows[r - 1] as i32, j as i32 - 1));
                    r = p;
                    j -= 1;
                }
                1 => {
                    out.push((sc.rows[r - 1] as i32, -1));
                    r = p;
                }
                _ => {
                    out.push((-1, j as i32 - 1));
                    j -= 1;
                }
            }
        }
        sc.js = js;
        out.reverse();
        Some(out)
    }

    /// Merge an alignment into the graph.
    fn add_alignment(&mut self, seq: &[u8], aln: &[(i32, i32)]) {
        let mut prev: Option<u32> = None;
        for &(node, rp) in aln {
            if rp < 0 {
                continue;
            }
            let b = seq[rp as usize];
            let cur = if node >= 0 {
                let v = node as u32;
                if self.base[v as usize] == b {
                    v
                } else if let Some(&s) = self.sib[v as usize]
                    .iter()
                    .find(|&&s| self.base[s as usize] == b)
                {
                    s
                } else {
                    let nn = self.new_node(b, self.anchor[v as usize]);
                    let mut group = self.sib[v as usize].clone();
                    group.push(v);
                    for &g in &group {
                        self.sib[g as usize].push(nn);
                    }
                    self.sib[nn as usize] = group;
                    nn
                }
            } else {
                let a = prev.map_or(0, |p| self.anchor[p as usize]);
                self.new_node(b, a)
            };
            self.cov[cur as usize] += 1;
            if let Some(p) = prev {
                self.add_edge(p, cur);
            }
            prev = Some(cur);
        }
    }

    /// Heaviest path through the graph.
    /// `span[a]` is the number of read segments spanning backbone position `a`.
    fn consensus(&self, order: &[u32], min_support: f64, span: &[u32]) -> Vec<(u8, u32)> {
        let n = self.base.len();
        let mut score = vec![0u64; n];
        let mut back = vec![NONE; n];
        let mut best_node = order[0];
        let mut best_score = 0u64;
        for &v in order {
            let v = v as usize;
            let mut bs = 0u64;
            let mut bp = NONE;
            for &(p, w) in &self.inp[v] {
                let c = score[p as usize] + w as u64;
                if c > bs {
                    bs = c;
                    bp = p;
                }
            }
            score[v] = bs;
            back[v] = bp;
            if bs > best_score {
                best_score = bs;
                best_node = v as u32;
            }
        }
        // Trace back the heaviest path, remembering how many reads support
        // each node.
        let mut nodes: Vec<(u8, u32, u32)> = Vec::new();
        let mut v = best_node;
        while v != NONE {
            nodes.push((self.base[v as usize], self.cov[v as usize], self.anchor[v as usize]));
            v = back[v as usize];
        }
        nodes.reverse();
        // Heaviest-path scoring favours long paths and nodes supported by
        // only ~1/3 of the reads, which keeps spurious bases (including
        // unsupported draft bases at the window ends).  Enforce a majority
        // rule: drop path nodes supported by less than `min_support` times
        // the local (median) read depth.
        let n = nodes.len();
        let half = 20usize;
        let mut out = Vec::with_capacity(n);
        let mut buf: Vec<u32> = Vec::with_capacity(2 * half + 1);
        for i in 0..n {
            buf.clear();
            buf.extend(
                nodes[i.saturating_sub(half)..(i + half + 1).min(n)]
                    .iter()
                    .map(|x| x.1),
            );
            let median = buf.len() / 2;
            let depth = *buf.select_nth_unstable(median).1;
            // The reference depth is the larger of the local median coverage of the path and
            // the number of read segments spanning the node's backbone position. The median
            // alone is 0 inside any run of 21+ unsupported nodes (the window holds 41 nodes),
            // and `0 >= 0.5 * 0` kept every base of such a run: poly-A/T stretches and short
            // repeats of the draft that no read supports survived the consensus.
            let spanning = span.get(nodes[i].2 as usize).copied().unwrap_or(0);
            let depth = depth.max(spanning);
            if (nodes[i].1 as f64) >= min_support * depth as f64 {
                out.push((nodes[i].0, nodes[i].2));
            }
        }
        out
    }
}

/// Build and solve the POA of one window.  `None` keeps the draft.
fn poa_window(poa: &mut Poa, sc: &mut DpScratch, draft: &[u8], layers: &[Layer]) -> Option<Vec<(u8, u32)>> {
    if draft.is_empty() || layers.is_empty() {
        return None;
    }
    poa.reset(draft);
    let mut indeg = std::mem::take(&mut sc.indeg);
    for layer in layers {
        if layer.seq.len() < 10 {
            continue;
        }
        let order = poa.topo(&mut indeg)?;
        let lo = layer.b.saturating_sub(SLACK);
        let hi = layer.e + SLACK;
        let start_max = layer.b + SLACK;
        let end_min = layer.e.saturating_sub(1).saturating_sub(SLACK);
        let aln = poa.align(&layer.seq, &order, (lo, hi, start_max, end_min), (layer.b, layer.e - layer.b), sc)?;
        poa.add_alignment(&layer.seq, &aln);
    }
    let order = poa.topo(&mut indeg)?;
    sc.indeg = indeg;
    // number of read segments spanning each backbone position
    let mut diff = vec![0i32; draft.len() + 1];
    for layer in layers.iter().filter(|l| l.seq.len() >= 10) {
        let (b, e) = (layer.b as usize, (layer.e as usize).min(draft.len()));
        if b < e {
            diff[b] += 1;
            diff[e] -= 1;
        }
    }
    let mut run = 0i32;
    let span: Vec<u32> = diff
        .iter()
        .map(|d| {
            run += d;
            run.max(0) as u32
        })
        .collect();
    let cons = poa.consensus(&order, MIN_SUPPORT, &span);
    // sanity guard: keep the draft if the consensus is implausible. The expected length is the
    // median length of the read segments that span the whole window (they differ from the draft
    // length when the draft has an indel); without such reads the draft length is used.
    let mut full: Vec<usize> = layers
        .iter()
        .filter(|l| l.b <= SLACK && l.e + SLACK >= draft.len() as u32)
        .map(|l| l.seq.len())
        .collect();
    let expected = if full.len() >= 3 {
        full.sort_unstable();
        full[full.len() / 2]
    } else {
        draft.len()
    };
    let tolerance = expected / 8 + 20;
    if cons.len() * 2 < draft.len().min(expected)
        || (cons.len() > draft.len() * 3 / 2 + 20 && cons.len().abs_diff(expected) > tolerance)
    {
        return None;
    }
    Some(cons)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn f(&mut self) -> f64 {
            (self.next() % 1_000_000) as f64 / 1_000_000.0
        }
    }

    fn mutate(rng: &mut Lcg, truth: &[u8], rate: f64) -> Vec<u8> {
        let mut out = Vec::new();
        for &b in truth {
            let r = rng.f();
            if r < rate / 3.0 {
                // deletion
            } else if r < 2.0 * rate / 3.0 {
                out.push(b);
                out.push((rng.next() % 4) as u8);
            } else if r < rate {
                out.push(((b as u64 + 1 + rng.next() % 3) % 4) as u8);
            } else {
                out.push(b);
            }
        }
        out
    }

    #[test]
    fn poa_recovers_truth_from_noisy_reads() {
        let mut rng = Lcg(7);
        let truth: Vec<u8> = (0..400).map(|_| (rng.next() % 4) as u8).collect();
        let draft = mutate(&mut rng, &truth, 0.08);
        let layers: Vec<Layer> = (0..25)
            .map(|_| {
                let s = mutate(&mut rng, &truth, 0.08);
                Layer {
                    b: 0,
                    e: draft.len() as u32,
                    seq: s,
                }
            })
            .collect();
        let mut poa = Poa::new();
        let mut sc = DpScratch::default();
        let cons: Vec<u8> = poa_window(&mut poa, &mut sc, &draft, &layers).unwrap().into_iter().map(|x| x.0).collect();
        assert_eq!(cons, truth, "consensus should equal the truth at 25x");
    }

    /// A draft that lacks (or has extra) hundreds of bases relative to the reads is repaired.
    #[test]
    fn poa_repairs_large_indels_in_draft() {
        let mut rng = Lcg(11);
        let truth: Vec<u8> = (0..800).map(|_| (rng.next() % 4) as u8).collect();
        for (missing, extra) in [(450usize, 0usize), (0, 300), (200, 0)] {
            let mut draft: Vec<u8> = truth.clone();
            if missing > 0 {
                draft.drain(300..300 + missing);
            }
            for _ in 0..extra {
                draft.insert(300, (rng.next() % 4) as u8);
            }
            let layers: Vec<Layer> = (0..25)
                .map(|_| Layer { b: 0, e: draft.len() as u32, seq: mutate(&mut rng, &truth, 0.06) })
                .collect();
            let mut poa = Poa::new();
            let mut sc = DpScratch::default();
            let cons: Vec<u8> = poa_window(&mut poa, &mut sc, &draft, &layers).expect("window re-called").into_iter().map(|x| x.0).collect();
            let diff = cons.iter().zip(&truth).filter(|(a, b)| a != b).count();
            assert!(cons.len().abs_diff(truth.len()) <= 2 && diff <= 3, "missing {missing}, extra {extra}: len {} vs {}, {diff} differ", cons.len(), truth.len());
        }
    }

    fn edit_distance(a: &[u8], b: &[u8]) -> usize {
        let mut prev: Vec<usize> = (0..=b.len()).collect();
        for i in 1..=a.len() {
            let mut cur = vec![i; b.len() + 1];
            for j in 1..=b.len() {
                cur[j] = (prev[j - 1] + (a[i - 1] != b[j - 1]) as usize)
                    .min(prev[j] + 1)
                    .min(cur[j - 1] + 1);
            }
            prev = cur;
        }
        prev[b.len()]
    }

    /// Runs of low-complexity bases (poly-A, short repeats) that sit in the draft but in none
    /// of the reads must be removed, whatever their length. Before, any run of 21 or more
    /// draft-only bases survived: the local-median coverage filter is 0 inside such a run, so
    /// `cov >= 0.5 * 0` kept every node.
    #[test]
    fn poa_removes_unsupported_runs_from_draft() {
        for (kind, motif) in [("poly-A", &[0u8][..]), ("poly-T", &[3u8][..]), ("ATTCC", &[0u8, 3, 3, 1, 1][..])] {
            for depth in [8usize, 25] {
                for junk in [10usize, 21, 22, 25, 30, 60, 100, 300, 450] {
                    let mut rng = Lcg(5 + junk as u64);
                    let truth: Vec<u8> = (0..500).map(|_| (rng.next() % 4) as u8).collect();
                    let mut draft = truth.clone();
                    for i in 0..junk {
                        draft.insert(250 + i, motif[i % motif.len()]);
                    }
                    let layers: Vec<Layer> = (0..depth)
                        .map(|_| Layer { b: 0, e: draft.len() as u32, seq: mutate(&mut rng, &truth, 0.06) })
                        .collect();
                    let mut poa = Poa::new();
                    let mut sc = DpScratch::default();
                    let cons: Vec<u8> = poa_window(&mut poa, &mut sc, &draft, &layers)
                        .unwrap_or_else(|| panic!("{kind} x{junk}, depth {depth}: window rejected"))
                        .into_iter()
                        .map(|x| x.0)
                        .collect();
                    let d = edit_distance(&cons, &truth);
                    assert!(
                        d <= 4,
                        "{kind} x{junk}, depth {depth}: consensus {} bp vs truth {} bp (edit distance {d})",
                        cons.len(),
                        truth.len()
                    );
                }
            }
        }
    }

    #[test]
    fn split_cigar_cuts_at_bounds() {
        let mut segs = Vec::new();
        // 10M2I10M3D10M starting at t=95 / q=0, borders at 100 and 120
        split_cigar("10M2I10M3D10M", 95, 0, &[100, 120], |a, b, c, d| {
            segs.push((a, b, c, d))
        });
        assert_eq!(
            segs,
            vec![(95, 100, 0, 5), (100, 120, 5, 24), (120, 128, 24, 32)]
        );
    }

    #[test]
    fn circular_layout_maps_flanks_back() {
        let lay = Layout::new(10_000, true, 500, 0);
        assert_eq!(lay.flank, CIRC_FLANK);
        // left flank position 0 == contig base 8000
        assert_eq!(lay.to_contig(0, 100), Some((8000, 8100)));
        // right flank == contig base 0
        assert_eq!(lay.to_contig(lay.flank + 10_000, lay.flank + 10_100), Some((0, 100)));
    }
}