//! Consensus (self-polishing) step
//!
//! Method:
//!   1. Map the assembly reads back to the unitigs with `minimap2 -c`.
//!   2. Cut every unitig into windows (default ~500 bp) and use each CIGAR to
//!      split every read alignment at the window borders.
//!   3. Per window, build a partial-order alignment (POA) graph: the draft
//!      window is the backbone, each read segment is aligned to the graph
//!      and merged in.
//!   4. The heaviest path through the POA graph is the new window sequence.
//!   5. The windows are concatenated; the whole procedure can be repeated
//!
//! Circular unitigs (`topology == 'c'`) are mapped with a flank of the
//! opposite end on each side so that reads crossing the origin contribute to
//! the windows on both sides of it, and the windows are tiled around the ring.
//!

use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::Command;
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
/// Minimum identity of a whole read alignment to be used.
const MIN_IDENTITY: f64 = 0.6;
/// A path node is kept only if at least this fraction of the local median
/// read depth supports it (majority rule).
const MIN_SUPPORT: f64 = 0.5;
/// Flank (bp) added to each side of circular unitigs while mapping.
const CIRC_FLANK: usize = 2000;
/// Unitigs shorter than this are left untouched.
const MIN_CONTIG_LEN: usize = 100;

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Polish all unitigs of `graph` in place and rewrite `fasta_path`.
///
/// `reads_path` is the FASTQ(.gz) that was used for assembly (it is only read,
/// never modified).  Intermediate files are written to `work_dir` and removed.
pub fn polish_graph(
    graph: &mut CompressedGraph,
    fasta_path: &Path,
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
        println!("Consensus round {}/{}...", round, params.rounds);
        let round_stats = polish_round(&mut contigs, reads_path, work_dir, params, round)?;
        stats.rounds_run += 1;
        stats.windows_total = round_stats.windows_total;
        stats.windows_changed += round_stats.windows_changed;
        stats.windows_unsupported = round_stats.windows_unsupported;
        println!(
            "  windows: {} total, {} changed, {} without enough support",
            round_stats.windows_total, round_stats.windows_changed, round_stats.windows_unsupported
        );
        if round_stats.windows_changed == 0 {
            break; // converged
        }
    }
    stats.bases_after = contigs.iter().map(|c| c.seq.len()).sum();

    for (unitig, contig) in graph.unitigs.iter_mut().zip(contigs.into_iter()) {
        unitig.fasta_seq = Some(String::from_utf8(contig.seq)?);
    }

    let mut fasta = std::fs::File::create(fasta_path)?;
    for unitig in &graph.unitigs {
        let sequence = unitig.fasta_seq.as_deref().unwrap_or("");
        writeln!(
            fasta,
            ">unitig_{} len={}bp topology={}",
            unitig.id,
            sequence.len(),
            unitig.topology
        )?;
        writeln!(fasta, "{sequence}")?;
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
    /// Window borders in contig coordinates: bounds[0]=0, bounds[K]=len.
    bounds: Vec<usize>,
    /// Border positions in coordinates of the (flanked) mapping reference.
    ext_bounds: Vec<usize>,
}

impl Layout {
    fn new(len: usize, circular: bool, window_len: usize) -> Layout {
        let flank = if circular { CIRC_FLANK.min(len / 2) } else { 0 };
        let k = ((len as f64 / window_len as f64).round() as usize).max(1);
        let bounds: Vec<usize> = (0..=k).map(|i| i * len / k).collect();
        let ext_len = len + 2 * flank;
        let mut ext_bounds = Vec::new();
        let shifts: &[isize] = if circular { &[-1, 0, 1] } else { &[0] };
        for &t in shifts {
            for &b in &bounds {
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
            ext_bounds,
        }
    }

    fn windows(&self) -> usize {
        self.bounds.len() - 1
    }

    /// Map an ext-coordinate interval (not crossing a border) to
    /// (window index, begin offset, end offset) in contig coordinates.
    fn canonical(&self, ts: usize, te: usize) -> Option<(usize, usize, usize)> {
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
        if y >= self.len {
            return None;
        }
        let k = self.bounds.partition_point(|&b| b <= y) - 1;
        if y_end > self.bounds[k + 1] {
            return None; // should not happen
        }
        Some((k, y - self.bounds[k], y_end - self.bounds[k]))
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
            let lay = Layout::new(c.seq.len(), c.circular, p.window_len.max(50));
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
    let status = Command::new("minimap2")
        .arg("-c")
        .arg("--secondary=no")
        .arg("-t")
        .arg(threads.to_string())
        .args(p.read_type.mapping_args())
        .arg(&ref_path)
        .arg(reads_path)
        .arg("-o")
        .arg(&paf_path)
        .stderr(std::process::Stdio::null())
        .status()?;
    if !status.success() {
        return Err("minimap2 failed during consensus step".into());
    }

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
            split_cigar(cigar, tstart, q0, &lay.ext_bounds, |ts, te, qs, qe| {
                let Some((k, b, e)) = lay.canonical(ts, te) else {
                    return;
                };
                let win_len = lay.bounds[k + 1] - lay.bounds[k];
                if qe <= qs || (e - b) * 2 < win_len {
                    return;
                }
                let id = *read_id.get_or_insert_with(|| {
                    let next = names.len() as u32;
                    *names.entry(f[0].to_string()).or_insert(next)
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
            });
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
            let wl = (lay.bounds[k + 1] - lay.bounds[k]) as u32;
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
                    let draft: Vec<u8> = contigs_ro[ci].seq[lay.bounds[k]..lay.bounds[k + 1]]
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
                        poa_window(&mut poa, &mut scratch, &draft, &layers)
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
    indeg: Vec<u32>,
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
    /// the anchor limits).  Returns (node or -1, read position or -1) pairs.
    fn align(
        &self,
        seq: &[u8],
        order: &[u32],
        lo: u32,
        hi: u32,
        start_max_anchor: u32,
        end_min_anchor: u32,
        sc: &mut DpScratch,
    ) -> Option<Vec<(i32, i32)>> {
        let m = seq.len();
        let w = m + 1;
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
        sc.h.clear();
        sc.h.resize((nr + 1) * w, i32::MIN / 2);
        sc.tb.clear();
        sc.tb.resize((nr + 1) * w, 0);
        // op: 0 diag, 1 del (node only), 2 ins (read only)
        for j in 0..w {
            sc.h[j] = j as i32 * GAP;
            sc.tb[j] = 2;
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
            let bv = self.base[v];
            // column 0
            {
                let mut best = i32::MIN / 2;
                let mut bp = 0u32;
                for &p in &preds {
                    let c = sc.h[p as usize * w] + GAP;
                    if c > best {
                        best = c;
                        bp = p;
                    }
                }
                sc.h[r * w] = best;
                sc.tb[r * w] = (bp << 2) | 1;
            }
            for j in 1..w {
                let s = if bv == seq[j - 1] && bv < 4 {
                    MATCH
                } else {
                    MISMATCH
                };
                let mut best = i32::MIN / 2;
                let mut code = 0u32;
                for &p in &preds {
                    let d = sc.h[p as usize * w + j - 1] + s;
                    if d > best {
                        best = d;
                        code = p << 2;
                    }
                }
                for &p in &preds {
                    let d = sc.h[p as usize * w + j] + GAP;
                    if d > best {
                        best = d;
                        code = (p << 2) | 1;
                    }
                }
                let ins = sc.h[r * w + j - 1] + GAP;
                if ins > best {
                    best = ins;
                    code = ((r as u32) << 2) | 2;
                }
                sc.h[r * w + j] = best;
                sc.tb[r * w + j] = code;
            }
        }
        // best end
        let mut best_r = 0usize;
        let mut best_s = i32::MIN;
        for r in 1..=nr {
            let v = sc.rows[r - 1] as usize;
            if self.anchor[v] < end_min_anchor {
                continue;
            }
            let s = sc.h[r * w + m];
            if s > best_s {
                best_s = s;
                best_r = r;
            }
        }
        if best_r == 0 {
            for r in 1..=nr {
                let s = sc.h[r * w + m];
                if s > best_s {
                    best_s = s;
                    best_r = r;
                }
            }
        }
        // traceback
        let mut out: Vec<(i32, i32)> = Vec::with_capacity(m + 16);
        let (mut r, mut j) = (best_r, m);
        while r != 0 || j != 0 {
            let code = sc.tb[r * w + j];
            let op = code & 3;
            let p = (code >> 2) as usize;
            match op {
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
    fn consensus(&self, order: &[u32], min_support: f64) -> Vec<u8> {
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
        let mut nodes: Vec<(u8, u32)> = Vec::new();
        let mut v = best_node;
        while v != NONE {
            nodes.push((self.base[v as usize], self.cov[v as usize]));
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
            buf.sort_unstable();
            let depth = buf[buf.len() / 2];
            if (nodes[i].1 as f64) >= min_support * depth as f64 {
                out.push(nodes[i].0);
            }
        }
        out
    }
}

/// Build and solve the POA of one window.  `None` keeps the draft.
fn poa_window(poa: &mut Poa, sc: &mut DpScratch, draft: &[u8], layers: &[Layer]) -> Option<Vec<u8>> {
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
        let aln = poa.align(&layer.seq, &order, lo, hi, start_max, end_min, sc)?;
        poa.add_alignment(&layer.seq, &aln);
    }
    let order = poa.topo(&mut indeg)?;
    sc.indeg = indeg;
    let cons = poa.consensus(&order, MIN_SUPPORT);
    // sanity guard: keep the draft if the consensus is implausible
    if cons.len() * 2 < draft.len() || cons.len() > draft.len() * 3 / 2 + 20 {
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
        let cons = poa_window(&mut poa, &mut sc, &draft, &layers).unwrap();
        assert_eq!(cons, truth, "consensus should equal the truth at 25x");
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
        let lay = Layout::new(10_000, true, 500);
        assert_eq!(lay.flank, CIRC_FLANK);
        // left flank position 0 == contig base 8000
        let (k, b, _) = lay.canonical(0, 100).unwrap();
        assert_eq!(lay.bounds[k] + b, 8000);
        // right flank == contig base 0
        let (k, b, _) = lay.canonical(lay.flank + 10_000, lay.flank + 10_100).unwrap();
        assert_eq!(lay.bounds[k] + b, 0);
    }
}