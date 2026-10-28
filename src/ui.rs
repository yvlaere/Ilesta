//! Console output: numbered stages, nested timed steps and short key figures.
//!
//! `NO_COLOR` is unset and `TERM` is not `dumb`.
//!
//! ```text
//! ━━ [2/6] Overlap filtering ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//! ▸ Filter read overlaps
//!   │ alignments: 1,204,331 kept
//!   └ 2.8 s
//! ```

use std::fmt::Display;
use std::io::IsTerminal;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

const RULE_WIDTH: usize = 64;

static DEPTH: AtomicUsize = AtomicUsize::new(0);
static STAGE_NUMBER: AtomicUsize = AtomicUsize::new(0);
static STAGE_TOTAL: AtomicUsize = AtomicUsize::new(0);
static CURRENT_STAGE: Mutex<Option<(String, Instant)>> = Mutex::new(None);
static FINISHED_STAGES: Mutex<Vec<(String, f64)>> = Mutex::new(Vec::new());

fn use_colour(is_terminal: bool) -> bool {
    is_terminal
        && std::env::var_os("NO_COLOR").is_none_or(|v| v.is_empty())
        && std::env::var("TERM").is_ok_and(|t| t != "dumb")
}

fn paint(code: &str, text: &str, stderr: bool) -> String {
    static COLOUR: OnceLock<(bool, bool)> = OnceLock::new();
    let (out, err) = *COLOUR.get_or_init(|| (use_colour(std::io::stdout().is_terminal()), use_colour(std::io::stderr().is_terminal())));
    if if stderr { err } else { out } { format!("\x1b[{code}m{text}\x1b[0m") } else { text.to_string() }
}

pub fn bold(text: &str) -> String {
    paint("1", text, false)
}

pub fn dim(text: &str) -> String {
    paint("2", text, false)
}

pub fn green(text: &str) -> String {
    paint("32", text, false)
}

fn indent(extra: usize) -> String {
    " ".repeat(2 * DEPTH.load(Relaxed) + extra)
}

/// Set how many stages the run has (shown as `[n/total]`) and start numbering from 1.
pub fn plan(total: usize) {
    STAGE_TOTAL.store(total, Relaxed);
    STAGE_NUMBER.store(0, Relaxed);
    DEPTH.store(0, Relaxed);
    FINISHED_STAGES.lock().unwrap().clear();
}

/// A free-standing block (the run header, the summary).
pub fn banner(title: &str) {
    finish_stage();
    println!("\n{}", rule(title));
}

fn rule(title: &str) -> String {
    let used = title.chars().count() + 4;
    paint("1;36", &format!("━━ {title} {}", "━".repeat(RULE_WIDTH.saturating_sub(used))), false)
}

/// Close the running stage (recording its time) and open the next one.
pub fn stage(title: &str) {
    finish_stage();
    let number = STAGE_NUMBER.fetch_add(1, Relaxed) + 1;
    let total = STAGE_TOTAL.load(Relaxed);
    println!("\n{}", rule(&format!("[{number}/{total}] {title}")));
    *CURRENT_STAGE.lock().unwrap() = Some((title.to_string(), Instant::now()));
}

pub fn finish_stage() {
    if let Some((title, started)) = CURRENT_STAGE.lock().unwrap().take() {
        FINISHED_STAGES.lock().unwrap().push((title, started.elapsed().as_secs_f64()));
    }
}

/// Time per stage, in pipeline order, with the total.
pub fn stage_times() -> Vec<(String, f64)> {
    finish_stage();
    FINISHED_STAGES.lock().unwrap().clone()
}

/// Run `operation` as a named step. Lines printed while it runs are indented under it, and it
/// ends with its elapsed time (or the error, on stderr).
pub fn step<T, E: Display>(label: &str, operation: impl FnOnce() -> Result<T, E>) -> Result<T, E> {
    println!("{}{} {}", indent(0), paint("36", "▸", false), bold(label));
    let started = Instant::now();
    DEPTH.fetch_add(1, Relaxed);
    let result = operation();
    DEPTH.fetch_sub(1, Relaxed);
    let elapsed = duration(started.elapsed().as_secs_f64());
    match &result {
        Ok(_) => println!("{}{} {}", indent(2), dim("└"), green(&elapsed)),
        Err(error) => eprintln!("{}{} {}", indent(2), paint("1;31", "✗", true), paint("31", &format!("{label} failed after {elapsed}: {error}"), true)),
    }
    result
}

/// A key figure or note belonging to the current step (or stage).
pub fn detail(text: impl Display) {
    if DEPTH.load(Relaxed) == 0 {
        println!("  {text}");
    } else {
        println!("{}{} {text}", indent(0), dim("│"));
    }
}

pub fn warn(text: impl Display) {
    println!("{}{} {text}", indent(0), paint("33", "!", false));
}

/// 1234567 -> "1,234,567"
pub fn count(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(' ');
        }
        out.push(c);
    }
    out
}

/// 5_962_000 -> "5.96 Mb", 18_200 -> "18.2 kb", 950 -> "950 bp"
pub fn bp(n: u64) -> String {
    match n {
        0..1_000 => format!("{n} bp"),
        1_000..1_000_000 => format!("{:.1} kb", n as f64 / 1e3),
        1_000_000..1_000_000_000 => format!("{:.2} Mb", n as f64 / 1e6),
        _ => format!("{:.2} Gb", n as f64 / 1e9),
    }
}

/// 3.2 -> "3.2 s", 410 -> "6m 50s", 3900 -> "1h 05m"
pub fn duration(seconds: f64) -> String {
    let whole = seconds.round() as u64;
    match seconds {
        s if s < 60.0 => format!("{s:.1} s"),
        s if s < 3600.0 => format!("{}m {:02}s", whole / 60, whole % 60),
        _ => format!("{}h {:02}m", whole / 3600, whole % 3600 / 60),
    }
}

/// N50 of a set of lengths.
pub fn n50(lengths: impl IntoIterator<Item = u64>) -> u64 {
    let mut lengths: Vec<u64> = lengths.into_iter().collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let half = lengths.iter().sum::<u64>().div_ceil(2);
    let mut sum = 0;
    for length in lengths {
        sum += length;
        if sum >= half {
            return length;
        }
    }
    0
}