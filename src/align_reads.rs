// module for filtering and anigning reads before assembly

use std::process::{Command};



pub fn run_minimap2(query: &str, threads: usize, out_path: &str) -> std::io::Result<()> {
    let mut child = Command::new("minimap2")
        .args(&["-x", "ava-ont", "-t", &threads.to_string(), query, query, "-o", out_path])
        .spawn()?;

    let status = child.wait()?;
    assert!(status.success());

    Ok(())
}