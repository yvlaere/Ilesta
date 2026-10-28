use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

const SMOKE_READ_COUNT: usize = 10_000;

#[test]
#[ignore = "requires the local evaluation dataset and minimap2 on PATH"]
fn assembles_a_sample_from_the_evaluation_dataset() {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let dataset = manifest_dir.join("evaluation/test_data/SRR28262566.fastq");
    assert!(
        dataset.is_file(),
        "missing test dataset: {}",
        dataset.display()
    );

    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock is before Unix epoch")
        .as_nanos();
    let temp_dir = std::env::temp_dir().join(format!(
        "ilesta-assembler-smoke-{}-{timestamp}",
        std::process::id()
    ));
    let output_dir = temp_dir.join("output");
    fs::create_dir_all(&output_dir).expect("create smoke-test output directory");
    let sample_path = temp_dir.join("sample.fastq");
    let sampled_records = copy_first_records(&dataset, &sample_path, SMOKE_READ_COUNT)
        .expect("extract complete FASTQ records from dataset");
    assert!(sampled_records > 0, "dataset contains no FASTQ records");

    let result = Command::new(env!("CARGO_BIN_EXE_Ilesta"))
        .arg("assemble")
        .arg("--reads-fq")
        .arg(&sample_path)
        .arg("--output-dir")
        .arg(&output_dir)
        .arg("--threads")
        .arg("1")
        .arg("--genome-size")
        .arg("10000000")
        .arg("--min-read-length")
        .arg("100")
        .arg("--min-base-quality")
        .arg("0")
        .arg("--min-overlap-length")
        .arg("200")
        .arg("--min-overlap-count")
        .arg("1")
        .output()
        .expect("start Ilesta assembler");

    assert!(
        result.status.success(),
        "assembler exited with {}\nstdout:\n{}\nstderr:\n{}",
        result.status,
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.is_empty(),
        "minimap2 diagnostics should not be forwarded to stderr: {stderr}"
    );
    let stdout = String::from_utf8_lossy(&result.stdout);
    assert!(
        stdout.contains("command   Ilesta assemble"),
        "the command line was not reported"
    );
    assert!(
        stdout.lines().any(|line| line.contains("[1/")) && stdout.contains("Done"),
        "pipeline stage headings were not reported"
    );
    assert!(
        !stdout.contains("[M::"),
        "minimap2 diagnostics should not be forwarded to stdout"
    );
    // by default only the assembly itself is left behind
    let mut files: Vec<_> = fs::read_dir(&output_dir)
        .expect("list output directory")
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    files.sort();
    assert_eq!(files, ["unitigs.fa", "unitigs.gfa"]);
    let unitigs =
        fs::read_to_string(output_dir.join("unitigs.fa")).expect("read primary assembly FASTA");
    assert!(
        unitigs.lines().any(|line| line.starts_with(">unitig_")),
        "assembler completed without producing a unitig"
    );
    let fasta_names: std::collections::HashSet<_> = unitigs
        .lines()
        .filter(|line| line.starts_with('>'))
        .filter_map(|line| line[1..].split_whitespace().next())
        .collect();
    let fasta_topologies: std::collections::HashMap<_, _> = unitigs
        .lines()
        .filter(|line| line.starts_with('>'))
        .filter_map(|line| {
            let mut fields = line[1..].split_whitespace();
            let name = fields.next()?;
            let topology = fields.find_map(|field| field.strip_prefix("topology="))?;
            Some((name, topology))
        })
        .collect();
    let gfa = fs::read_to_string(output_dir.join("unitigs.gfa")).expect("read combined GFA");
    let gfa_names: Vec<_> = gfa
        .lines()
        .filter(|line| line.starts_with("S\t"))
        .filter_map(|line| line.split('\t').nth(1))
        .collect();
    let unique_gfa_names: std::collections::HashSet<_> = gfa_names.iter().copied().collect();
    assert_eq!(
        gfa_names.len(),
        unique_gfa_names.len(),
        "GFA segment names repeat"
    );
    assert_eq!(fasta_names.len(), gfa_names.len());
    for gfa_name in gfa_names {
        let (fasta_name, gfa_topology) = gfa_name
            .rsplit_once('_')
            .expect("GFA segment name includes topology");
        assert!(fasta_names.contains(fasta_name));
        assert_eq!(
            fasta_topologies.get(fasta_name).copied(),
            Some(gfa_topology)
        );
    }
    let rescued_unitig_count = stdout
        .lines()
        .find_map(|line| {
            line.split("unitigs containing rescued reads: ")
                .nth(1)
                .and_then(|rest| rest.trim().parse::<usize>().ok())
        })
        .expect("rescue unitig count was not reported");
    assert!(
        rescued_unitig_count > 0,
        "no rescued unitigs reached the combined assembly"
    );

    fs::remove_dir_all(temp_dir).expect("remove smoke-test temporary directory");
}

fn copy_first_records(input: &Path, output: &Path, max_records: usize) -> std::io::Result<usize> {
    let mut reader = BufReader::new(File::open(input)?);
    let mut writer = BufWriter::new(File::create(output)?);
    let mut copied = 0;

    while copied < max_records {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if !line.starts_with('@') {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "expected FASTQ header line",
            ));
        }
        writer.write_all(line.as_bytes())?;

        let mut sequence_bases = 0;
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "FASTQ record ended before plus line",
                ));
            }
            let content = line.trim_end_matches(['\r', '\n']);
            writer.write_all(line.as_bytes())?;
            if content.starts_with('+') {
                break;
            }
            sequence_bases += content.len();
        }

        let mut quality_bases = 0;
        while quality_bases < sequence_bases {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "FASTQ record ended before quality was complete",
                ));
            }
            quality_bases += line.trim_end_matches(['\r', '\n']).len();
            writer.write_all(line.as_bytes())?;
        }
        if quality_bases != sequence_bases {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "FASTQ sequence and quality lengths differ",
            ));
        }
        copied += 1;
    }

    writer.flush()?;
    Ok(copied)
}
