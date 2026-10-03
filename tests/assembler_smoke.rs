use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

const SMOKE_READ_COUNT: usize = 10_000;

#[test]
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
    for artifact in [
        "filtered_all.fq",
        "filtered.fq",
        "unitigs.fa",
        "unitigs.gfa",
    ] {
        let path = output_dir.join(artifact);
        assert!(
            path.is_file(),
            "assembler did not create {}",
            path.display()
        );
    }
    let unitigs =
        fs::read_to_string(output_dir.join("unitigs.fa")).expect("read primary assembly FASTA");
    assert!(
        unitigs.lines().any(|line| line.starts_with(">unitig_")),
        "assembler completed without producing a unitig"
    );
    assert!(
        fs::metadata(output_dir.join("filtered_all.fq"))
            .expect("read filtered FASTQ metadata")
            .len()
            > 0,
        "assembler filtered out every sampled read"
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
