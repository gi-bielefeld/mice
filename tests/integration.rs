use mice::compression::*;
use std::error::Error;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

#[test]
fn compare_all_gfas_random() -> Result<(), Box<dyn Error>> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let gfa_dir = root.join("tests/data/random/gfa");
    let out_dir = root.join("tests/data/random/size");
    compare_gfa_dir(&gfa_dir, &out_dir)
}

#[test]
fn compare_all_gfas_real() -> Result<(), Box<dyn Error>> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let gfa_dir = root.join("tests/data/real/gfa");
    let out_dir = root.join("tests/data/real/size");
    compare_gfa_dir(&gfa_dir, &out_dir)
}

fn compare_lines(stem: &str, got: &[String], expected: &[String]) -> Vec<String> {
    let mut failures: Vec<String> = Vec::new();
    let mut diff_count = 0usize;
    let max_len = got.len().max(expected.len());

    for i in 0..max_len {
        let miss = "<MISSING>".to_string();
        let g = got.get(i).unwrap_or(&miss);
        let e = expected.get(i).unwrap_or(&miss);

        if g != e {
            diff_count += 1;
            if diff_count <= 10 {
                failures.push(format!(
                    "[DIFF] {}: line {}: got {:?}  !=  expected {:?}",
                    stem,
                    i + 1,
                    g,
                    e
                ));
            }
        }
    }

    if diff_count > 10 {
        failures.push(format!(
            "[DIFF] {}: +{} more differing lines omitted",
            stem,
            diff_count - 10
        ));
    }

    failures
}

fn compare_gfa_dir(gfa_dir: &Path, out_dir: &Path) -> Result<(), Box<dyn Error>> {
    let gfa_files: Vec<PathBuf> = fs::read_dir(gfa_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map(|e| e == "gfa").unwrap_or(false))
        .collect();
    //gfa_files.sort();

    if gfa_files.is_empty() {
        panic!("No .gfa files found in {}", gfa_dir.display());
    }

    let mut failures: Vec<String> = Vec::new();

    for gfa_path in gfa_files {
        let stem = gfa_path
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .to_string();
        let out_path = out_dir.join(format!("{}.txt", stem));

        // Run mice
        let num_parts = run_mice_test(gfa_path.to_str().unwrap(), Some("gfa")).unwrap();
        let got = vec![num_parts.to_string()];

        // Load expected
        let expected = match read_expected_lines(&out_path) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!(
                    "[ERROR] {}: cannot read expected file {}: {}",
                    gfa_path.display(),
                    out_path.display(),
                    e
                ));
                continue;
            }
        };

        failures.extend(compare_lines(&stem, &got, &expected));
    }

    if failures.is_empty() {
        Ok(())
    } else {
        panic!("\n{}\n", failures.join("\n"));
    }
}

fn read_expected_lines(path: &Path) -> Result<Vec<String>, Box<dyn Error>> {
    let f = fs::File::open(path)?;
    let rdr = BufReader::new(f);
    let mut lines: Vec<String> = Vec::new();
    for line in rdr.lines() {
        let mut s = line?;
        if s.ends_with('\r') {
            s.pop(); // normalize CRLF
        }
        lines.push(s);
    }
    Ok(lines)
}

