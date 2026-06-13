use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::{Command, Stdio};

use crate::DynResult;

pub(crate) fn extract_interval_pair_with_samtools(
    source_url: &str,
    region: &str,
    r1_all: &Path,
    r2_all: &Path,
    other: &Path,
    singleton: &Path,
) -> DynResult<bool> {
    let mut view = Command::new("samtools")
        .args(["view", "-u", "-f", "3", "-F", "2304", source_url, region])
        .stdout(Stdio::piped())
        .spawn()?;
    let view_stdout = view
        .stdout
        .take()
        .ok_or("samtools view stdout pipe unavailable")?;
    let mut collate = Command::new("samtools")
        .args(["collate", "-u", "-O", "-"])
        .stdin(Stdio::from(view_stdout))
        .stdout(Stdio::piped())
        .spawn()?;
    let collate_stdout = collate
        .stdout
        .take()
        .ok_or("samtools collate stdout pipe unavailable")?;
    let fastq_status = Command::new("samtools")
        .arg("fastq")
        .arg("-n")
        .arg("-1")
        .arg(r1_all)
        .arg("-2")
        .arg(r2_all)
        .arg("-0")
        .arg(other)
        .arg("-s")
        .arg(singleton)
        .arg("-")
        .stdin(Stdio::from(collate_stdout))
        .status()?;
    let collate_status = collate.wait()?;
    let view_status = view.wait()?;
    Ok(view_status.success() && collate_status.success() && fastq_status.success())
}

pub(crate) fn curl_to_file(url: &str, out_path: &Path) -> DynResult<bool> {
    let out = fs::File::create(out_path)?;
    let status = Command::new("curl")
        .args(["-fsSL", url])
        .stdout(Stdio::from(out))
        .status()?;
    Ok(status.success())
}

pub(crate) fn curl_gzip_to_file(url: &str, out_path: &Path) -> DynResult<bool> {
    let mut curl = Command::new("curl")
        .args(["-fsSL", url])
        .stdout(Stdio::piped())
        .spawn()?;
    let curl_stdout = curl.stdout.take().ok_or("curl stdout pipe unavailable")?;
    let out = fs::File::create(out_path)?;
    let mut gzip = Command::new("gzip")
        .arg("-dc")
        .stdin(Stdio::from(curl_stdout))
        .stdout(Stdio::from(out))
        .spawn()?;
    let gzip_status = gzip.wait()?;
    let curl_status = curl.wait()?;
    Ok(curl_status.success() && gzip_status.success())
}

pub(crate) fn curl_gzip_prefix_lines_to_file(
    url: &str,
    lines: usize,
    out_path: &Path,
) -> DynResult<bool> {
    let mut curl = Command::new("curl")
        .args(["-fsSL", url])
        .stdout(Stdio::piped())
        .spawn()?;
    let curl_stdout = curl.stdout.take().ok_or("curl stdout pipe unavailable")?;
    let mut gzip = Command::new("gzip")
        .arg("-dc")
        .stdin(Stdio::from(curl_stdout))
        .stdout(Stdio::piped())
        .spawn()?;
    let gzip_stdout = gzip.stdout.take().ok_or("gzip stdout pipe unavailable")?;
    let mut reader = BufReader::new(gzip_stdout);
    let mut out = fs::File::create(out_path)?;
    let mut copied = 0usize;
    let mut line = Vec::new();
    while copied < lines {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        out.write_all(&line)?;
        copied += 1;
    }
    out.flush()?;
    drop(reader);
    let complete = copied == lines;
    if complete {
        let _ = gzip.kill();
        let _ = curl.kill();
    }
    let gzip_status = gzip.wait()?;
    let curl_status = curl.wait()?;
    Ok(complete || (curl_status.success() && gzip_status.success() && copied == lines))
}
