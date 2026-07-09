//! agua-telemetry — rastreabilidade total de cada run.
//!
//! Todo experimento produz `runs/<exp>/run-<ts>-s<seed>/` com:
//!   manifest.json   — quem/quando/config/seed/git (reprodutível)
//!   config.toml     — cópia exata da config usada
//!   metrics.jsonl   — série temporal (1 linha JSON por log)
//!   summary.json    — resultado final + baselines + wall time
//!   population.json — a codificação em si (os seres vivos)
//!   snapshots/*.png — recon/resíduo ao longo do tempo
//! E appenda 1 linha honesta no JOURNAL.md do workspace.

use agua_core::Field2D;
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub fn unix_now() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

/// Hash do git HEAD do diretório, se for repo (traça código ↔ resultado).
pub fn git_hash(dir: &Path) -> Option<String> {
    let out = std::process::Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .current_dir(dir)
        .output()
        .ok()?;
    if out.status.success() {
        Some(String::from_utf8_lossy(&out.stdout).trim().to_string())
    } else {
        None
    }
}

pub struct RunDir {
    pub path: PathBuf,
    metrics: BufWriter<File>,
}

impl RunDir {
    pub fn create(root: &Path, exp: &str, seed: u64) -> std::io::Result<Self> {
        let path = root
            .join("runs")
            .join(exp)
            .join(format!("run-{}-s{}", unix_now(), seed));
        fs::create_dir_all(path.join("snapshots"))?;
        let metrics = BufWriter::new(File::create(path.join("metrics.jsonl"))?);
        Ok(RunDir { path, metrics })
    }

    pub fn write_manifest(&self, manifest: &serde_json::Value) -> std::io::Result<()> {
        fs::write(
            self.path.join("manifest.json"),
            serde_json::to_string_pretty(manifest)?,
        )
    }

    pub fn copy_config(&self, config_toml: &str) -> std::io::Result<()> {
        fs::write(self.path.join("config.toml"), config_toml)
    }

    /// Uma linha da série temporal.
    pub fn log(&mut self, v: &serde_json::Value) -> std::io::Result<()> {
        writeln!(self.metrics, "{v}")?;
        self.metrics.flush()
    }

    pub fn write_summary(&self, v: &serde_json::Value) -> std::io::Result<()> {
        fs::write(
            self.path.join("summary.json"),
            serde_json::to_string_pretty(v)?,
        )
    }

    pub fn write_population(&self, v: &serde_json::Value) -> std::io::Result<()> {
        fs::write(
            self.path.join("population.json"),
            serde_json::to_string_pretty(v)?,
        )
    }

    /// Campo [0,1] → PNG (pra alvo/reconstrução com DC re-somado).
    pub fn save_field_png(&self, name: &str, f: &Field2D) -> std::io::Result<()> {
        save_png(&self.path.join("snapshots").join(name), f, 0.0, 1.0)
    }

    /// Campo com sinais (resíduo): mapeia [-half, +half] → [0,255], 0 = cinza.
    pub fn save_signed_png(&self, name: &str, f: &Field2D, half: f32) -> std::io::Result<()> {
        save_png(&self.path.join("snapshots").join(name), f, -half, half)
    }
}

fn save_png(path: &Path, f: &Field2D, lo: f32, hi: f32) -> std::io::Result<()> {
    let span = (hi - lo).max(1e-6);
    let raw: Vec<u8> = f
        .data
        .iter()
        .map(|v| (((v - lo) / span).clamp(0.0, 1.0) * 255.0).round() as u8)
        .collect();
    let img = image::GrayImage::from_raw(f.w as u32, f.h as u32, raw).expect("dims");
    img.save(path)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
}

/// Uma linha por run no diário do workspace — visibilidade sem abrir JSON.
pub fn journal_append(root: &Path, line: &str) -> std::io::Result<()> {
    let path = root.join("JOURNAL.md");
    if !path.exists() {
        fs::write(
            &path,
            "# JOURNAL — iara-agua (1 linha por run, honesto)\n\n\
             | quando | exp | seed | ticks | psnr | vivos | bytes | vs MP | vs random | nota |\n\
             |---|---|---|---|---|---|---|---|---|---|\n",
        )?;
    }
    let mut f = OpenOptions::new().append(true).open(&path)?;
    writeln!(f, "{line}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rundir_creates_and_writes() {
        let tmp = std::env::temp_dir().join(format!("agua-test-{}", unix_now()));
        let mut rd = RunDir::create(&tmp, "unit", 1).unwrap();
        rd.write_manifest(&serde_json::json!({"seed": 1})).unwrap();
        rd.log(&serde_json::json!({"tick": 1, "psnr": 10.0})).unwrap();
        rd.write_summary(&serde_json::json!({"ok": true})).unwrap();
        let f = Field2D::from_fn(8, 8, |x, _| x as f32 / 8.0);
        rd.save_field_png("t.png", &f).unwrap();
        assert!(rd.path.join("manifest.json").exists());
        assert!(rd.path.join("metrics.jsonl").exists());
        assert!(rd.path.join("summary.json").exists());
        assert!(rd.path.join("snapshots/t.png").exists());
        journal_append(&tmp, "| x | unit | 1 | 1 | 10 | 0 | 0 | - | - | teste |").unwrap();
        assert!(tmp.join("JOURNAL.md").exists());
        fs::remove_dir_all(&tmp).ok();
    }
}
