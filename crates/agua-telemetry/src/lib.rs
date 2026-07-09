//! agua-telemetry — rastreabilidade total de cada run + IO por modalidade.
//!
//! Todo experimento produz `runs/<exp>/run-<ts>-s<seed>/` com:
//!   manifest.json   — quem/quando/config/seed/git (reprodutível)
//!   config.toml     — cópia exata da config usada
//!   metrics.jsonl   — série temporal (1 linha JSON por log)
//!   summary.json    — resultado final + baselines + wall time
//!   population.json — a codificação em si (os seres vivos)
//!   snapshots/*     — recon/resíduo (PNG imagem/vídeo, PNG+WAV áudio)
//! E appenda 1 linha honesta no JOURNAL.md do workspace.

use agua_core::Field;
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

    pub fn snapshot_path(&self, name: &str) -> PathBuf {
        self.path.join("snapshots").join(name)
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

    /// Campo 2D [0,1] → PNG (pra alvo/reconstrução com DC re-somado).
    pub fn save_field_png(&self, name: &str, f: &Field) -> std::io::Result<()> {
        save_png(&self.snapshot_path(name), f, 0.0, 1.0)
    }

    /// Campo 2D com sinais (resíduo): mapeia [-half, +half] → [0,255], 0 = cinza.
    pub fn save_signed_png(&self, name: &str, f: &Field, half: f32) -> std::io::Result<()> {
        save_png(&self.snapshot_path(name), f, -half, half)
    }
}

fn save_png(path: &Path, f: &Field, lo: f32, hi: f32) -> std::io::Result<()> {
    assert_eq!(f.dims.len(), 2, "save_png é pra campo 2D [h,w]");
    let (h, w) = (f.dims[0], f.dims[1]);
    let span = (hi - lo).max(1e-6);
    let raw: Vec<u8> = f
        .data
        .iter()
        .map(|v| (((v - lo) / span).clamp(0.0, 1.0) * 255.0).round() as u8)
        .collect();
    let img = image::GrayImage::from_raw(w as u32, h as u32, raw).expect("dims");
    img.save(path)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
}

// ─────────────────────────── áudio: WAV + waveform ───────────────────────────

/// Escreve WAV PCM16 mono. `samples` em [-1,1] (clamp).
pub fn save_wav(path: &Path, samples: &[f32], sample_rate: u32) -> std::io::Result<()> {
    let n = samples.len() as u32;
    let data_len = n * 2;
    let mut out = Vec::with_capacity(44 + data_len as usize);
    out.extend_from_slice(b"RIFF");
    out.extend_from_slice(&(36 + data_len).to_le_bytes());
    out.extend_from_slice(b"WAVEfmt ");
    out.extend_from_slice(&16u32.to_le_bytes());
    out.extend_from_slice(&1u16.to_le_bytes()); // PCM
    out.extend_from_slice(&1u16.to_le_bytes()); // mono
    out.extend_from_slice(&sample_rate.to_le_bytes());
    out.extend_from_slice(&(sample_rate * 2).to_le_bytes()); // byte rate
    out.extend_from_slice(&2u16.to_le_bytes()); // block align
    out.extend_from_slice(&16u16.to_le_bytes()); // bits
    out.extend_from_slice(b"data");
    out.extend_from_slice(&data_len.to_le_bytes());
    for s in samples {
        let v = (s.clamp(-1.0, 1.0) * 32767.0).round() as i16;
        out.extend_from_slice(&v.to_le_bytes());
    }
    fs::write(path, out)
}

/// Lê WAV PCM16 (mono ou multi-canal → média dos canais) → ([-1,1], sample_rate).
pub fn load_wav(path: &Path) -> std::io::Result<(Vec<f32>, u32)> {
    let b = fs::read(path)?;
    let err = |m: &str| std::io::Error::new(std::io::ErrorKind::InvalidData, m.to_string());
    if b.len() < 44 || &b[0..4] != b"RIFF" || &b[8..12] != b"WAVE" {
        return Err(err("não é WAV RIFF"));
    }
    let u16le = |s: &[u8]| u16::from_le_bytes([s[0], s[1]]);
    let u32le = |s: &[u8]| u32::from_le_bytes([s[0], s[1], s[2], s[3]]);
    let mut pos = 12usize;
    let mut channels = 1u16;
    let mut sample_rate = 0u32;
    let mut bits = 0u16;
    let mut data: Option<(usize, usize)> = None;
    while pos + 8 <= b.len() {
        let id = &b[pos..pos + 4];
        let sz = u32le(&b[pos + 4..pos + 8]) as usize;
        let body_start = pos + 8;
        let body_end = (body_start + sz).min(b.len());
        match id {
            b"fmt " if sz >= 16 => {
                channels = u16le(&b[body_start + 2..body_start + 4]).max(1);
                sample_rate = u32le(&b[body_start + 4..body_start + 8]);
                bits = u16le(&b[body_start + 14..body_start + 16]);
            }
            b"data" => data = Some((body_start, body_end)),
            _ => {}
        }
        pos = body_start + sz + (sz & 1);
    }
    if bits != 16 {
        return Err(err("só PCM16 suportado (converta: ffmpeg -i in -acodec pcm_s16le out.wav)"));
    }
    let (ds, de) = data.ok_or_else(|| err("sem chunk data"))?;
    let frame = 2 * channels as usize;
    let n = (de - ds) / frame;
    let mut samples = Vec::with_capacity(n);
    for i in 0..n {
        let mut acc = 0.0f32;
        for c in 0..channels as usize {
            let o = ds + i * frame + c * 2;
            acc += i16::from_le_bytes([b[o], b[o + 1]]) as f32 / 32768.0;
        }
        samples.push(acc / channels as f32);
    }
    Ok((samples, sample_rate))
}

/// Desenha a forma de onda (min-max por coluna) num PNG w×h — visibilidade
/// de sinais 1D nos snapshots.
pub fn waveform_png(path: &Path, samples: &[f32], w: usize, h: usize) -> std::io::Result<()> {
    let mut img = image::GrayImage::from_pixel(w as u32, h as u32, image::Luma([16u8]));
    let n = samples.len().max(1);
    let peak = samples.iter().fold(1e-6f32, |m, v| m.max(v.abs()));
    for x in 0..w {
        let a = x * n / w;
        let b = (((x + 1) * n / w).max(a + 1)).min(n);
        let (mut lo, mut hi) = (f32::MAX, f32::MIN);
        for v in &samples[a..b] {
            lo = lo.min(*v);
            hi = hi.max(*v);
        }
        let to_row = |v: f32| -> usize {
            let norm = (v / peak).clamp(-1.0, 1.0);
            (((1.0 - norm) * 0.5) * (h - 1) as f32).round() as usize
        };
        let (r0, r1) = (to_row(hi), to_row(lo));
        for y in r0..=r1.min(h - 1) {
            img.put_pixel(x as u32, y as u32, image::Luma([230u8]));
        }
    }
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
        let f = Field::from_fn_2d(8, 8, |x, _| x as f32 / 8.0);
        rd.save_field_png("t.png", &f).unwrap();
        assert!(rd.path.join("manifest.json").exists());
        assert!(rd.path.join("metrics.jsonl").exists());
        assert!(rd.path.join("summary.json").exists());
        assert!(rd.path.join("snapshots/t.png").exists());
        journal_append(&tmp, "| x | unit | 1 | 1 | 10 | 0 | 0 | - | - | teste |").unwrap();
        assert!(tmp.join("JOURNAL.md").exists());
        fs::remove_dir_all(&tmp).ok();
    }

    #[test]
    fn wav_roundtrip() {
        let tmp = std::env::temp_dir().join(format!("agua-wav-{}.wav", unix_now()));
        let samples: Vec<f32> = (0..1000)
            .map(|i| 0.5 * (std::f32::consts::TAU * 0.01 * i as f32).sin())
            .collect();
        save_wav(&tmp, &samples, 16000).unwrap();
        let (back, sr) = load_wav(&tmp).unwrap();
        assert_eq!(sr, 16000);
        assert_eq!(back.len(), 1000);
        let max_err = samples
            .iter()
            .zip(&back)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f32, f32::max);
        assert!(max_err < 1e-3, "erro quantização PCM16: {max_err}");
        fs::remove_file(&tmp).ok();
    }

    #[test]
    fn waveform_png_writes() {
        let tmp = std::env::temp_dir().join(format!("agua-wf-{}.png", unix_now()));
        let samples: Vec<f32> = (0..500)
            .map(|i| (std::f32::consts::TAU * 0.02 * i as f32).sin())
            .collect();
        waveform_png(&tmp, &samples, 256, 96).unwrap();
        assert!(tmp.exists());
        fs::remove_file(&tmp).ok();
    }
}
