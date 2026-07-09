//! agua-lab — runner de experimentos do solver-por-sobrevivência.
//!
//! Uso:
//!   agua-lab run --config configs/exp001_image_survival.toml [--out .] [--seed N] [--ticks N]
//!   agua-lab demo                        # smoke rápido com alvo sintético
//!
//! Cada run gera runs/<exp>/run-<ts>-s<seed>/ (manifest, metrics.jsonl,
//! snapshots, population.json, summary.json) e appenda 1 linha no JOURNAL.md.

use agua_atoms::AtomBounds;
use agua_life::{LifeConfig, World};
use agua_metrics::{jpeg_curve, matching_pursuit, random_atoms, RdPoint};
use agua_telemetry::{git_hash, journal_append, unix_now, RunDir};
use serde::Deserialize;
use std::path::Path;

#[derive(Deserialize, Clone)]
struct Config {
    experiment: ExperimentCfg,
    world: WorldCfg,
    life: LifeConfig,
    atoms: AtomBounds,
}

#[derive(Deserialize, Clone)]
struct ExperimentCfg {
    name: String,
    seeds: Vec<u64>,
    max_ticks: u64,
    log_every: u64,
    snapshot_every: u64,
}

#[derive(Deserialize, Clone)]
struct WorldCfg {
    /// "synthetic" ou caminho de um PNG (vira grayscale, redimensiona pra size)
    target: String,
    size: usize,
}

/// Bandas de frequência (ciclos/px) pro histograma da pirâmide.
const BAND_EDGES: [f32; 4] = [0.03, 0.06, 0.12, 0.24];

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("run") => {
            let config_path = flag_value(&args, "--config").unwrap_or_else(|| {
                eprintln!("faltou --config <path>");
                std::process::exit(2);
            });
            let out = flag_value(&args, "--out").unwrap_or_else(|| ".".into());
            let toml_src = std::fs::read_to_string(&config_path)
                .unwrap_or_else(|e| panic!("lendo {config_path}: {e}"));
            let mut cfg: Config =
                toml::from_str(&toml_src).unwrap_or_else(|e| panic!("config inválida: {e}"));
            if let Some(s) = flag_value(&args, "--seed") {
                cfg.experiment.seeds = vec![s.parse().expect("--seed inteiro")];
            }
            if let Some(t) = flag_value(&args, "--ticks") {
                cfg.experiment.max_ticks = t.parse().expect("--ticks inteiro");
            }
            for seed in cfg.experiment.seeds.clone() {
                run_one(&cfg, &toml_src, seed, Path::new(&out));
            }
        }
        Some("demo") => {
            let cfg = demo_config();
            let toml_src = "# demo builtin".to_string();
            run_one(&cfg, &toml_src, 42, Path::new("."));
        }
        _ => {
            eprintln!("uso: agua-lab run --config <toml> [--out dir] [--seed N] [--ticks N] | agua-lab demo");
            std::process::exit(2);
        }
    }
}

fn flag_value(args: &[String], flag: &str) -> Option<String> {
    args.iter()
        .position(|a| a == flag)
        .and_then(|i| args.get(i + 1).cloned())
}

fn demo_config() -> Config {
    Config {
        experiment: ExperimentCfg {
            name: "demo".into(),
            seeds: vec![42],
            max_ticks: 800,
            log_every: 50,
            snapshot_every: 200,
        },
        world: WorldCfg { target: "synthetic".into(), size: 64 },
        life: LifeConfig::default(),
        atoms: AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 },
    }
}

/// Alvo sintético multi-escala determinístico: gradiente suave (baixa freq) +
/// blob gaussiano (média) + listras finas num quadrante (alta freq) + quadrado
/// de borda dura. Testa exatamente a claim "a pirâmide emerge".
fn synthetic_target(size: usize) -> agua_core::Field2D {
    let s = size as f32;
    agua_core::Field2D::from_fn(size, size, |x, y| {
        let (xf, yf) = (x as f32, y as f32);
        let grad = 0.25 * (std::f32::consts::TAU * xf / s * 0.5).sin()
            * (std::f32::consts::TAU * yf / s * 0.5).cos();
        let blob = 0.45
            * (-((xf - s * 0.3).powi(2) + (yf - s * 0.32).powi(2))
                / (2.0 * (s * 0.12).powi(2)))
            .exp();
        let stripes = if xf > s * 0.55 && yf > s * 0.55 {
            0.22 * (std::f32::consts::TAU * 0.22 * xf).sin()
        } else {
            0.0
        };
        let square = if (xf - s * 0.7).abs() < s * 0.09 && (yf - s * 0.25).abs() < s * 0.09 {
            0.3
        } else {
            0.0
        };
        (0.45 + grad + blob + stripes + square).clamp(0.0, 1.0)
    })
}

fn load_target(wcfg: &WorldCfg) -> agua_core::Field2D {
    if wcfg.target == "synthetic" {
        return synthetic_target(wcfg.size);
    }
    let img = image::open(&wcfg.target)
        .unwrap_or_else(|e| panic!("abrindo alvo {}: {e}", wcfg.target))
        .resize_exact(
            wcfg.size as u32,
            wcfg.size as u32,
            image::imageops::FilterType::CatmullRom,
        )
        .to_luma8();
    agua_core::Field2D {
        w: wcfg.size,
        h: wcfg.size,
        data: img.as_raw().iter().map(|v| *v as f32 / 255.0).collect(),
    }
}

fn rd_json(p: &RdPoint) -> serde_json::Value {
    serde_json::json!({"label": p.label, "atoms": p.atoms, "bytes": p.bytes, "psnr_db": p.psnr_db})
}

fn run_one(cfg: &Config, toml_src: &str, seed: u64, out_root: &Path) {
    let exp = &cfg.experiment;
    let t_start = std::time::Instant::now();

    // alvo com DC (pra JPEG/snapshots) e versão zero-média (pro mundo)
    let target_dc = load_target(&cfg.world);
    let dc = target_dc.mean();
    let mut target = target_dc.clone();
    for v in target.data.iter_mut() {
        *v -= dc;
    }

    let mut rd = RunDir::create(out_root, &exp.name, seed).expect("criando run dir");
    rd.write_manifest(&serde_json::json!({
        "exp": exp.name,
        "seed": seed,
        "created_unix": unix_now(),
        "git": git_hash(out_root),
        "world": {"target": cfg.world.target, "size": cfg.world.size, "dc": dc},
        "max_ticks": exp.max_ticks,
    }))
    .unwrap();
    rd.copy_config(toml_src).unwrap();
    rd.save_field_png("target.png", &target_dc).unwrap();

    let mut world = World::new(target.clone(), cfg.life.clone(), cfg.atoms.clone(), seed);
    let mut best_psnr = world.psnr(1.0);

    println!(
        "▶ {} seed={} | alvo {}x{} | pop inicial {} | psnr0 {:.2} dB",
        exp.name, seed, cfg.world.size, cfg.world.size, world.beings.len(), best_psnr
    );

    for tick in 1..=exp.max_ticks {
        world.step();
        if tick % exp.log_every == 0 {
            let p = world.psnr(1.0);
            best_psnr = best_psnr.max(p);
            let bands = world.freq_bands(&BAND_EDGES);
            rd.log(&serde_json::json!({
                "tick": tick,
                "psnr_db": p,
                "alive": world.beings.len(),
                "births": world.births,
                "deaths": world.deaths,
                "bytes": world.bytes(),
                "residual_energy": world.residual.energy(),
                "freq_bands": bands,
            }))
            .unwrap();
        }
        if tick % exp.snapshot_every == 0 {
            let mut recon_dc = world.recon.clone();
            for v in recon_dc.data.iter_mut() {
                *v += dc;
            }
            rd.save_field_png(&format!("recon_{tick:06}.png"), &recon_dc).unwrap();
            rd.save_signed_png(&format!("residual_{tick:06}.png"), &world.residual, 0.5)
                .unwrap();
        }
    }

    world.rebuild(); // números finais sem deriva float
    let final_psnr = world.psnr(1.0);
    let alive = world.beings.len();
    let bytes = world.bytes();
    let wall = t_start.elapsed().as_secs_f32();

    // ---- baselines honestos com o MESMO número de átomos ----
    let n = alive.max(1);
    let cps = [n / 4, n / 2, n]
        .iter()
        .copied()
        .filter(|v| *v > 0)
        .collect::<Vec<_>>();
    let rnd = random_atoms(&target, n, &cfg.atoms, seed ^ 0xABCD, &cps);
    let mp = matching_pursuit(&target, n, 24, &cfg.atoms, seed ^ 0x1234, &cps);
    let jpeg = jpeg_curve(&target_dc, &[10, 30, 50, 70, 90]);
    let rnd_final = rnd.last().map(|p| p.psnr_db).unwrap_or(0.0);
    let mp_final = mp.last().map(|p| p.psnr_db).unwrap_or(0.0);

    rd.write_population(&serde_json::to_value(world.population_records()).unwrap())
        .unwrap();
    rd.write_summary(&serde_json::json!({
        "final_psnr_db": final_psnr,
        "best_psnr_db": best_psnr,
        "alive": alive,
        "bytes": bytes,
        "births": world.births,
        "deaths": world.deaths,
        "ticks": exp.max_ticks,
        "wall_secs": wall,
        "invariant_error": world.invariant_error(),
        "baselines": {
            "random": rnd.iter().map(rd_json).collect::<Vec<_>>(),
            "matching_pursuit": mp.iter().map(rd_json).collect::<Vec<_>>(),
            "jpeg": jpeg.iter().map(rd_json).collect::<Vec<_>>(),
        },
    }))
    .unwrap();

    // snapshot final
    let mut recon_dc = world.recon.clone();
    for v in recon_dc.data.iter_mut() {
        *v += dc;
    }
    rd.save_field_png("recon_final.png", &recon_dc).unwrap();
    rd.save_signed_png("residual_final.png", &world.residual, 0.5).unwrap();

    println!(
        "✔ {} seed={} | {:.1}s | psnr {:.2} dB ({} seres, {} B) | random {:.2} | MP {:.2} | jpeg: {}",
        exp.name,
        seed,
        wall,
        final_psnr,
        alive,
        bytes,
        rnd_final,
        mp_final,
        jpeg.iter()
            .map(|p| format!("q?{}B={:.1}dB", p.bytes, p.psnr_db))
            .collect::<Vec<_>>()
            .join(" ")
    );

    let verdict = if final_psnr > rnd_final + 0.5 {
        format!("bate random (+{:.2} dB), gap MP {:.2} dB", final_psnr - rnd_final, mp_final - final_psnr)
    } else {
        "NÃO bate random — ruído".to_string()
    };
    journal_append(
        out_root,
        &format!(
            "| {} | {} | {} | {} | {:.2} | {} | {} | {:.2} | {:.2} | {} |",
            unix_now(),
            exp.name,
            seed,
            exp.max_ticks,
            final_psnr,
            alive,
            bytes,
            mp_final,
            rnd_final,
            verdict
        ),
    )
    .unwrap();
}
