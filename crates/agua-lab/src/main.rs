//! agua-lab — runner de experimentos do solver-por-sobrevivência.
//!
//! Uso:
//!   agua-lab run --config configs/exp001_image_survival.toml [--out .] [--seed N] [--ticks N]
//!   agua-lab demo          # smoke imagem sintética
//!   agua-lab demo-audio    # smoke áudio sintético
//!   agua-lab demo-video    # smoke vídeo sintético
//!
//! Modalidades ([world].modality): "image" (PNG/sintético), "audio"
//! (WAV PCM16/sintético), "video" (sintético). A física é a MESMA
//! (`World<A: Atom>`); muda só o átomo e o IO.
//!
//! Cada run gera runs/<exp>/run-<ts>-s<seed>/ (manifest, metrics.jsonl,
//! snapshots, population.json, summary.json) e appenda 1 linha no JOURNAL.md.

use agua_atoms::{Atom, AtomBounds, Gabor1D, Gabor2D, Gabor3D};
use agua_core::Field;
use agua_life::{LifeConfig, World};
use agua_metrics::{jpeg_curve, matching_pursuit, random_atoms, RdPoint};
use agua_telemetry::{
    git_hash, journal_append, load_wav, save_wav, unix_now, waveform_png, RunDir,
};
use serde::Deserialize;
use std::path::Path;

mod compare;
mod denoise;
mod law;
mod law2;
mod law3;
mod law4;
mod law5;
mod law6;
mod law7;
mod law8;
mod law9;
mod generate;
mod occlusion;
mod probe;
mod resonance;
mod targets;
use targets::{synthetic_audio, synthetic_image, synthetic_video};

#[derive(Deserialize, Clone)]
struct Config {
    experiment: ExperimentCfg,
    world: WorldCfg,
    life: LifeConfig,
    atoms: AtomBounds,
}

fn default_mp_candidates() -> usize {
    24
}

#[derive(Deserialize, Clone)]
struct ExperimentCfg {
    name: String,
    seeds: Vec<u64>,
    max_ticks: u64,
    log_every: u64,
    snapshot_every: u64,
    /// candidatos por passo do matching pursuit (teto guloso). Mais = teto mais duro.
    /// 0 = ISO-BUDGET: dá ao MP o mesmo nº de avaliações que a sobrevivência gastou.
    #[serde(default = "default_mp_candidates")]
    mp_candidates: usize,
}

fn default_modality() -> String {
    "image".into()
}
fn default_frames() -> usize {
    16
}
fn default_sample_rate() -> u32 {
    16000
}

#[derive(Deserialize, Clone)]
struct WorldCfg {
    /// "image" | "audio" | "video"
    #[serde(default = "default_modality")]
    modality: String,
    /// "synthetic" ou caminho (.png pra image, .wav PCM16 pra audio)
    target: String,
    /// image/video: lado; audio: nº de amostras
    size: usize,
    /// video: nº de frames
    #[serde(default = "default_frames")]
    frames: usize,
    /// audio: taxa pro WAV exportado (e truncagem do carregado)
    #[serde(default = "default_sample_rate")]
    sample_rate: u32,
    /// ADAPTAÇÃO: segundo alvo. Se presente + switch_at>0, o mundo TROCA de
    /// sinal no tick switch_at (a população tem que se re-harmonizar).
    #[serde(default)]
    target_b: Option<String>,
    #[serde(default)]
    switch_at: u64,
}

/// Bandas de frequência (ciclos/amostra) pro histograma da pirâmide.
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
                dispatch(&cfg, &toml_src, seed, Path::new(&out));
            }
        }
        Some("compare") => {
            let target = flag_value(&args, "--target").unwrap_or_else(|| "synthetic".into());
            let out = flag_value(&args, "--out").unwrap_or_else(|| "/tmp/agua-compare".into());
            let seed: u64 = flag_value(&args, "--seed").map_or(11, |v| v.parse().unwrap());
            let ticks: u64 = flag_value(&args, "--ticks").map_or(3000, |v| v.parse().unwrap());
            let met: f32 = flag_value(&args, "--metabolism").map_or(0.8, |v| v.parse().unwrap());
            let size: usize = flag_value(&args, "--size").map_or(64, |v| v.parse().unwrap());
            let wcfg = WorldCfg {
                modality: "image".into(),
                target,
                size,
                frames: 0,
                sample_rate: 16000,
                target_b: None,
                switch_at: 0,
            };
            let mut life = LifeConfig::default();
            life.cap = 3000;
            life.metabolism = met;
            life.rain_every = 1;
            life.rain_count = 6; // rainHeavy: melhor variante da bateria F1
            compare::run(compare::CompareArgs {
                target_dc: load_image_target(&wcfg),
                life,
                bounds: AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 },
                ticks,
                seed,
                out: Path::new(&out),
                checkpoints: &[10, 50, 200, 800, 3000],
            });
        }
        Some("denoise") => {
            let target = flag_value(&args, "--target").unwrap_or_else(|| "synthetic".into());
            let modality = flag_value(&args, "--modality").unwrap_or_else(|| "image".into());
            let seed: u64 = flag_value(&args, "--seed").map_or(11, |v| v.parse().unwrap());
            let ticks: u64 = flag_value(&args, "--ticks").map_or(1500, |v| v.parse().unwrap());
            let sigma: f32 = flag_value(&args, "--sigma").map_or(0.1, |v| v.parse().unwrap());
            let mut met: f32 = flag_value(&args, "--metabolism").map_or(0.8, |v| v.parse().unwrap());
            let size: usize = flag_value(&args, "--size").map_or(96, |v| v.parse().unwrap());
            let grad: bool = flag_value(&args, "--gradient").map_or(true, |v| v == "true");
            let out = flag_value(&args, "--out");
            // --auto: aplica a LEI met = k·σ̂², com σ̂ estimado do sinal ruidoso.
            // Nenhum hiperparâmetro, nenhum acesso ao sinal limpo.
            let auto: bool = args.iter().any(|a| a == "--auto");
            let k_law: f32 = flag_value(&args, "--k").map_or(40.0, |v| v.parse().unwrap());
            let wcfg = WorldCfg {
                modality: modality.clone(),
                target,
                size,
                frames: 0,
                sample_rate: 16000,
                target_b: None,
                switch_at: 0,
            };
            let mut life = LifeConfig::default();
            life.cap = 3000;
            life.metabolism = met;
            life.rain_every = 1;
            life.rain_count = 6;
            life.rain_competition = 16;
            life.tune_every = 1;
            life.gradient = grad;
            let bounds = if modality == "audio" {
                AtomBounds { freq_min: 0.002, freq_max: 0.45, sigma_min: 2.0, sigma_max: 400.0 }
            } else {
                AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 }
            };
            if auto {
                // estima σ do MESMO sinal ruidoso que o solver vai ver
                let clean_probe = if modality == "audio" { load_audio_target(&wcfg) } else { load_image_target(&wcfg) };
                let noisy_probe = denoise::add_noise(&clean_probe, sigma, seed ^ 0x0157);
                let sig_hat = denoise::estimate_sigma(&noisy_probe);
                met = k_law * sig_hat * sig_hat;
                life.metabolism = met;
                println!("[AUTO] σ real={sigma:.4}  σ̂ estimado={sig_hat:.4}  (erro {:+.1}%)  → met = {k_law}·σ̂² = {met:.4}",
                    100.0 * (sig_hat - sigma) / sigma);
            }
            let (rows, recon, noisy) = if modality == "audio" {
                let clean = load_audio_target(&wcfg);
                denoise::run::<Gabor1D>(denoise::DenoiseArgs {
                    clean, sigma, life, bounds, ticks, seed,
                    mp_candidates: 64, n_max: 1200, blur: false, label: "audio",
                })
            } else {
                let clean = load_image_target(&wcfg);
                denoise::run::<Gabor2D>(denoise::DenoiseArgs {
                    clean, sigma, life, bounds, ticks, seed,
                    mp_candidates: 64, n_max: 1200, blur: true, label: "image",
                })
            };
            println!("\n== DENOISE {modality} sigma={sigma} seed={seed} met={met} grad={grad} ==");
            println!("{:<18} {:>10} {:>8}  {}", "metodo", "psnr_limpo", "atomos", "nota");
            for r in &rows {
                println!("{:<18} {:>10.2} {:>8}  {}", r.metodo, r.psnr_limpo, r.atomos, r.nota);
            }
            println!("{}", serde_json::to_string(&rows).unwrap());
            if let Some(dir) = out {
                let d = Path::new(&dir);
                std::fs::create_dir_all(d).unwrap();
                if modality != "audio" {
                    let save = |name: &str, f: &Field| {
                        let (h, w) = (f.dims[0], f.dims[1]);
                        let raw: Vec<u8> = f.data.iter().map(|v| (v.clamp(0.0,1.0)*255.0).round() as u8).collect();
                        image::GrayImage::from_raw(w as u32, h as u32, raw).unwrap().save(d.join(name)).unwrap();
                    };
                    save("clean.png", &load_image_target(&wcfg));
                    save("noisy.png", &noisy);
                    save("survival.png", &recon);
                }
            }
        }
        Some("law") => {
            let seed: u64 = flag_value(&args, "--seed").map_or(11, |v| v.parse().unwrap());
            let sigmas = [0.02f32, 0.04, 0.08, 0.16, 0.32];
            println!("=== LEI DO METABOLISMO CRÍTICO: met_c = c·σ^α ===");
            println!("teoria: átomo que se ajusta a ruído captura c·σ² de energia espúria.");
            println!("P1: α = 2   P2: c independe do sinal   P3: c cresce com os graus de liberdade\n");

            let b2 = AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
            println!("-- Gabor2D (6 params de geometria), campo 64x64 --");
            let r2d = law::run_for::<Gabor2D>("Gabor2D", &[64, 64], &b2, &sigmas, seed);
            let (a2, k2, q2) = law::fit_exponent(&r2d);
            println!("   ajuste: met_c = {k2:.3}·σ^{a2:.3}   (R²={q2:.4})");

            let b1 = AtomBounds { freq_min: 0.002, freq_max: 0.45, sigma_min: 2.0, sigma_max: 400.0 };
            println!("\n-- Gabor1D (4 params de geometria), campo 4096 --");
            let r1d = law::run_for::<Gabor1D>("Gabor1D", &[4096], &b1, &sigmas, seed);
            let (a1, k1, q1) = law::fit_exponent(&r1d);
            println!("   ajuste: met_c = {k1:.3}·σ^{a1:.3}   (R²={q1:.4})");

            let mean = |rs: &[law::LawRow], f: fn(&law::LawRow) -> f32| -> f32 {
                rs.iter().map(f).sum::<f32>() / rs.len() as f32
            };
            println!("\n=== VEREDITO ===");
            println!("P1 expoente α:  2D {a2:.3}  |  1D {a1:.3}   (predito 2.000)");
            println!("P2 c_crítico:   2D {:.2} ± spread  |  1D {:.2}", mean(&r2d, |r| r.c_critico), mean(&r1d, |r| r.c_critico));
            println!("P3 c_busca:     2D {:.2} (6 geom)  |  1D {:.2} (4 geom)  → 1D deve ser MENOR",
                mean(&r2d, |r| r.c_com_busca), mean(&r1d, |r| r.c_com_busca));
            println!("âncora c_sem_busca (deve ser ≈1): 2D {:.2} | 1D {:.2}",
                mean(&r2d, |r| r.c_sem_busca), mean(&r1d, |r| r.c_sem_busca));
            let all: Vec<_> = r2d.iter().chain(r1d.iter()).collect();
            println!("\n{}", serde_json::to_string(&all).unwrap());
        }
        Some("law2") => {
            let seed: u64 = flag_value(&args, "--seed").map_or(11, |v| v.parse().unwrap());
            law2::run(seed);
        }
        Some("law3") => {
            let seed: u64 = flag_value(&args, "--seed").map_or(11, |v| v.parse().unwrap());
            law3::run(seed);
        }
        Some("law4") => {
            let seed: u64 = flag_value(&args, "--seed").map_or(0, |v| v.parse().unwrap());
            law4::run(seed);
        }
        Some("law5") => {
            let seed: u64 = flag_value(&args, "--seed").map_or(11, |v| v.parse().unwrap());
            law5::run(seed);
        }
        Some("law6") => {
            let seed: u64 = flag_value(&args, "--seed").map_or(11, |v| v.parse().unwrap());
            law6::run(seed);
        }
        Some("law7") => {
            let seed: u64 = flag_value(&args, "--seed").map_or(11, |v| v.parse().unwrap());
            law7::run(seed);
        }
        Some("law8") => {
            let seed: u64 = flag_value(&args, "--seed").map_or(0, |v| v.parse().unwrap());
            law8::run(seed);
        }
        Some("law9") => {
            let seed: u64 = flag_value(&args, "--seed").map_or(11, |v| v.parse().unwrap());
            law9::run(seed);
        }
        Some("probe") => probe::run(&args),
        Some("occlusion") => occlusion::run(&args),
        Some("generate") => generate::run(&args),
        Some("resonance") => resonance::run(&args),
        Some("demo") => demo("image", 800),
        Some("demo-audio") => demo("audio", 800),
        Some("demo-video") => demo("video", 400),
        _ => {
            eprintln!(
                "uso: agua-lab run --config <toml> [--out dir] [--seed N] [--ticks N] | demo | demo-audio | demo-video"
            );
            std::process::exit(2);
        }
    }
}

fn flag_value(args: &[String], flag: &str) -> Option<String> {
    args.iter()
        .position(|a| a == flag)
        .and_then(|i| args.get(i + 1).cloned())
}

fn demo(modality: &str, ticks: u64) {
    let (size, frames) = match modality {
        "audio" => (4096, 0),
        "video" => (48, 12),
        _ => (64, 0),
    };
    let cfg = Config {
        experiment: ExperimentCfg {
            name: format!("demo-{modality}"),
            seeds: vec![42],
            max_ticks: ticks,
            log_every: 50,
            snapshot_every: ticks / 4,
            mp_candidates: 24,
        },
        world: WorldCfg {
            modality: modality.into(),
            target: "synthetic".into(),
            size,
            frames,
            sample_rate: 16000,
            target_b: None,
            switch_at: 0,
        },
        life: LifeConfig::default(),
        atoms: AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 },
    };
    dispatch(&cfg, "# demo builtin", 42, Path::new("."));
}

/// Dispatch de modalidade: escolhe átomo + IO; a física é a mesma.
fn dispatch(cfg: &Config, toml_src: &str, seed: u64, out_root: &Path) {
    match cfg.world.modality.as_str() {
        "image" => {
            let target_dc = load_image_target(&cfg.world);
            run_generic::<Gabor2D>(cfg, toml_src, seed, out_root, target_dc, Viz::Image);
        }
        "audio" => {
            let target_dc = load_audio_target(&cfg.world);
            run_generic::<Gabor1D>(
                cfg,
                toml_src,
                seed,
                out_root,
                target_dc,
                Viz::Audio { sample_rate: cfg.world.sample_rate },
            );
        }
        "video" => {
            let target_dc = load_video_target(&cfg.world);
            run_generic::<Gabor3D>(cfg, toml_src, seed, out_root, target_dc, Viz::Video);
        }
        m => panic!("modalidade desconhecida: {m} (use image|audio|video)"),
    }
}

// ─────────────────────────── carregamento de alvos ───────────────────────────

fn load_image_target(wcfg: &WorldCfg) -> Field {
    if wcfg.target == "synthetic" {
        return synthetic_image(wcfg.size);
    }
    let img = image::open(&wcfg.target)
        .unwrap_or_else(|e| panic!("abrindo alvo {}: {e}", wcfg.target))
        .resize_exact(
            wcfg.size as u32,
            wcfg.size as u32,
            image::imageops::FilterType::CatmullRom,
        )
        .to_luma8();
    Field {
        dims: vec![wcfg.size, wcfg.size],
        data: img.as_raw().iter().map(|v| *v as f32 / 255.0).collect(),
    }
}

fn load_audio_target(wcfg: &WorldCfg) -> Field {
    if wcfg.target == "synthetic" {
        return synthetic_audio(wcfg.size);
    }
    let (samples, _sr) =
        load_wav(Path::new(&wcfg.target)).unwrap_or_else(|e| panic!("abrindo {}: {e}", wcfg.target));
    let n = wcfg.size.min(samples.len());
    Field { dims: vec![n], data: samples[..n].to_vec() }
}

fn load_video_target(wcfg: &WorldCfg) -> Field {
    if wcfg.target == "synthetic" {
        return synthetic_video(wcfg.frames, wcfg.size);
    }
    panic!("vídeo v0 só suporta target = \"synthetic\" (carregar frames PNG entra na bateria de vídeo)");
}

// ─────────────────────────── visualização por modalidade ───────────────────────────

enum Viz {
    Image,
    Audio { sample_rate: u32 },
    Video,
}

impl Viz {
    /// Salva a visualização de um campo (com DC) nos snapshots do run.
    fn save(&self, rd: &RunDir, name: &str, f: &Field, signed: bool) {
        match self {
            Viz::Image => {
                if signed {
                    rd.save_signed_png(&format!("{name}.png"), f, 0.5).unwrap();
                } else {
                    rd.save_field_png(&format!("{name}.png"), f).unwrap();
                }
            }
            Viz::Audio { .. } => {
                waveform_png(&rd.snapshot_path(&format!("{name}.png")), &f.data, 512, 128)
                    .unwrap();
            }
            Viz::Video => {
                // primeiro / meio / último frame
                let t = f.dims[0];
                for (tag, ti) in [("f0", 0usize), ("fm", t / 2), ("fl", t - 1)] {
                    let frame = video_frame(f, ti);
                    if signed {
                        rd.save_signed_png(&format!("{name}_{tag}.png"), &frame, 0.5).unwrap();
                    } else {
                        rd.save_field_png(&format!("{name}_{tag}.png"), &frame).unwrap();
                    }
                }
            }
        }
    }

    /// Export final extra (áudio: WAVs ouvíveis do alvo e da reconstrução).
    fn export_final(&self, rd: &RunDir, target_dc: &Field, recon_dc: &Field) {
        if let Viz::Audio { sample_rate } = self {
            save_wav(&rd.snapshot_path("target.wav"), &target_dc.data, *sample_rate).unwrap();
            save_wav(&rd.snapshot_path("recon.wav"), &recon_dc.data, *sample_rate).unwrap();
        }
    }
}

fn video_frame(f: &Field, t: usize) -> Field {
    let (h, w) = (f.dims[1], f.dims[2]);
    let start = t * h * w;
    Field { dims: vec![h, w], data: f.data[start..start + h * w].to_vec() }
}

// ─────────────────────────── o run genérico (física única) ───────────────────────────

fn rd_json(p: &RdPoint) -> serde_json::Value {
    serde_json::json!({"label": p.label, "atoms": p.atoms, "bytes": p.bytes, "psnr_db": p.psnr_db})
}

fn run_generic<A: Atom>(
    cfg: &Config,
    toml_src: &str,
    seed: u64,
    out_root: &Path,
    target_dc: Field,
    viz: Viz,
) {
    let exp = &cfg.experiment;
    let t_start = std::time::Instant::now();

    // alvo com DC (pra viz/JPEG) e versão zero-média (pro mundo)
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
        "world": {
            "modality": cfg.world.modality,
            "target": cfg.world.target,
            "dims": target.dims,
            "dc": dc,
        },
        "max_ticks": exp.max_ticks,
        "mp_candidates": exp.mp_candidates,
    }))
    .unwrap();
    rd.copy_config(toml_src).unwrap();
    viz.save(&rd, "target", &target_dc, false);

    // alvo B (adaptação): mesmo DC-shift do alvo A pra comparar em pé de igualdade
    let target_b: Option<Field> = cfg.world.target_b.as_ref().map(|p| {
        let mut wb = cfg.world.clone();
        wb.target = p.clone();
        let mut tb = match cfg.world.modality.as_str() {
            "image" => load_image_target(&wb),
            "audio" => load_audio_target(&wb),
            _ => panic!("target_b só image/audio"),
        };
        let dcb = tb.mean();
        for v in tb.data.iter_mut() {
            *v -= dcb;
        }
        tb
    });

    let mut world: World<A> = World::new(target.clone(), cfg.life.clone(), cfg.atoms.clone(), seed);
    let mut best_psnr = world.psnr(1.0);
    let mut evals_at_switch = 0u64;
    let mut psnr_on_switch = 0.0f32; // PSNR em B ANTES de adaptar (população de A crua)

    println!(
        "▶ {} seed={} [{}] | dims {:?} | pop inicial {} | psnr0 {:.2} dB",
        exp.name,
        seed,
        cfg.world.modality,
        target.dims,
        world.beings.len(),
        best_psnr
    );

    for tick in 1..=exp.max_ticks {
        // TROCA DE SINAL: o relevo muda; a água tem que reescoar
        if let (Some(tb), true) = (target_b.as_ref(), cfg.world.switch_at == tick) {
            evals_at_switch = world.evals;
            world.target = tb.clone();
            world.rebuild();
            psnr_on_switch = world.psnr(1.0);
            println!("  ↻ tick {tick}: sinal trocou; psnr sem adaptar = {psnr_on_switch:.2} dB");
        }
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
        if exp.snapshot_every > 0 && tick % exp.snapshot_every == 0 {
            let recon_dc = with_dc(&world.recon, dc);
            viz.save(&rd, &format!("recon_{tick:06}"), &recon_dc, false);
            viz.save(&rd, &format!("residual_{tick:06}"), &world.residual, true);
        }
    }

    world.rebuild(); // números finais sem deriva float
    let final_psnr = world.psnr(1.0);
    let alive = world.beings.len();
    let bytes = world.bytes();
    let wall = t_start.elapsed().as_secs_f32();

    // ---- baselines honestos com o MESMO número de átomos e MESMA família ----
    let n = alive.max(1);
    let cps: Vec<usize> = [n / 4, n / 2, n].iter().copied().filter(|v| *v > 0).collect();
    let rnd = random_atoms::<A>(&target, n, &cfg.atoms, seed ^ 0xABCD, &cps);
    // iso-budget: MP recebe o MESMO orçamento de avaliações que a sobrevivência gastou.
    // No modo adaptação, o orçamento justo é o da FASE B (resposta à mudança) —
    // o MP é offline e teria que resolver B do zero com esse orçamento.
    let budget = if target_b.is_some() {
        world.evals - evals_at_switch
    } else {
        world.evals
    };
    let mp_cands = if exp.mp_candidates == 0 {
        ((budget as usize) / n).max(1)
    } else {
        exp.mp_candidates
    };
    let mp = matching_pursuit::<A>(&target, n, mp_cands, &cfg.atoms, seed ^ 0x1234, &cps);
    let jpeg = if cfg.world.modality == "image" {
        jpeg_curve(&target_dc, &[10, 30, 50, 70, 90])
    } else {
        Vec::new()
    };
    let rnd_final = rnd.last().map(|p| p.psnr_db).unwrap_or(0.0);
    let mp_final = mp.last().map(|p| p.psnr_db).unwrap_or(0.0);

    rd.write_population(&serde_json::to_value(world.population_records()).unwrap())
        .unwrap();
    rd.write_summary(&serde_json::json!({
        "modality": cfg.world.modality,
        "final_psnr_db": final_psnr,
        "best_psnr_db": best_psnr,
        "alive": alive,
        "bytes": bytes,
        "births": world.births,
        "deaths": world.deaths,
        "evals": world.evals,
        "adapt": target_b.is_some(),
        "switch_at": cfg.world.switch_at,
        "evals_phase_b": world.evals - evals_at_switch,
        "psnr_on_switch_db": psnr_on_switch,
        "mp_candidates": mp_cands,
        "mp_evals": n * mp_cands,
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

    // visualização final
    let recon_dc = with_dc(&world.recon, dc);
    viz.save(&rd, "recon_final", &recon_dc, false);
    viz.save(&rd, "residual_final", &world.residual, true);
    viz.export_final(&rd, &target_dc, &recon_dc);

    let jpeg_str = if jpeg.is_empty() {
        "n/a".to_string()
    } else {
        jpeg.iter()
            .map(|p| format!("{}B={:.1}dB", p.bytes, p.psnr_db))
            .collect::<Vec<_>>()
            .join(" ")
    };
    println!(
        "✔ {} seed={} | {:.1}s | psnr {:.2} dB ({} seres, {} B) | random {:.2} | MP {:.2} | jpeg: {}",
        exp.name, seed, wall, final_psnr, alive, bytes, rnd_final, mp_final, jpeg_str
    );

    let verdict = if final_psnr > rnd_final + 0.5 {
        format!(
            "bate random (+{:.2} dB), gap MP {:.2} dB",
            final_psnr - rnd_final,
            mp_final - final_psnr
        )
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

fn with_dc(f: &Field, dc: f32) -> Field {
    let mut out = f.clone();
    for v in out.data.iter_mut() {
        *v += dc;
    }
    out
}
