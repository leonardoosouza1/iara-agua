//! `agua-lab probe` — medições primitivas pra baterias: 1 invocação = 1 JSON.
//!
//! Kinds:
//!   c           — fator de busca c = E[explicado]/σ² em ruído puro
//!   morozov     — denoise com met = k·σ̂²; PSNR vs limpo + discrepância
//!   equilibrium — população em sinal limpo: testa N·met ≈ E_capturada (Lei da População)
//!   spectrum    — histograma de freq da população vs energia por banda (Lei do Espelho)
//!   gradcheck   — N átomos aleatórios: gradiente analítico vs diferenças finitas
//!   drift       — run longo: deriva do invariante + estabilidade numérica
//!   adapt       — troca de sinal no meio: velocidade de re-harmonização
//!   perf        — wall time, evals/s, RSS de pico (VmHWM), por tamanho
//!   codec       — população real → .agua nos 4 perfis → decode → PSNR (armazenamento)

use crate::denoise::estimate_sigma;
use crate::targets::{synthetic_audio, synthetic_image};
use agua_atoms::{Atom, AtomBounds, Gabor1D, Gabor2D};
use agua_core::{psnr, Field, Rng};
use agua_life::{LifeConfig, World};
use serde_json::json;
use std::path::Path;

fn arg(args: &[String], k: &str) -> Option<String> {
    args.iter().position(|a| a == k).and_then(|i| args.get(i + 1).cloned())
}
fn argf(args: &[String], k: &str, d: f32) -> f32 {
    arg(args, k).map_or(d, |v| v.parse().unwrap())
}
fn argu(args: &[String], k: &str, d: u64) -> u64 {
    arg(args, k).map_or(d, |v| v.parse().unwrap())
}
fn args_(args: &[String], k: &str, d: &str) -> String {
    arg(args, k).unwrap_or_else(|| d.into())
}

fn bounds2() -> AtomBounds {
    AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 }
}
fn bounds1() -> AtomBounds {
    AtomBounds { freq_min: 0.002, freq_max: 0.45, sigma_min: 2.0, sigma_max: 400.0 }
}

/// ruído com média 0 e variância σ² — formas diferentes (lei 5)
fn ruido(dims: &[usize], sigma: f32, seed: u64, tipo: &str) -> Field {
    let mut rng = Rng::new(seed);
    let n: usize = dims.iter().product();
    let data: Vec<f32> = (0..n)
        .map(|_| match tipo {
            "laplace" => {
                let u = rng.f32() - 0.5;
                let b = sigma / std::f32::consts::SQRT_2;
                -b * u.signum() * (1.0 - 2.0 * u.abs()).max(1e-9).ln()
            }
            "uniforme" => {
                let a = sigma * 3.0f32.sqrt();
                rng.range(-a, a)
            }
            t if t.starts_with("sal") => {
                let p: f32 = t[4..].parse().unwrap_or(0.01);
                if rng.f32() < p {
                    let a = sigma / p.sqrt();
                    if rng.f32() < 0.5 { a } else { -a }
                } else {
                    0.0
                }
            }
            _ => rng.normal() * sigma,
        })
        .collect();
    Field { dims: dims.to_vec(), data }
}

fn add_noise(f: &Field, sigma: f32, seed: u64) -> Field {
    let mut rng = Rng::new(seed);
    let mut o = f.clone();
    for v in o.data.iter_mut() {
        *v += rng.normal() * sigma;
    }
    o
}

fn zero_mean(f: &Field) -> (Field, f32) {
    let dc = f.mean();
    let mut o = f.clone();
    for v in o.data.iter_mut() {
        *v -= dc;
    }
    (o, dc)
}

fn mse(a: &Field, b: &Field) -> f32 {
    a.data.iter().zip(&b.data).map(|(x, y)| (x - y) * (x - y)).sum::<f32>() / a.len() as f32
}

/// carrega alvo: "synthetic" | "synthetic:<seed>" | caminho .png | caminho .wav
fn load_target(modality: &str, target: &str, size: usize) -> Field {
    if target.starts_with("synthetic") {
        return match modality {
            "audio" => synthetic_audio(size),
            _ => synthetic_image(size),
        };
    }
    if modality == "audio" {
        let (s, _) = agua_telemetry::load_wav(Path::new(target)).expect("wav");
        let n = size.min(s.len());
        Field { dims: vec![n], data: s[..n].to_vec() }
    } else {
        let img = image::open(target)
            .expect("img")
            .resize_exact(size as u32, size as u32, image::imageops::FilterType::CatmullRom)
            .to_luma8();
        Field { dims: vec![size, size], data: img.as_raw().iter().map(|v| *v as f32 / 255.0).collect() }
    }
}

fn cfg_c(comp: usize, grad: bool, pred: f32, rain: bool) -> LifeConfig {
    let mut c = LifeConfig::default();
    c.cap = 400;
    c.init = if rain { 200 } else { 400 };
    c.rain_every = if rain { 1 } else { 0 };
    c.rain_count = 4;
    c.rain_competition = comp;
    c.tune_every = 1;
    c.gradient = grad;
    c.metabolism = 1e-6;
    c.initial_energy = 50.0;
    c.repro_threshold = 1e9;
    c.predation = pred;
    c.predation_every = 5;
    c
}

fn cfg_denoise(met: f32) -> LifeConfig {
    let mut c = LifeConfig::default();
    c.cap = 3000;
    c.metabolism = met;
    c.rain_every = 1;
    c.rain_count = 6;
    c.rain_competition = 16;
    c.tune_every = 1;
    c.gradient = true;
    c
}

fn run_world<A: Atom>(target: Field, cfg: LifeConfig, b: &AtomBounds, seed: u64, ticks: u64) -> World<A> {
    let mut w: World<A> = World::new(target, cfg, b.clone(), seed);
    for _ in 0..ticks {
        w.step();
    }
    w
}

fn kind_c(args: &[String]) {
    let modality = args_(args, "--modality", "image");
    let size = argu(args, "--size", 64) as usize;
    let sigma = argf(args, "--sigma", 0.1);
    let tipo = args_(args, "--noise", "gauss");
    let ticks = argu(args, "--ticks", 800);
    let comp = argu(args, "--comp", 8) as usize;
    let frozen: u32 = argu(args, "--frozen", 0) as u32;
    let pred = argf(args, "--pred", 0.0);
    let seed = argu(args, "--seed", 11);
    let rain = argu(args, "--rain", 1) == 1;

    let dims: Vec<usize> = if modality == "audio" { vec![size] } else { vec![size, size] };
    let campo = ruido(&dims, sigma, seed, &tipo);
    let var = campo.data.iter().map(|v| v * v).sum::<f32>() / campo.len() as f32;
    let mut cfg = cfg_c(comp, true, pred, rain);
    cfg.frozen = frozen;

    let (c, alive) = if modality == "audio" {
        let w = run_world::<Gabor1D>(campo, cfg, &bounds1(), seed, ticks);
        let e = w.explained_now();
        ((e.iter().sum::<f32>() / e.len().max(1) as f32) / var, w.beings.len())
    } else {
        let w = run_world::<Gabor2D>(campo, cfg, &bounds2(), seed, ticks);
        let e = w.explained_now();
        ((e.iter().sum::<f32>() / e.len().max(1) as f32) / var, w.beings.len())
    };
    println!(
        "{}",
        json!({"kind":"c","modality":modality,"size":size,"sigma":sigma,"noise":tipo,
               "ticks":ticks,"comp":comp,"frozen":frozen,"pred":pred,"seed":seed,
               "c":c,"alive":alive})
    );
}

fn kind_morozov(args: &[String]) {
    let modality = args_(args, "--modality", "image");
    let size = argu(args, "--size", 96) as usize;
    let sigma = argf(args, "--sigma", 0.15);
    let k = argf(args, "--k", 40.0);
    let ticks = argu(args, "--ticks", 800);
    let seed = argu(args, "--seed", 11);
    let target = args_(args, "--target", "synthetic");
    let oracle = argu(args, "--oracle", 0) == 1;

    let clean = load_target(&modality, &target, size);
    let noisy = add_noise(&clean, sigma, seed ^ 0x0157);
    let sig_hat = if oracle { sigma } else { estimate_sigma(&noisy) };
    let (noisy0, dc) = zero_mean(&noisy);
    let met = k * sig_hat * sig_hat;

    let (p, disc, alive) = if modality == "audio" {
        let mut w = run_world::<Gabor1D>(noisy0, cfg_denoise(met), &bounds1(), seed, ticks);
        w.rebuild();
        let mut rec = w.recon.clone();
        for v in rec.data.iter_mut() {
            *v += dc;
        }
        (psnr(mse(&rec, &clean), 1.0), w.residual.energy() / w.residual.len() as f32, w.beings.len())
    } else {
        let mut w = run_world::<Gabor2D>(noisy0, cfg_denoise(met), &bounds2(), seed, ticks);
        w.rebuild();
        let mut rec = w.recon.clone();
        for v in rec.data.iter_mut() {
            *v += dc;
        }
        (psnr(mse(&rec, &clean), 1.0), w.residual.energy() / w.residual.len() as f32, w.beings.len())
    };
    println!(
        "{}",
        json!({"kind":"morozov","modality":modality,"target":target,"size":size,"sigma":sigma,
               "k":k,"oracle":oracle,"ticks":ticks,"seed":seed,"sig_hat":sig_hat,
               "psnr":p,"disc_ratio":disc/(sigma*sigma),"alive":alive,
               "psnr_noisy":psnr(mse(&noisy,&clean),1.0)})
    );
}

fn kind_equilibrium(args: &[String]) {
    let size = argu(args, "--size", 64) as usize;
    let met = argf(args, "--met", 0.2);
    let ticks = argu(args, "--ticks", 1500);
    let seed = argu(args, "--seed", 11);
    let target = args_(args, "--target", "synthetic");

    let clean = load_target("image", &target, size);
    let (t0, _) = zero_mean(&clean);
    let e_signal = t0.energy();
    let mut w = run_world::<Gabor2D>(t0, cfg_denoise(met), &bounds2(), seed, ticks);
    w.rebuild();
    let e = w.explained_now();
    let total_expl: f32 = e.iter().sum();
    let e_captured = w.recon.energy();
    println!(
        "{}",
        json!({"kind":"equilibrium","target":target,"size":size,"met":met,"ticks":ticks,"seed":seed,
               "alive":w.beings.len(),"explained_mean":total_expl/e.len().max(1) as f32,
               "explained_total":total_expl,"e_signal":e_signal,"e_captured":e_captured,
               "psnr":w.psnr(1.0),
               // Lei da População: N ≈ E_capturada/met (da equipartição)
               "n_pred":e_captured/met})
    );
}

fn kind_spectrum(args: &[String]) {
    let size = argu(args, "--size", 64) as usize;
    let met = argf(args, "--met", 0.2);
    let ticks = argu(args, "--ticks", 1500);
    let seed = argu(args, "--seed", 11);
    let target = args_(args, "--target", "synthetic");

    let clean = load_target("image", &target, size);
    let (t0, _) = zero_mean(&clean);
    let mut w = run_world::<Gabor2D>(t0, cfg_denoise(met), &bounds2(), seed, ticks);
    w.rebuild();
    let edges = [0.03f32, 0.06, 0.12, 0.24];
    let pop = w.freq_bands(&edges);
    // energia REALMENTE capturada por banda (recon_band) — o "espectro do sinal
    // visto pela população"
    let mut band_energy = Vec::new();
    let mut lo = 0.0f32;
    for e in edges.iter().chain(std::iter::once(&10.0f32)) {
        band_energy.push(w.recon_band(lo, *e).energy());
        lo = *e;
    }
    println!(
        "{}",
        json!({"kind":"spectrum","target":target,"size":size,"met":met,"seed":seed,
               "band_pop":pop,"band_energy":band_energy,"alive":w.beings.len(),"psnr":w.psnr(1.0)})
    );
}

/// trava 2.5σ em x.5 (senão a diferença finita mede o salto da janela de
/// suporte, não o gradiente — artefato que já derrubou testes válidos)
fn snap<A: Atom + 'static>(g: A, bd: &AtomBounds) -> A {
    use std::any::Any;
    // 1) σ longe dos LIMITES (senão apply_delta clampa um lado da diferença
    //    central e a medida vira one-sided); 2) 2.5σ termina em .5 (janela de
    //    suporte não salta com ±eps)
    let fix = |sig: f32| -> f32 {
        let s = sig.clamp(bd.sigma_min + 0.4, bd.sigma_max - 0.4);
        ((2.5 * s).floor() + 0.5) / 2.5
    };
    let any: &dyn Any = &g;
    if let Some(g2) = any.downcast_ref::<Gabor2D>() {
        let mut c = g2.clone();
        c.sigma = fix(c.sigma);
        let b: &dyn Any = &c;
        return b.downcast_ref::<A>().unwrap().clone();
    }
    if let Some(g1) = any.downcast_ref::<Gabor1D>() {
        let mut c = g1.clone();
        c.sigma = fix(c.sigma);
        let b: &dyn Any = &c;
        return b.downcast_ref::<A>().unwrap().clone();
    }
    g
}

fn kind_gradcheck(args: &[String]) {
    let n = argu(args, "--n", 100) as usize;
    let seed = argu(args, "--seed", 11);
    let modality = args_(args, "--modality", "image");

    fn rel_energy(unit: &agua_core::Patch, rw: &Field, amp: f32) -> f64 {
        -(amp as f64) * unit.dot_field(rw) as f64 + 0.5 * (amp as f64).powi(2) * unit.norm2() as f64
    }
    fn check<A: Atom + 'static>(dims: &[usize], b: &AtomBounds, n: usize, seed: u64, eps: &[f32]) -> (usize, usize, f32) {
        let mut rng = Rng::new(seed);
        let tot: usize = dims.iter().product();
        let (mut pass, mut fail, mut worst) = (0usize, 0usize, 0.0f32);
        for _ in 0..n {
            let rw = Field { dims: dims.to_vec(), data: (0..tot).map(|_| rng.normal() * 0.3).collect() };
            let pos: Vec<f32> = dims.iter().map(|d| rng.range(0.25 * *d as f32, 0.75 * *d as f32) + 0.37).collect();
            let g = A::rain_at(&pos, b, &mut rng);
            let g = snap::<A>(g, b);
            let amp = rng.range(-1.0, 1.0);
            let mut resid = rw.clone();
            resid.add_scaled_patch(&g.render_unit(dims), -amp);
            let ana = g.grad(dims, &resid, amp);
            // critério VETORIAL: o que importa pro Adam é a DIREÇÃO do gradiente.
            // Componentes minúsculas somem no ruído f32 do render (ΔE ~1e-6 sobre
            // energias ~1) — compará-las uma a uma mede ruído, não erro.
            let mut num = vec![0.0f32; A::GEOM];
            for k in 0..A::GEOM {
                let mut dp = vec![0.0; A::GEOM];
                dp[k] = eps[k];
                let mut gp = g.clone();
                gp.apply_delta(&dp, b, dims);
                dp[k] = -eps[k];
                let mut gm = g.clone();
                gm.apply_delta(&dp, b, dims);
                num[k] = ((rel_energy(&gp.render_unit(dims), &rw, amp)
                    - rel_energy(&gm.render_unit(dims), &rw, amp))
                    / (2.0 * eps[k] as f64)) as f32;
            }
            let dot: f32 = ana.iter().zip(&num).map(|(a, n)| a * n).sum();
            let na: f32 = ana.iter().map(|v| v * v).sum::<f32>().sqrt();
            let nn: f32 = num.iter().map(|v| v * v).sum::<f32>().sqrt();
            let cos = if na * nn > 1e-9 { dot / (na * nn) } else { 1.0 };
            let mag = if nn > 1e-6 { (na / nn - 1.0).abs() } else { 0.0 };
            worst = worst.max(1.0 - cos);
            if (cos > 0.95 && mag < 0.25) || na < 1e-3 {
                pass += 1
            } else {
                fail += 1
            }
        }
        (pass, fail, worst)
    }
    let (pass, fail, worst) = if modality == "audio" {
        check::<Gabor1D>(&[2048], &bounds1(), n, seed, &[2e-3, 2e-6, 1e-3, 2e-3])
    } else {
        check::<Gabor2D>(&[48, 48], &bounds2(), n, seed, &[2e-3, 2e-3, 2e-5, 1e-3, 1e-3, 2e-3])
    };
    println!(
        "{}",
        json!({"kind":"gradcheck","modality":modality,"n":n,"seed":seed,
               "pass":pass,"fail":fail,"worst_rel":worst})
    );
}

fn kind_drift(args: &[String]) {
    let size = argu(args, "--size", 64) as usize;
    let ticks = argu(args, "--ticks", 20000);
    let seed = argu(args, "--seed", 11);
    let met = argf(args, "--met", 0.2);
    let clean = load_target("image", "synthetic", size);
    let (t0, _) = zero_mean(&clean);
    let t = std::time::Instant::now();
    let w = run_world::<Gabor2D>(t0, cfg_denoise(met), &bounds2(), seed, ticks);
    println!(
        "{}",
        json!({"kind":"drift","size":size,"ticks":ticks,"seed":seed,"met":met,
               "invariant_error":w.invariant_error(),"psnr":w.psnr(1.0),
               "alive":w.beings.len(),"evals":w.evals,"wall_s":t.elapsed().as_secs_f32()})
    );
}

fn kind_adapt(args: &[String]) {
    let size = argu(args, "--size", 96) as usize;
    let ticks_a = argu(args, "--ticks-a", 2000);
    let ticks_b = argu(args, "--ticks-b", 600);
    let seed = argu(args, "--seed", 11);
    let met = argf(args, "--met", 0.5);
    let ta = args_(args, "--target", "assets/photo_texture.png");
    let tb = args_(args, "--target-b", "assets/photo_texture_b.png");

    let a = load_target("image", &ta, size);
    let b_ = load_target("image", &tb, size);
    let (a0, _) = zero_mean(&a);
    let (b0, _) = zero_mean(&b_);
    let mut w = run_world::<Gabor2D>(a0, cfg_denoise(met), &bounds2(), seed, ticks_a);
    w.rebuild();
    let psnr_a = w.psnr(1.0);
    // troca o relevo
    w.target = b0;
    w.rebuild();
    let psnr_on = w.psnr(1.0);
    let evals0 = w.evals;
    let mut curve = Vec::new();
    for i in 0..ticks_b {
        w.step();
        if (i + 1) % (ticks_b / 6).max(1) == 0 {
            curve.push((i + 1, w.psnr(1.0)));
        }
    }
    w.rebuild();
    println!(
        "{}",
        json!({"kind":"adapt","size":size,"met":met,"seed":seed,
               "psnr_a":psnr_a,"psnr_on_switch":psnr_on,"psnr_final":w.psnr(1.0),
               "recovery_db":w.psnr(1.0)-psnr_on,"evals_b":w.evals-evals0,
               "curve":curve.iter().map(|(t,p)| json!([t,p])).collect::<Vec<_>>()})
    );
}

fn vm_hwm_kb() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|v| v.parse().ok())
        })
        .unwrap_or(0)
}

fn kind_perf(args: &[String]) {
    let size = argu(args, "--size", 128) as usize;
    let ticks = argu(args, "--ticks", 1000);
    let seed = argu(args, "--seed", 11);
    let cap = argu(args, "--cap", 3000) as usize;
    let met = argf(args, "--met", 0.1);
    let clean = load_target("image", "synthetic", size);
    let (t0, _) = zero_mean(&clean);
    let mut cfg = cfg_denoise(met);
    cfg.cap = cap;
    let t = std::time::Instant::now();
    let w = run_world::<Gabor2D>(t0, cfg, &bounds2(), seed, ticks);
    let wall = t.elapsed().as_secs_f32();
    println!(
        "{}",
        json!({"kind":"perf","size":size,"ticks":ticks,"cap":cap,"met":met,"seed":seed,
               "wall_s":wall,"evals":w.evals,"evals_per_s":w.evals as f32/wall,
               "alive":w.beings.len(),"psnr":w.psnr(1.0),"rss_peak_mb":vm_hwm_kb() as f32/1024.0,
               "ticks_per_s":ticks as f32/wall})
    );
}

fn kind_codec(args: &[String]) {
    use agua_codec::{decode, encode, render, Decoded, Profile};
    let size = argu(args, "--size", 96) as usize;
    let ticks = argu(args, "--ticks", 1500);
    let seed = argu(args, "--seed", 11);
    let met = argf(args, "--met", 0.3);
    let target = args_(args, "--target", "synthetic");

    let clean = load_target("image", &target, size);
    let (t0, dc) = zero_mean(&clean);
    let mut w = run_world::<Gabor2D>(t0.clone(), cfg_denoise(met), &bounds2(), seed, ticks);
    w.rebuild();
    let atoms: Vec<(Gabor2D, f32)> = w.beings.iter().map(|b| (b.atom.clone(), b.amp)).collect();
    let base = w.psnr(1.0);
    let mut rows = Vec::new();
    for (nome, p) in [("f32", Profile::F32), ("fp16", Profile::Fp16), ("q12", Profile::Q12), ("q8", Profile::Q8)] {
        let enc = encode(&atoms, &t0.dims, dc, &bounds2(), p);
        let dec: Decoded<Gabor2D> = decode(&enc.bytes).unwrap();
        let rec_dc = render(&dec);
        let mut rec = rec_dc.clone();
        for v in rec.data.iter_mut() {
            *v -= dc;
        }
        let ps = psnr(mse(&rec, &t0), 1.0);
        rows.push(json!({"profile":nome,"bytes":enc.bytes.len(),"psnr":ps}));
    }
    println!(
        "{}",
        json!({"kind":"codec","target":target,"size":size,"met":met,"seed":seed,
               "alive":atoms.len(),"psnr_f32_world":base,"profiles":rows})
    );
}


fn kind_gamma_fdr(args: &[String]) {
    // γ = met_ótimo / met_c corresponde a um QUANTIL FIXO da distribuição de
    // energia explicada pelo ruído? Se sim, γ é um controle de falsa-sobrevivência.
    let modality = args_(args, "--modality", "image");
    let size = argu(args, "--size", 96) as usize;
    let sigma = argf(args, "--sigma", 0.15);
    let ticks = argu(args, "--ticks", 800);
    let seed = argu(args, "--seed", 11);
    let target = args_(args, "--target", "synthetic");

    let dims: Vec<usize> = if modality == "audio" { vec![size] } else { vec![size, size] };

    // 1) distribuição de energia explicada em RUÍDO PURO (normalizada por σ²)
    let noise = ruido(&dims, sigma, seed, "gauss");
    let var = noise.data.iter().map(|v| v*v).sum::<f32>() / noise.len() as f32;
    let mut cfg = cfg_c(8, true, 0.0, true);
    cfg.cap = 1000; cfg.init = 400;
    let dist: Vec<f32> = if modality == "audio" {
        let w = run_world::<Gabor1D>(noise, cfg, &bounds1(), seed, ticks);
        w.explained_now().iter().map(|e| e/var).collect()
    } else {
        let w = run_world::<Gabor2D>(noise, cfg, &bounds2(), seed, ticks);
        w.explained_now().iter().map(|e| e/var).collect()
    };
    let mut d = dist.clone();
    d.sort_by(|a,b| a.partial_cmp(b).unwrap());
    let met_c = if d.is_empty() {0.0} else { d.iter().sum::<f32>()/d.len() as f32 }; // média = c

    // 2) met ótimo para denoise no MESMO σ (varre k = met/σ²)
    let clean = load_target(&modality, &target, size);
    let noisy = add_noise(&clean, sigma, seed ^ 0x0157);
    let (noisy0, dc) = zero_mean(&noisy);
    let ks = [4.0f32, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 40.0, 55.0, 75.0, 100.0];
    let (mut best_k, mut best_p) = (0.0f32, f32::MIN);
    for &k in &ks {
        let met = k*sigma*sigma;
        let p = if modality == "audio" {
            let mut w = run_world::<Gabor1D>(noisy0.clone(), cfg_denoise(met), &bounds1(), seed, ticks);
            w.rebuild(); let mut r=w.recon.clone(); for v in r.data.iter_mut(){*v+=dc;}
            psnr(mse(&r,&clean),1.0)
        } else {
            let mut w = run_world::<Gabor2D>(noisy0.clone(), cfg_denoise(met), &bounds2(), seed, ticks);
            w.rebuild(); let mut r=w.recon.clone(); for v in r.data.iter_mut(){*v+=dc;}
            psnr(mse(&r,&clean),1.0)
        };
        if p>best_p { best_p=p; best_k=k; }
    }
    // 3) em que QUANTIL da distribuição de ruído cai met_ótimo? (fração morta)
    let met_opt = best_k; // já em unidades de σ² (= energia explicada / σ²)
    let killed = d.iter().filter(|&&e| e < met_opt).count() as f32 / d.len().max(1) as f32;
    let gamma = if met_c>0.0 { met_opt/met_c } else {0.0};
    println!("{}", json!({"kind":"gamma_fdr","modality":modality,"target":target,"sigma":sigma,
        "seed":seed,"met_c":met_c,"met_opt":met_opt,"gamma":gamma,
        "frac_ruido_morto":killed,"frac_ruido_vivo":1.0-killed,"best_psnr":best_p,"n_dist":d.len()}));
}

pub fn run(args: &[String]) {
    match args_(args, "--kind", "c").as_str() {
        "c" => kind_c(args),
        "morozov" => kind_morozov(args),
        "equilibrium" => kind_equilibrium(args),
        "spectrum" => kind_spectrum(args),
        "gradcheck" => kind_gradcheck(args),
        "drift" => kind_drift(args),
        "adapt" => kind_adapt(args),
        "perf" => kind_perf(args),
        "codec" => kind_codec(args),
        "gamma_fdr" => kind_gamma_fdr(args),
        k => {
            eprintln!("kind desconhecido: {k}");
            std::process::exit(2);
        }
    }
}
