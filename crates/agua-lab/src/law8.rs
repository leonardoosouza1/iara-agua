//! `agua-lab law8` — O N CRÍTICO DA PREDAÇÃO.
//!
//! Lei 7 (iso-N): predação ganha +0.29 dB com N=120 (3/3) e PERDE com N=250
//! (−1.88, 0/3) e N=500 (−1.59, 0/3). Existe um N* onde o sinal inverte.
//!
//! HIPÓTESE: a predação paga enquanto o resíduo é dominado por SINAL ainda não
//! explicado (aí redundância é desperdício). Quando N cresce o bastante, o que
//! falta é detalhe fino e a redundância vira CAPACIDADE (um frame redundante
//! representa melhor que uma base ortogonal do mesmo tamanho).
//!
//! PREDIÇÃO: N* é MAIOR num sinal mais complexo — ele precisa de mais átomos
//! antes que a capacidade vire o gargalo. Se N* for igual nos dois, a hipótese
//! da "complexidade" morre e sobra alguma constante geométrica do átomo.
//!
//! Medimos a complexidade independentemente: nº de átomos que o matching pursuit
//! precisa pra atingir 30 dB (proxy honesto e externo à população).

use agua_atoms::{AtomBounds, Gabor2D};
use agua_core::{psnr, Field, Rng};
use agua_life::{LifeConfig, World};
use agua_metrics::matching_pursuit;

/// dois sinais, complexidade deliberadamente diferente
fn alvo(nome: &str, size: usize, seed: u64) -> Field {
    let s = size as f32;
    let mut rng = Rng::new(seed);
    let jx = rng.range(0.28, 0.42);
    let mut f = match nome {
        // SIMPLES: um blob suave + uma faixa de listras
        "simples" => Field::from_fn_2d(size, size, |x, y| {
            let (xf, yf) = (x as f32, y as f32);
            let blob = 0.5
                * (-((xf - s * jx).powi(2) + (yf - s * 0.35).powi(2)) / (2.0 * (s * 0.16).powi(2)))
                    .exp();
            let stripes = if xf > s * 0.6 { 0.22 * (std::f32::consts::TAU * 0.12 * xf).sin() } else { 0.0 };
            blob + stripes
        }),
        // COMPLEXO: 6 blobs + 3 bandas de frequência + textura fina em todo lugar
        _ => Field::from_fn_2d(size, size, |x, y| {
            let (xf, yf) = (x as f32, y as f32);
            let mut v = 0.0;
            for k in 0..6 {
                let a = k as f32 * 1.7 + jx;
                let (bx, by) = (s * (0.2 + 0.6 * (a.sin() * 0.5 + 0.5)), s * (0.2 + 0.6 * (a.cos() * 0.5 + 0.5)));
                v += 0.35 * (-((xf - bx).powi(2) + (yf - by).powi(2)) / (2.0 * (s * 0.07).powi(2))).exp();
            }
            v += 0.15 * (std::f32::consts::TAU * 0.10 * xf).sin();
            v += 0.12 * (std::f32::consts::TAU * 0.22 * yf).sin();
            v += 0.10 * (std::f32::consts::TAU * 0.35 * (xf + yf)).sin();
            v
        }),
    };
    let m = f.mean();
    for v in f.data.iter_mut() {
        *v -= m;
    }
    f
}

fn cfg(met: f32, pred: f32) -> LifeConfig {
    let mut c = LifeConfig::default();
    c.cap = 2000;
    c.init = 64;
    c.metabolism = met;
    c.rain_every = 1;
    c.rain_count = 6;
    c.rain_competition = 16;
    c.tune_every = 1;
    c.gradient = true;
    c.predation = pred;
    c.predation_every = 5;
    c
}

fn run_one(t: &Field, b: &AtomBounds, met: f32, pred: f32, seed: u64, ticks: u64) -> (usize, f32) {
    let mut w: World<Gabor2D> = World::new(t.clone(), cfg(met, pred), b.clone(), seed);
    for _ in 0..ticks {
        w.step();
    }
    w.rebuild();
    (w.beings.len(), w.psnr(1.0))
}

fn met_for_n(t: &Field, b: &AtomBounds, n_alvo: usize, pred: f32, seed: u64, ticks: u64) -> f32 {
    let (mut lo, mut hi) = (0.01f32, 8.0f32);
    let mut best = 0.0f32;
    for _ in 0..8 {
        let m = (lo * hi).sqrt();
        let (n, p) = run_one(t, b, m, pred, seed, ticks);
        best = p;
        if n > n_alvo {
            lo = m;
        } else {
            hi = m;
        }
        if (n as f32 - n_alvo as f32).abs() / (n_alvo as f32) < 0.06 {
            break;
        }
    }
    best
}

/// complexidade = nº de átomos que o MP (externo à população) gasta pra 30 dB
fn complexidade(t: &Field, b: &AtomBounds, seed: u64) -> usize {
    let cps: Vec<usize> = (1..=12).map(|k| k * 50).collect();
    let pts = matching_pursuit::<Gabor2D>(t, 600, 24, b, seed, &cps);
    pts.iter().find(|p| p.psnr_db >= 30.0).map(|p| p.atoms).unwrap_or(600)
}

pub fn run(seed_base: u64) {
    let size = 64usize;
    let b = AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
    let ticks = 700u64;

    println!("=== O N CRÍTICO DA PREDAÇÃO ===");
    println!("Predação paga no esparso e prejudica no denso. Onde está N*?");
    println!("PREDIÇÃO: N* maior no sinal MAIS COMPLEXO.\n");

    for nome in ["simples", "complexo"] {
        let t0 = alvo(nome, size, 11);
        let cx = complexidade(&t0, &b, 7);
        println!("--- sinal {nome} (complexidade externa: {cx} átomos de MP pra 30 dB) ---");
        println!("{:>8} {:>14} {:>14} {:>10} {:>10}", "N", "sem predação", "com predação", "Δ (dB)", "vitórias");
        for &n in &[60usize, 120, 180, 250, 400] {
            let mut sem = Vec::new();
            let mut com = Vec::new();
            for seed in [11u64, 22, 33] {
                let t = alvo(nome, size, seed);
                sem.push(met_for_n(&t, &b, n, 0.0, seed ^ seed_base, ticks));
                com.push(met_for_n(&t, &b, n, 0.85, seed ^ seed_base, ticks));
            }
            let ms = sem.iter().sum::<f32>() / 3.0;
            let mc = com.iter().sum::<f32>() / 3.0;
            let w = sem.iter().zip(&com).filter(|(a, b)| b > a).count();
            let marca = if mc > ms { "  ⟵ predação ganha" } else { "" };
            println!("{n:>8} {ms:>14.2} {mc:>14.2} {:>+10.2} {w:>7}/3{marca}", mc - ms);
        }
        println!();
    }
    println!("Se N* cresce com a complexidade, o gargalo é CAPACIDADE.");
    println!("Se N* é o mesmo, é uma constante geométrica do átomo (empacotamento).");
    let _ = psnr;
}
