//! `agua-lab law4` — A PREDAÇÃO PAGA? (teste ISO-N, o único justo)
//!
//! H4 mostrou: predação por nicho reduz N e sobe PSNR/átomo. Mas isso sozinho
//! não prova nada — menos átomos com menos PSNR pode ser só andar na curva
//! taxa-distorção. O teste honesto é ISO-N: force o MESMO número de átomos com
//! e sem predação (calibrando o metabolismo por bisseção) e compare o PSNR.
//!
//! Predição (exclusão competitiva de Gause): com N fixo, a população que expulsa
//! redundantes cobre mais nichos distintos ⇒ explica mais sinal ⇒ PSNR maior.
//!
//! Se der empate, a predação é só um atalho pro mesmo ponto da curva — útil
//! (converge mais rápido) mas não é ganho de representação.

use agua_atoms::{AtomBounds, Gabor2D};
use agua_core::{Field, Rng};
use agua_life::{LifeConfig, World};

fn synth(size: usize, seed: u64) -> Field {
    let s = size as f32;
    let mut rng = Rng::new(seed);
    let jx = rng.range(0.25, 0.45);
    let mut f = Field::from_fn_2d(size, size, |x, y| {
        let (xf, yf) = (x as f32, y as f32);
        let blob = 0.5
            * (-((xf - s * jx).powi(2) + (yf - s * 0.35).powi(2)) / (2.0 * (s * 0.13).powi(2))).exp();
        let stripes = if xf > s * 0.55 {
            0.25 * (std::f32::consts::TAU * 0.15 * xf).sin()
        } else {
            0.0
        };
        let fine = if yf > s * 0.7 {
            0.15 * (std::f32::consts::TAU * 0.3 * yf).sin()
        } else {
            0.0
        };
        blob + stripes + fine
    });
    let m = f.mean();
    for v in f.data.iter_mut() {
        *v -= m;
    }
    f
}

fn base(met: f32, pred: f32) -> LifeConfig {
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

/// roda e devolve (N, psnr)
fn run_one(target: &Field, b: &AtomBounds, met: f32, pred: f32, seed: u64, ticks: u64) -> (usize, f32) {
    let mut w: World<Gabor2D> = World::new(target.clone(), base(met, pred), b.clone(), seed);
    for _ in 0..ticks {
        w.step();
    }
    w.rebuild();
    (w.beings.len(), w.psnr(1.0))
}

/// bisseção no metabolismo (log) até a população bater n_alvo (±5%)
fn met_for_n(target: &Field, b: &AtomBounds, n_alvo: usize, pred: f32, seed: u64, ticks: u64) -> (f32, usize, f32) {
    let (mut lo, mut hi) = (0.01f32, 8.0f32); // lo = muitos seres, hi = poucos
    let mut best = (0.0f32, 0usize, 0.0f32);
    for _ in 0..10 {
        let m = (lo * hi).sqrt();
        let (n, p) = run_one(target, b, m, pred, seed, ticks);
        best = (m, n, p);
        if n > n_alvo {
            lo = m; // população grande demais ⇒ subir o custo
        } else {
            hi = m;
        }
        if (n as f32 - n_alvo as f32).abs() / (n_alvo as f32) < 0.05 {
            break;
        }
    }
    best
}

pub fn run(seed_base: u64) {
    let size = 64usize;
    let b = AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
    let ticks = 900u64;

    println!("=== PREDAÇÃO POR NICHO, TESTE ISO-N ===");
    println!("O metabolismo é calibrado por bisseção pra que ambos cheguem ao MESMO N.");
    println!("Se a predação vence a N igual, ela melhora a REPRESENTAÇÃO (não só a taxa).\n");

    for &n_alvo in &[120usize, 250, 500] {
        println!("--- N alvo = {n_alvo} ---");
        println!("{:>8} {:>8} {:>10} {:>10} {:>10}", "seed", "θ", "N obtido", "met", "PSNR");
        let mut sem = Vec::new();
        let mut com = Vec::new();
        for seed in [11u64, 22, 33] {
            let target = synth(size, seed);
            let (m0, n0, p0) = met_for_n(&target, &b, n_alvo, 0.0, seed ^ seed_base, ticks);
            let (m1, n1, p1) = met_for_n(&target, &b, n_alvo, 0.85, seed ^ seed_base, ticks);
            println!("{seed:>8} {:>8} {n0:>10} {m0:>10.4} {p0:>10.2}", "off");
            println!("{seed:>8} {:>8.2} {n1:>10} {m1:>10.4} {p1:>10.2}", 0.85);
            sem.push(p0);
            com.push(p1);
        }
        let ms: f32 = sem.iter().sum::<f32>() / sem.len() as f32;
        let mc: f32 = com.iter().sum::<f32>() / com.len() as f32;
        let wins = sem.iter().zip(&com).filter(|(a, b)| b > a).count();
        println!(
            "  média: sem predação {ms:.2} dB | com predação {mc:.2} dB | Δ = {:+.2} dB | vitórias {wins}/3\n",
            mc - ms
        );
    }
}
