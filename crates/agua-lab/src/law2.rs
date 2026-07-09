//! `agua-lab law2` — A LEI DO METABOLISMO CRÍTICO, versão corrigida.
//!
//! ERRO DA v1 (law.rs): medi "população viva após T ticks". Mas o tempo até a
//! morte é energia_inicial/(met − explicado) — com met pequeno o ser demora a
//! morrer e PARECE sobreviver. Aquilo media tempo de morte, não o limiar.
//!
//! O limiar verdadeiro é E[explicado]: um ser vive sse a energia que arranca
//! do resíduo supera o custo de existir. Em ruído puro esse valor define
//!
//!         met_c = c(M) · σ²
//!
//! TEORIA DOS VALORES EXTREMOS. Se o ser fica com a MELHOR de M geometrias
//! testadas contra aquela realização de ruído, então ⟨r,ψ̂⟩²/σ² é o máximo de
//! M variáveis χ²₁, e E[max] ≈ 2·ln(M) para M grande. Logo:
//!
//!         met_c ≈ 2σ²·ln(M_eff)                                    (LEI)
//!
//! PREDIÇÕES:
//!   Q1 — c é LINEAR em ln(orçamento de busca). Coeficiente angular ≈ 2.
//!   Q2 — c NÃO depende de σ (a divisão por σ² já absorve tudo).
//!   Q3 — met = c·σ² é o ponto de extinção: com 0.8·met_c a população
//!        persiste em ruído puro; com 1.25·met_c ela colapsa.

use agua_atoms::{Atom, AtomBounds, Gabor1D, Gabor2D};
use agua_core::{Field, Rng};
use agua_life::{LifeConfig, World};

fn noise_field(dims: &[usize], sigma: f32, seed: u64) -> Field {
    let mut rng = Rng::new(seed);
    let n: usize = dims.iter().product();
    Field { dims: dims.to_vec(), data: (0..n).map(|_| rng.normal() * sigma).collect() }
}

fn base_cfg(comp: usize, grad: bool) -> LifeConfig {
    let mut c = LifeConfig::default();
    c.cap = 600;
    c.init = 300;
    c.rain_every = 1;
    c.rain_count = 4;
    c.rain_competition = comp;
    c.tune_every = 1;
    c.gradient = grad;
    c.metabolism = 1e-6; // ~zero: ninguém morre de fome, medimos o que capturam
    c.initial_energy = 50.0;
    c.repro_threshold = 1e9; // sem reprodução: isola a captura individual
    c
}

/// c = E[explicado dos vivos]/σ² em RUÍDO PURO, após `ticks` de busca.
fn c_dynamic<A: Atom>(
    dims: &[usize],
    sigma: f32,
    bounds: &AtomBounds,
    seed: u64,
    ticks: u64,
    comp: usize,
) -> f32 {
    c_dynamic_pop::<A>(dims, sigma, bounds, seed, ticks, comp, None)
}

/// Igual, mas permite fixar a DENSIDADE de seres (pop por pixel). Sem isso,
/// campos pequenos com população fixa têm o resíduo devorado pela competição
/// e o "explicado" cai — confound que arruína a medida de escala.
fn c_dynamic_pop<A: Atom>(
    dims: &[usize],
    sigma: f32,
    bounds: &AtomBounds,
    seed: u64,
    ticks: u64,
    comp: usize,
    pop: Option<usize>,
) -> f32 {
    let noise = noise_field(dims, sigma, seed);
    let mut cfg = base_cfg(comp, A::GEOM > 0);
    if let Some(p) = pop {
        cfg.cap = p.max(8);
        cfg.init = (p / 2).max(4);
        cfg.rain_count = (p / 150).max(1);
    }
    let mut w: World<A> = World::new(noise, cfg, bounds.clone(), seed);
    for _ in 0..ticks {
        w.step();
    }
    let e = w.explained_now();
    if e.is_empty() {
        return f32::NAN;
    }
    let mean = e.iter().sum::<f32>() / e.len() as f32;
    mean / (sigma * sigma)
}

/// Q3: teste de extinção. Devolve (pop com 0.8·met_c, pop com 1.25·met_c).
fn extinction_test<A: Atom>(
    dims: &[usize],
    sigma: f32,
    bounds: &AtomBounds,
    seed: u64,
    met_c: f32,
    comp: usize,
) -> (usize, usize) {
    let noise = noise_field(dims, sigma, seed);
    let run = |met: f32| -> usize {
        let mut cfg = base_cfg(comp, A::GEOM > 0);
        cfg.metabolism = met;
        cfg.initial_energy = 1.0; // energia inicial baixa: não mascara o limiar
        cfg.max_energy = 5.0;
        let mut w: World<A> = World::new(noise.clone(), cfg, bounds.clone(), seed);
        for _ in 0..600 {
            w.step();
        }
        w.beings.len()
    };
    (run(0.8 * met_c), run(1.25 * met_c))
}

fn linfit(xs: &[f32], ys: &[f32]) -> (f32, f32, f32) {
    let n = xs.len() as f32;
    let mx = xs.iter().sum::<f32>() / n;
    let my = ys.iter().sum::<f32>() / n;
    let sxy: f32 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f32 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let syy: f32 = ys.iter().map(|y| (y - my) * (y - my)).sum();
    let a = sxy / sxx;
    let b = my - a * mx;
    let r2 = if syy > 0.0 { sxy * sxy / (sxx * syy) } else { 1.0 };
    (a, b, r2)
}

pub fn run(seed: u64) {
    println!("=== LEI: met_c = c(M)·σ² , com c ≈ 2·ln(M_eff) ===\n");

    let b2 = AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
    let b1 = AtomBounds { freq_min: 0.002, freq_max: 0.45, sigma_min: 2.0, sigma_max: 400.0 };
    let d2 = [64usize, 64];
    let d1 = [4096usize];

    // ── Q2: c depende de σ? ──────────────────────────────────────────────
    println!("Q2 — c depende de σ? (predição: NÃO)");
    println!("{:<8} {:>8} {:>10} {:>10}", "atomo", "sigma", "c", "met_c");
    let sigmas = [0.02f32, 0.05, 0.10, 0.20, 0.40];
    let mut c2s = Vec::new();
    for &s in &sigmas {
        let c = c_dynamic::<Gabor2D>(&d2, s, &b2, seed, 200, 8);
        c2s.push(c);
        println!("{:<8} {:>8.3} {:>10.2} {:>10.5}", "Gabor2D", s, c, c * s * s);
    }
    let mean2 = c2s.iter().sum::<f32>() / c2s.len() as f32;
    let sd2 = (c2s.iter().map(|c| (c - mean2).powi(2)).sum::<f32>() / c2s.len() as f32).sqrt();
    println!("   c = {mean2:.2} ± {sd2:.2}  (CV = {:.1}%)\n", 100.0 * sd2 / mean2);

    // ── Q1: c é linear em ln(orçamento de busca)? ────────────────────────
    println!("Q1 — c vs ln(orçamento de busca)  (predição: reta, inclinação ≈ 2)");
    println!("{:<10} {:>8} {:>8} {:>10}", "atomo", "ticks", "comp", "c");
    let sigma = 0.1f32;
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for &ticks in &[5u64, 15, 50, 150, 400] {
        for &comp in &[1usize, 8, 64] {
            let c = c_dynamic::<Gabor2D>(&d2, sigma, &b2, seed, ticks, comp);
            // orçamento efetivo por ser ≈ (nº de propostas de geometria) × (competidores na vaga)
            let m_eff = (ticks as f32) * (comp as f32);
            xs.push(m_eff.ln());
            ys.push(c);
            println!("{:<10} {:>8} {:>8} {:>10.2}", "Gabor2D", ticks, comp, c);
        }
    }
    let (a, b, r2) = linfit(&xs, &ys);
    println!("   ajuste: c = {a:.3}·ln(M) + {b:.3}   (R² = {r2:.4})");
    println!("   → M_eff implícito em c={mean2:.1}: M ≈ {:.0}\n", (mean2 / 2.0).exp());

    // ── P3 (v1): graus de liberdade ──────────────────────────────────────
    println!("P3 — c cresce com os graus de liberdade do átomo?");
    let c2 = c_dynamic::<Gabor2D>(&d2, sigma, &b2, seed, 200, 8);
    let c1 = c_dynamic::<Gabor1D>(&d1, sigma, &b1, seed, 200, 8);
    println!("   Gabor2D (6 geom): c = {c2:.2}");
    println!("   Gabor1D (4 geom): c = {c1:.2}");
    println!("   razão c2/c1 = {:.2}  (graus 6/4 = 1.50)\n", c2 / c1);

    // ── Q4: LIGAÇÃO COM DONOHO–JOHNSTONE ────────────────────────────────
    // Limiar universal de wavelet: λ = σ·√(2·ln N). Em ENERGIA: λ² = 2σ²·ln N.
    // Nossa lei: met_c = c·σ². Se c = 2·ln(N_pixels), então o metabolismo
    // crítico É o limiar universal. Predição: c cresce com ln(N) com
    // inclinação 2 (N = nº de lugares onde o ruído pode conspirar).
    println!("Q4 — c vs ln(tamanho do campo)   (Donoho-Johnstone: inclinacao = 2)");
    println!("{:<10} {:>10} {:>10} {:>12} {:>12}", "lado", "N_pixels", "c medido", "2*ln(N)", "razao");
    let mut xs4 = Vec::new();
    let mut ys4 = Vec::new();
    // DENSIDADE FIXA: 1 ser a cada 32 pixels. Sem isso, campo pequeno = seres
    // demais por pixel = resíduo devorado = "explicado" artificialmente baixo.
    for &lado in &[24usize, 32, 48, 64, 96, 128] {
        let dims = [lado, lado];
        let pop = (lado * lado) / 32;
        let c = c_dynamic_pop::<Gabor2D>(&dims, 0.1, &b2, seed, 200, 8, Some(pop));
        let npx = (lado * lado) as f32;
        let dj = 2.0 * npx.ln();
        xs4.push(npx.ln());
        ys4.push(c);
        println!("{:<10} {:>10.0} {:>10.2} {:>12.2} {:>12.2}  (pop={pop})", lado, npx, c, dj, c / dj);
    }
    let (a4, b4, r24) = linfit(&xs4, &ys4);
    println!("   ajuste: c = {a4:.3}*ln(N) + {b4:.3}   (R2 = {r24:.4})");
    println!("   predicao Donoho-Johnstone: inclinacao 2.000\n");

    // ── Q3: o ponto de extinção é mesmo em met_c? ────────────────────────
    println!("Q3 — extinção em ruído puro (predição: vive a 0.8·met_c, morre a 1.25·met_c)");
    println!("{:<8} {:>8} {:>10} {:>14} {:>14}", "atomo", "sigma", "met_c", "pop@0.8·met_c", "pop@1.25·met_c");
    for &s in &[0.05f32, 0.1, 0.2] {
        let c = c_dynamic::<Gabor2D>(&d2, s, &b2, seed, 200, 8);
        let met_c = c * s * s;
        let (lo, hi) = extinction_test::<Gabor2D>(&d2, s, &b2, seed, met_c, 8);
        println!("{:<8} {:>8.3} {:>10.5} {:>14} {:>14}", "Gabor2D", s, met_c, lo, hi);
    }
}
