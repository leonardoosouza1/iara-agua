//! `agua-lab law3` — AS LEIS DA HARMONIA E DA CAÇA AO RUÍDO.
//!
//! Três leis candidatas, cada uma com um análogo conhecido na literatura:
//!
//! H1 — EQUIPARTIÇÃO (preço sombra / KKT). No equilíbrio, a utilidade marginal
//!      iguala o custo marginal: todo ser vivo explica ≈ met. Se explicasse
//!      muito mais, vizinhos nasceriam e dividiriam o resíduo; se menos, morre.
//!      Predição: E[explicado] = met, com inclinação 1 na reta E[expl] vs met.
//!
//! H2 — HARMONIA = DECOERÊNCIA (coerência mútua, Donoho & Elad 2003). A
//!      "harmonia" entre seres é literalmente a coerência μ = |⟨ψ̂ᵢ,ψ̂ⱼ⟩| entre
//!      átomos normalizados. Átomos que se pisam competem pelo mesmo resíduo.
//!      Predição: a população SOBREVIVENTE tem coerência menor que uma
//!      população aleatória de mesmo tamanho. Ninguém programou ortogonalidade.
//!
//! H3 — CAÇA AO RUÍDO = MOROZOV. O princípio da discrepância diz: pare de
//!      ajustar quando ‖r‖²/n ≈ σ². Predição: com met = met_ótimo, a população
//!      para SOZINHA no nível de ruído — sem ninguém lhe dizer σ.

use agua_atoms::{Atom, AtomBounds, Gabor2D};
use agua_core::{Field, Rng};
use agua_life::{LifeConfig, World};

fn synth(size: usize) -> Field {
    let s = size as f32;
    let mut f = Field::from_fn_2d(size, size, |x, y| {
        let (xf, yf) = (x as f32, y as f32);
        let blob = 0.5
            * (-((xf - s * 0.35).powi(2) + (yf - s * 0.35).powi(2)) / (2.0 * (s * 0.13).powi(2)))
                .exp();
        let stripes = if xf > s * 0.55 {
            0.25 * (std::f32::consts::TAU * 0.15 * xf).sin()
        } else {
            0.0
        };
        blob + stripes
    });
    let m = f.mean();
    for v in f.data.iter_mut() {
        *v -= m;
    }
    f
}

fn add_noise(f: &Field, sigma: f32, seed: u64) -> Field {
    let mut rng = Rng::new(seed);
    let mut o = f.clone();
    for v in o.data.iter_mut() {
        *v += rng.normal() * sigma;
    }
    o
}

fn cfg(met: f32) -> LifeConfig {
    let mut c = LifeConfig::default();
    c.cap = 1500;
    c.init = 64;
    c.metabolism = met;
    c.rain_every = 1;
    c.rain_count = 6;
    c.rain_competition = 16;
    c.tune_every = 1;
    c.gradient = true;
    c
}

/// Coerência mútua média entre pares de átomos com suporte sobreposto.
/// μ = |⟨ψ̂ᵢ,ψ̂ⱼ⟩|. Só conta pares que de fato se tocam (senão é 0 trivial).
fn coherence<A: Atom>(units: &[agua_core::Patch]) -> (f32, f32, usize) {
    let mut sum = 0.0f64;
    let mut max: f32 = 0.0;
    let mut pares = 0usize;
    for i in 0..units.len() {
        let ni = units[i].norm2().sqrt();
        if ni < 1e-6 {
            continue;
        }
        for j in (i + 1)..units.len() {
            let nj = units[j].norm2().sqrt();
            if nj < 1e-6 {
                continue;
            }
            let d = units[i].dot_patch(&units[j]).abs() / (ni * nj);
            if d > 1e-6 {
                sum += d as f64;
                max = max.max(d);
                pares += 1;
            }
        }
    }
    let media = if pares > 0 { (sum / pares as f64) as f32 } else { 0.0 };
    (media, max, pares)
}

pub fn run(seed: u64) {
    let size = 64usize;
    let dims = [size, size];
    let b = AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
    let clean = synth(size);

    // ── H1: equipartição ────────────────────────────────────────────────
    println!("=== H1 — EQUIPARTIÇÃO: no equilíbrio, todo ser explica ≈ met ===");
    println!("(utilidade marginal = custo marginal; análogo do preço sombra / KKT)\n");
    println!("{:>10} {:>14} {:>10} {:>12} {:>8}", "met", "E[explicado]", "razão", "desvio/média", "vivos");
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for &met in &[0.02f32, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6] {
        let mut w: World<Gabor2D> = World::new(clean.clone(), cfg(met), b.clone(), seed);
        for _ in 0..900 {
            w.step();
        }
        let e = w.explained_now();
        if e.is_empty() {
            println!("{met:>10.3} {:>14} {:>10} {:>12} {:>8}", "-", "-", "-", 0);
            continue;
        }
        let n = e.len() as f32;
        let mean = e.iter().sum::<f32>() / n;
        let sd = (e.iter().map(|v| (v - mean).powi(2)).sum::<f32>() / n).sqrt();
        xs.push(met);
        ys.push(mean);
        println!(
            "{met:>10.3} {mean:>14.4} {:>10.2} {:>12.2} {:>8}",
            mean / met,
            sd / mean,
            e.len()
        );
    }
    // inclinação de E[expl] vs met (predição: 1.0)
    let n = xs.len() as f32;
    let mx = xs.iter().sum::<f32>() / n;
    let my = ys.iter().sum::<f32>() / n;
    let sxy: f32 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f32 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    println!("\n  inclinação E[explicado] vs met = {:.3}   (predição: 1.000)\n", sxy / sxx);

    // ── H2: harmonia = decoerência ──────────────────────────────────────
    println!("=== H2 — HARMONIA: a população sobrevivente é menos COERENTE? ===");
    println!("(coerência mútua μ = |⟨ψ̂ᵢ,ψ̂ⱼ⟩|, Donoho & Elad; menor = menos se pisam)\n");
    for &met in &[0.05f32, 0.2, 0.8] {
        let mut w: World<Gabor2D> = World::new(clean.clone(), cfg(met), b.clone(), seed);
        for _ in 0..900 {
            w.step();
        }
        let units: Vec<_> = w.beings.iter().map(|x| x.unit.clone()).collect();
        let (mu_vivo, max_vivo, p1) = coherence::<Gabor2D>(&units);

        // controle: MESMO número de átomos, geometrias aleatórias
        let mut rng = Rng::new(seed ^ 0xC0FFEE);
        let rnd: Vec<_> = (0..units.len())
            .map(|_| {
                let pos = [rng.range(0.0, size as f32 - 1.0), rng.range(0.0, size as f32 - 1.0)];
                Gabor2D::rain_at(&pos, &b, &mut rng).render_unit(&dims)
            })
            .collect();
        let (mu_rnd, max_rnd, p2) = coherence::<Gabor2D>(&rnd);
        println!(
            "met={met:<5} N={:<5} μ_médio vivo={mu_vivo:.4} (max {max_vivo:.3}, {p1} pares)  |  μ_médio aleatório={mu_rnd:.4} (max {max_rnd:.3}, {p2} pares)  →  {:+.1}%",
            units.len(),
            100.0 * (mu_vivo - mu_rnd) / mu_rnd
        );
    }

    // ── H3: caça ao ruído = Morozov ─────────────────────────────────────
    println!("\n=== H3 — CAÇA AO RUÍDO: a população para no nível do ruído? ===");
    println!("(princípio da discrepância de Morozov: pare quando ‖r‖²/n ≈ σ²)");
    println!("met vem da LEI (40·σ²); ninguém informa σ à população.\n");
    println!("{:>8} {:>10} {:>16} {:>16} {:>10} {:>8}", "sigma", "met=40σ²", "‖r‖²/n final", "σ² (alvo)", "razão", "vivos");
    for &s in &[0.05f32, 0.10, 0.15, 0.20, 0.30] {
        let noisy = add_noise(&clean, s, seed ^ 0x0157);
        let met = 40.0 * s * s;
        let mut w: World<Gabor2D> = World::new(noisy, cfg(met), b.clone(), seed);
        for _ in 0..900 {
            w.step();
        }
        w.rebuild();
        let disc = w.residual.energy() / w.residual.len() as f32;
        println!(
            "{s:>8.2} {met:>10.4} {disc:>16.5} {:>16.5} {:>10.2} {:>8}",
            s * s,
            disc / (s * s),
            w.beings.len()
        );
    }
    println!("\n  Morozov: razão ≈ 1.0 significa que a população parou exatamente no ruído.");

    // ── H4: PREDAÇÃO POR NICHO (exclusão competitiva de Gause) ──────────
    println!("\n=== H4 — PREDAÇÃO: seres do mesmo nicho se devoram. A harmonia paga? ===");
    println!("(H2 mostrou que a população NÃO se ortogonaliza sozinha: μ_vivo > μ_aleatório)");
    println!("θ = limiar de coerência acima do qual o mais fraco é devorado.\n");
    println!("--- reconstrução (sinal limpo, met=0.2) ---");
    println!("{:>8} {:>8} {:>10} {:>10} {:>10} {:>12}", "θ", "vivos", "PSNR", "μ_médio", "μ_max", "PSNR/átomo");
    for &theta in &[0.0f32, 0.95, 0.85, 0.7, 0.5] {
        let mut c = cfg(0.2);
        c.predation = theta;
        c.predation_every = 5;
        let mut w: World<Gabor2D> = World::new(clean.clone(), c, b.clone(), seed);
        for _ in 0..900 {
            w.step();
        }
        w.rebuild();
        let units: Vec<_> = w.beings.iter().map(|x| x.unit.clone()).collect();
        let (mu, mx, _) = coherence::<Gabor2D>(&units);
        let p = w.psnr(1.0);
        let lbl = if theta == 0.0 { "off".to_string() } else { format!("{theta:.2}") };
        println!(
            "{lbl:>8} {:>8} {p:>10.2} {mu:>10.4} {mx:>10.3} {:>12.4}",
            w.beings.len(),
            p / w.beings.len().max(1) as f32
        );
    }

    println!("\n--- denoise (σ=0.15, met=40σ²=0.9; PSNR contra o LIMPO) ---");
    println!("{:>8} {:>8} {:>10} {:>10} {:>16}", "θ", "vivos", "PSNR", "μ_médio", "‖r‖²/n ÷ σ²");
    let sg = 0.15f32;
    let noisy = add_noise(&clean, sg, seed ^ 0x0157);
    for &theta in &[0.0f32, 0.95, 0.85, 0.7, 0.5] {
        let mut c = cfg(40.0 * sg * sg);
        c.predation = theta;
        c.predation_every = 5;
        let mut w: World<Gabor2D> = World::new(noisy.clone(), c, b.clone(), seed);
        for _ in 0..900 {
            w.step();
        }
        w.rebuild();
        let mse: f32 = w
            .recon
            .data
            .iter()
            .zip(&clean.data)
            .map(|(a, c)| (a - c) * (a - c))
            .sum::<f32>()
            / clean.len() as f32;
        let units: Vec<_> = w.beings.iter().map(|x| x.unit.clone()).collect();
        let (mu, _, _) = coherence::<Gabor2D>(&units);
        let disc = w.residual.energy() / w.residual.len() as f32;
        let lbl = if theta == 0.0 { "off".to_string() } else { format!("{theta:.2}") };
        println!(
            "{lbl:>8} {:>8} {:>10.2} {mu:>10.4} {:>16.2}",
            w.beings.len(),
            agua_core::psnr(mse, 1.0),
            disc / (sg * sg)
        );
    }
}
