//! `agua-lab law` — MEDINDO A LEI DO METABOLISMO CRÍTICO.
//!
//! TEORIA. Um átomo ajustado sobre ruído branco N(0,σ²) captura energia
//! espúria. A energia explicada por um átomo é E = ⟨r,ψ̂⟩² com ψ̂ unitário.
//! Para ruído puro e ψ̂ FIXO (sem busca), ⟨r,ψ̂⟩ ~ N(0,σ²) ⇒ E[E] = σ².
//! Mas o ser OTIMIZA sua geometria contra o ruído: ele escolhe, entre muitas
//! geometrias, a que mais casa com aquela realização do ruído. Isso infla a
//! energia capturada por um FATOR DE BUSCA c:
//!
//!         E[explicado | ruído, após busca] = c · σ²
//!
//! Um ser sobrevive sse explicado > metabolismo. Logo existe um metabolismo
//! CRÍTICO abaixo do qual até o ruído paga a conta de existir:
//!
//!         met_c = c · σ²                                        (LEI 1)
//!
//! PREDIÇÕES FALSIFICÁVEIS:
//!   P1 — expoente: met_c ∝ σ^α com α = 2 (medimos α por regressão log-log).
//!   P2 — universalidade: c NÃO depende do sinal (é propriedade do átomo+busca).
//!   P3 — c cresce com os graus de liberdade da busca: Gabor1D (4 geom) deve
//!        ter c MENOR que Gabor2D (6 geom). Se c for igual, a teoria do fator
//!        de busca está errada.
//!
//! Medimos c de duas formas independentes (têm que bater):
//!   (a) DIRETO: energia explicada média de átomos em ruído puro, sem busca
//!       (deve dar ≈ σ², ancorando a escala) e COM busca (dá c·σ²).
//!   (b) CRÍTICO: varre o metabolismo em ruído puro e acha onde a população
//!       colapsa. met_c/σ² = c.

use agua_atoms::{fit_amp, Atom, AtomBounds};
use agua_core::{Field, Rng};
use agua_life::{LifeConfig, World};
use serde::Serialize;

#[derive(Serialize)]
pub struct LawRow {
    pub atomo: String,
    pub sigma: f32,
    /// E[explicado]/σ² para átomos SEM busca (âncora: deve ser ≈ 1)
    pub c_sem_busca: f32,
    /// E[explicado]/σ² após refino por gradiente (o fator de busca)
    pub c_com_busca: f32,
    /// met_c medido pelo colapso da população, dividido por σ²
    pub c_critico: f32,
    pub met_c: f32,
    pub pop_em_met_c: f32,
}

fn noise_field(dims: &[usize], sigma: f32, seed: u64) -> Field {
    let mut rng = Rng::new(seed);
    let n: usize = dims.iter().product();
    Field { dims: dims.to_vec(), data: (0..n).map(|_| rng.normal() * sigma).collect() }
}

/// (a) energia explicada média por átomos aleatórios em ruído puro,
/// sem busca e depois de `refine` passos de gradiente.
fn measured_c<A: Atom>(
    dims: &[usize],
    sigma: f32,
    bounds: &AtomBounds,
    seed: u64,
    samples: usize,
    refine: usize,
) -> (f32, f32) {
    let noise = noise_field(dims, sigma, seed);
    let mut rng = Rng::new(seed ^ 0xBEEF);
    let (mut sum_raw, mut sum_ref, mut n_ok) = (0.0f64, 0.0f64, 0usize);

    for _ in 0..samples {
        let pos: Vec<f32> = dims
            .iter()
            .map(|d| rng.range(0.0, (*d - 1) as f32))
            .collect();
        let mut g = A::rain_at(&pos, bounds, &mut rng);
        let unit = g.render_unit(dims);
        let Some((_, expl0)) = fit_amp(&unit, &noise) else { continue };
        sum_raw += expl0 as f64;

        // refino por gradiente contra o MESMO ruído (é isso que infla a captura)
        if A::GEOM > 0 {
            let (mut m, mut v) = (vec![0.0f32; A::GEOM], vec![0.0f32; A::GEOM]);
            for t in 1..=refine {
                let unit = g.render_unit(dims);
                let Some((amp, _)) = fit_amp(&unit, &noise) else { break };
                // resíduo-sem-ele = ruído (o átomo não está subtraído do campo)
                let mut resid = noise.clone();
                resid.add_scaled_patch(&unit, -amp);
                let gr = g.grad(dims, &resid, amp);
                let sc = g.param_scale(dims);
                let (b1, b2, eps) = (0.9f32, 0.999f32, 1e-8f32);
                let tf = t as f32;
                let mut delta = vec![0.0f32; A::GEOM];
                for k in 0..A::GEOM {
                    m[k] = b1 * m[k] + (1.0 - b1) * gr[k];
                    v[k] = b2 * v[k] + (1.0 - b2) * gr[k] * gr[k];
                    let mh = m[k] / (1.0 - b1.powf(tf));
                    let vh = v[k] / (1.0 - b2.powf(tf));
                    delta[k] = -0.03 * sc[k] * mh / (vh.sqrt() + eps);
                }
                g.apply_delta(&delta, bounds, dims);
            }
        }
        let unit = g.render_unit(dims);
        if let Some((_, expl1)) = fit_amp(&unit, &noise) {
            sum_ref += expl1 as f64;
            n_ok += 1;
        }
    }
    let s2 = (sigma * sigma) as f64;
    let n = n_ok.max(1) as f64;
    ((sum_raw / n / s2) as f32, (sum_ref / n / s2) as f32)
}

/// (b) metabolismo crítico: em RUÍDO PURO, varre met e acha onde a população
/// colapsa (< 5% do que sobrevive com met→0). Busca binária em log.
fn critical_met<A: Atom>(
    dims: &[usize],
    sigma: f32,
    bounds: &AtomBounds,
    seed: u64,
    ticks: u64,
    base: &LifeConfig,
) -> (f32, f32) {
    let noise = noise_field(dims, sigma, seed);
    let pop_at = |met: f32| -> f32 {
        let mut cfg = base.clone();
        cfg.metabolism = met;
        let mut w: World<A> = World::new(noise.clone(), cfg, bounds.clone(), seed);
        for _ in 0..ticks {
            w.step();
        }
        w.beings.len() as f32
    };
    // colapso = população cai abaixo de 5% do cap
    let alvo = base.cap as f32 * 0.05;
    let (mut lo, mut hi) = (0.02f32 * sigma * sigma, 400.0f32 * sigma * sigma);
    // garante bracket: lo sobrevive, hi colapsa
    if pop_at(hi) > alvo {
        return (hi, pop_at(hi));
    }
    for _ in 0..9 {
        let m = (lo * hi).sqrt(); // bisseção na média GEOMÉTRICA (a lei é em potência)
        if pop_at(m) > alvo {
            lo = m;
        } else {
            hi = m;
        }
    }
    let met_c = (lo * hi).sqrt();
    (met_c, pop_at(met_c))
}

pub fn run_for<A: Atom>(
    nome: &str,
    dims: &[usize],
    bounds: &AtomBounds,
    sigmas: &[f32],
    seed: u64,
) -> Vec<LawRow> {
    let mut base = LifeConfig::default();
    base.cap = 800;
    base.init = 200;
    base.rain_every = 1;
    base.rain_count = 4;
    base.rain_competition = 8;
    base.tune_every = 1;
    base.gradient = A::GEOM > 0;

    let mut rows = Vec::new();
    for &sigma in sigmas {
        let (c_raw, c_ref) = measured_c::<A>(dims, sigma, bounds, seed, 220, 30);
        let (met_c, pop) = critical_met::<A>(dims, sigma, bounds, seed, 260, &base);
        rows.push(LawRow {
            atomo: nome.into(),
            sigma,
            c_sem_busca: c_raw,
            c_com_busca: c_ref,
            c_critico: met_c / (sigma * sigma),
            met_c,
            pop_em_met_c: pop,
        });
        let r = rows.last().unwrap();
        println!(
            "{:<8} σ={:<6.3} c_sem_busca={:>6.2}  c_com_busca={:>7.2}  met_c={:<9.4} c_crit={:>7.2}",
            r.atomo, r.sigma, r.c_sem_busca, r.c_com_busca, r.met_c, r.c_critico
        );
    }
    rows
}

/// Regressão log-log: met_c = A·σ^α. Devolve (α, A, R²).
pub fn fit_exponent(rows: &[LawRow]) -> (f32, f32, f32) {
    let n = rows.len() as f32;
    let xs: Vec<f32> = rows.iter().map(|r| r.sigma.ln()).collect();
    let ys: Vec<f32> = rows.iter().map(|r| r.met_c.ln()).collect();
    let mx = xs.iter().sum::<f32>() / n;
    let my = ys.iter().sum::<f32>() / n;
    let sxy: f32 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f32 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let alpha = sxy / sxx;
    let a = (my - alpha * mx).exp();
    let syy: f32 = ys.iter().map(|y| (y - my) * (y - my)).sum();
    let r2 = if syy > 0.0 { sxy * sxy / (sxx * syy) } else { 1.0 };
    (alpha, a, r2)
}
