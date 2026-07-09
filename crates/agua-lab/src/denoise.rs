//! `agua-lab denoise` — RESTAURAÇÃO POR SOBREVIVÊNCIA (zero-shot).
//!
//! Hipótese: ruído branco é INCOMPRESSÍVEL — nenhum átomo consegue explicar
//! energia suficiente de ruído pra pagar seu metabolismo. Estrutura, sim.
//! Logo a população adota a estrutura e deixa o ruído morrer de fome, SEM
//! treino, SEM dataset, SEM prior. O metabolismo vira um critério de parada
//! espacialmente adaptativo (mais átomos onde há estrutura, nenhum no ruído).
//!
//! É o único regime onde a água tem um ativo que o guloso não tem: o MP
//! precisa que você DIGA quantos átomos usar. Aqui isso é a pergunta.
//!
//! Baselines (todos contra o sinal LIMPO, que o solver nunca vê):
//!   - ruidoso        : não fazer nada (piso)
//!   - gauss(σ)       : filtro linear clássico, σ varrido, melhor = oráculo
//!   - MP@N-oráculo   : matching pursuit com o N que MAXIMIZA o PSNR limpo
//!                      (trapaça a favor dele: usa a imagem limpa pra escolher N)
//!   - MP@Morozov     : MP parando quando ‖r‖²/n ≤ σ² (o critério que um
//!                      engenheiro honesto usaria sem ver o limpo)
//!   - sobrevivência  : metabolismo fixo, escolhe N sozinha

use agua_atoms::{fit_amp, Atom, AtomBounds};
use agua_core::{psnr, unflatten, Field, HotspotSampler, Rng};
use agua_life::{LifeConfig, World};
use serde::Serialize;

#[derive(Serialize, Clone)]
pub struct DenoiseRow {
    pub metodo: String,
    pub psnr_limpo: f32,
    pub atomos: usize,
    pub nota: String,
}

fn mse_vs(a: &Field, b: &Field) -> f32 {
    a.data
        .iter()
        .zip(&b.data)
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f32>()
        / a.len() as f32
}

/// Estimador MAD do desvio-padrão do ruído, SEM ver o sinal limpo.
/// Clássico (Donoho & Johnstone): num detalhe de alta-frequência, a mediana
/// do valor absoluto é robusta à estrutura do sinal; σ̂ = MAD/0.6745.
/// Detalhe = diferença de vizinhos / √2 (2D: diagonal, que mata gradientes
/// lineares em ambos os eixos).
pub fn estimate_sigma(f: &Field) -> f32 {
    let mut d: Vec<f32> = Vec::new();
    if f.dims.len() == 2 {
        let (h, w) = (f.dims[0], f.dims[1]);
        for y in 0..h - 1 {
            for x in 0..w - 1 {
                // laplaciano diagonal: (a - b - c + d)/2 tem variância σ²
                let a = f.data[y * w + x];
                let b = f.data[y * w + x + 1];
                let c = f.data[(y + 1) * w + x];
                let e = f.data[(y + 1) * w + x + 1];
                d.push((a - b - c + e) * 0.5);
            }
        }
    } else {
        for i in 0..f.len() - 1 {
            d.push((f.data[i + 1] - f.data[i]) / std::f32::consts::SQRT_2);
        }
    }
    for v in d.iter_mut() {
        *v = v.abs();
    }
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = d[d.len() / 2];
    med / 0.6745
}

/// Ruído gaussiano determinístico.
pub fn add_noise(clean: &Field, sigma: f32, seed: u64) -> Field {
    let mut rng = Rng::new(seed);
    let mut f = clean.clone();
    for v in f.data.iter_mut() {
        *v += rng.normal() * sigma;
    }
    f
}

/// Blur gaussiano separável (baseline clássico), só 2D.
fn gaussian_blur(f: &Field, sigma: f32) -> Field {
    if f.dims.len() != 2 {
        return f.clone();
    }
    let (h, w) = (f.dims[0], f.dims[1]);
    let r = (3.0 * sigma).ceil() as i64;
    let k: Vec<f32> = (-r..=r)
        .map(|i| (-(i * i) as f32 / (2.0 * sigma * sigma)).exp())
        .collect();
    let sum: f32 = k.iter().sum();
    let k: Vec<f32> = k.iter().map(|v| v / sum).collect();
    let mut tmp = Field::zeros(&f.dims);
    for y in 0..h {
        for x in 0..w {
            let mut acc = 0.0;
            for (i, kv) in k.iter().enumerate() {
                let xx = (x as i64 + i as i64 - r).clamp(0, w as i64 - 1) as usize;
                acc += kv * f.data[y * w + xx];
            }
            tmp.data[y * w + x] = acc;
        }
    }
    let mut out = Field::zeros(&f.dims);
    for y in 0..h {
        for x in 0..w {
            let mut acc = 0.0;
            for (i, kv) in k.iter().enumerate() {
                let yy = (y as i64 + i as i64 - r).clamp(0, h as i64 - 1) as usize;
                acc += kv * tmp.data[yy * w + x];
            }
            out.data[y * w + x] = acc;
        }
    }
    out
}

/// MP no sinal RUIDOSO, avaliando cada passo contra o LIMPO.
/// Devolve (N_oráculo, psnr_oráculo, N_morozov, psnr_morozov).
fn mp_denoise<A: Atom>(
    noisy: &Field,
    clean: &Field,
    dc: f32,
    n_max: usize,
    candidates: usize,
    bounds: &AtomBounds,
    seed: u64,
    sigma: f32,
) -> (usize, f32, usize, f32) {
    let mut rng = Rng::new(seed);
    let mut residual = noisy.clone();
    let mut recon = Field::zeros(&noisy.dims);
    let eval = |r: &Field| {
        let mut o = r.clone();
        for v in o.data.iter_mut() {
            *v += dc;
        }
        psnr(mse_vs(&o, clean), 1.0)
    };
    let (mut best_n, mut best_p) = (0usize, eval(&recon));
    let (mut moro_n, mut moro_p) = (0usize, best_p);
    let mut moro_done = false;
    let sigma2 = sigma * sigma;

    for k in 1..=n_max {
        let sampler = HotspotSampler::new(&residual);
        let mut best: Option<(f32, f32, agua_core::Patch)> = None;
        for _ in 0..candidates {
            let idx = sampler.sample(&mut rng);
            let pos: Vec<f32> = unflatten(&noisy.dims, idx).into_iter().map(|v| v as f32).collect();
            let g = A::rain_at(&pos, bounds, &mut rng);
            let unit = g.render_unit(&noisy.dims);
            if let Some((amp, expl)) = fit_amp(&unit, &residual) {
                if best.as_ref().map_or(true, |(e, _, _)| expl > *e) {
                    best = Some((expl, amp, unit));
                }
            }
        }
        let Some((_, amp, unit)) = best else { continue };
        residual.add_scaled_patch(&unit, -amp);
        recon.add_scaled_patch(&unit, amp);

        let p = eval(&recon);
        if p > best_p {
            best_p = p;
            best_n = k;
        }
        // Morozov: parar quando o resíduo atinge o nível de ruído
        if !moro_done && residual.energy() / residual.len() as f32 <= sigma2 {
            moro_done = true;
            moro_n = k;
            moro_p = p;
        }
    }
    if !moro_done {
        moro_n = n_max;
        moro_p = eval(&recon);
    }
    (best_n, best_p, moro_n, moro_p)
}

pub struct DenoiseArgs<'a> {
    pub clean: Field,
    pub sigma: f32,
    pub life: LifeConfig,
    pub bounds: AtomBounds,
    pub ticks: u64,
    pub seed: u64,
    pub mp_candidates: usize,
    pub n_max: usize,
    /// só 2D: roda o baseline de blur
    pub blur: bool,
    pub label: &'a str,
}

/// Roda o experimento e devolve as linhas da tabela (+ a recon da vida).
pub fn run<A: Atom>(a: DenoiseArgs) -> (Vec<DenoiseRow>, Field, Field) {
    let noisy = add_noise(&a.clean, a.sigma, a.seed ^ 0x0157);
    // átomos têm média ~0 ⇒ o DC é side-info. Estimado do RUIDOSO (o solver
    // nunca vê o limpo); E[média(ruidoso)] = média(limpo), estimador correto.
    let dc = noisy.mean();
    let mut noisy0 = noisy.clone();
    for v in noisy0.data.iter_mut() {
        *v -= dc;
    }
    let with_dc = |f: &Field| {
        let mut o = f.clone();
        for v in o.data.iter_mut() {
            *v += dc;
        }
        o
    };
    let mut rows = Vec::new();

    rows.push(DenoiseRow {
        metodo: "ruidoso (nada)".into(),
        psnr_limpo: psnr(mse_vs(&noisy, &a.clean), 1.0),
        atomos: 0,
        nota: format!("sigma={}", a.sigma),
    });

    if a.blur && a.clean.dims.len() == 2 {
        let (mut bp, mut bs) = (f32::MIN, 0.0f32);
        for s in [0.5f32, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5] {
            let p = psnr(mse_vs(&gaussian_blur(&noisy, s), &a.clean), 1.0);
            if p > bp {
                bp = p;
                bs = s;
            }
        }
        rows.push(DenoiseRow {
            metodo: "gauss@oráculo".into(),
            psnr_limpo: bp,
            atomos: 0,
            nota: format!("sigma_blur={bs} escolhido vendo o limpo"),
        });
    }

    // sobrevivência: o mundo é o sinal RUIDOSO; ela nunca vê o limpo
    let mut world: World<A> = World::new(noisy0.clone(), a.life, a.bounds.clone(), a.seed);
    for _ in 0..a.ticks {
        world.step();
    }
    world.rebuild();
    let vida_recon = with_dc(&world.recon);
    rows.push(DenoiseRow {
        metodo: "SOBREVIVÊNCIA".into(),
        psnr_limpo: psnr(mse_vs(&vida_recon, &a.clean), 1.0),
        atomos: world.beings.len(),
        nota: format!("N emergiu do metabolismo; evals={}", world.evals),
    });

    let n_max = a.n_max.max(world.beings.len() * 2);
    let (n_or, p_or, n_mo, p_mo) =
        mp_denoise::<A>(&noisy0, &a.clean, dc, n_max, a.mp_candidates, &a.bounds, a.seed ^ 0x1234, a.sigma);
    rows.push(DenoiseRow {
        metodo: "MP@N-oráculo".into(),
        psnr_limpo: p_or,
        atomos: n_or,
        nota: "N escolhido VENDO o limpo (trapaça a favor dele)".into(),
    });
    rows.push(DenoiseRow {
        metodo: "MP@Morozov".into(),
        psnr_limpo: p_mo,
        atomos: n_mo,
        nota: "para quando ‖r‖²/n ≤ σ² (critério honesto sem ver o limpo)".into(),
    });

    (rows, vida_recon, noisy)
}
