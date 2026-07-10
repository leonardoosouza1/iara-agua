//! `agua-lab law9` — γ É UNIVERSAL? (predição cross-modal do k do denoise)
//!
//! O denoise usa met = k·σ̂², com k = 40 calibrado em IMAGEM. De onde vem 40?
//! Da Lei 1, met_c = c·σ², e do fato empírico de que o ótimo de PSNR fica acima
//! do crítico: met_ótimo = γ·met_c ⇒ k = γ·c.
//!
//! Se γ for uma constante UNIVERSAL (mesma física, modalidade diferente), então
//!
//!         k_áudio = k_imagem · (c_áudio / c_imagem)
//!
//! c é medido em RUÍDO PURO, na mesma config do denoise. Nada disso olha o sinal
//! limpo. Se o k ótimo medido em áudio bater a predição, γ é universal e a
//! calibração some: basta medir c uma vez por átomo.
//!
//! Também testamos se Morozov (Lei 4) vale em áudio.

use agua_atoms::{Atom, AtomBounds, Gabor1D, Gabor2D};
use agua_core::{psnr, Field, Rng};
use agua_life::{LifeConfig, World};

fn noise_field(dims: &[usize], sigma: f32, seed: u64) -> Field {
    let mut rng = Rng::new(seed);
    let n: usize = dims.iter().product();
    Field { dims: dims.to_vec(), data: (0..n).map(|_| rng.normal() * sigma).collect() }
}

fn add_noise(f: &Field, sigma: f32, seed: u64) -> Field {
    let mut rng = Rng::new(seed);
    let mut o = f.clone();
    for v in o.data.iter_mut() {
        *v += rng.normal() * sigma;
    }
    o
}

/// EXATAMENTE a config do `denoise --auto`
fn denoise_cfg(met: f32) -> LifeConfig {
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

/// c = E[explicado]/σ² em ruído puro, na config do denoise (sem morte nem repro)
fn c_denoise<A: Atom>(dims: &[usize], b: &AtomBounds, seed: u64, ticks: u64) -> f32 {
    let sigma = 0.1f32;
    let campo = noise_field(dims, sigma, seed);
    let mut cfg = denoise_cfg(1e-6);
    cfg.cap = 500;
    cfg.init = 300;
    cfg.initial_energy = 50.0;
    cfg.repro_threshold = 1e9;
    let mut w: World<A> = World::new(campo, cfg, b.clone(), seed);
    for _ in 0..ticks {
        w.step();
    }
    let e = w.explained_now();
    (e.iter().sum::<f32>() / e.len().max(1) as f32) / (sigma * sigma)
}

fn mse(a: &Field, b: &Field) -> f32 {
    a.data.iter().zip(&b.data).map(|(x, y)| (x - y) * (x - y)).sum::<f32>() / a.len() as f32
}

/// varre k, devolve (k_ótimo, psnr_ótimo, razão de Morozov no ótimo)
fn melhor_k<A: Atom>(
    clean: &Field,
    b: &AtomBounds,
    sigma: f32,
    seed: u64,
    ticks: u64,
    ks: &[f32],
) -> (f32, f32, f32) {
    let noisy = add_noise(clean, sigma, seed ^ 0x0157);
    let dc = noisy.mean();
    let mut noisy0 = noisy.clone();
    for v in noisy0.data.iter_mut() {
        *v -= dc;
    }
    let (mut bk, mut bp, mut bm) = (0.0f32, f32::MIN, 0.0f32);
    for &k in ks {
        let mut w: World<A> = World::new(noisy0.clone(), denoise_cfg(k * sigma * sigma), b.clone(), seed);
        for _ in 0..ticks {
            w.step();
        }
        w.rebuild();
        let mut rec = w.recon.clone();
        for v in rec.data.iter_mut() {
            *v += dc;
        }
        let p = psnr(mse(&rec, clean), 1.0);
        let moro = (w.residual.energy() / w.residual.len() as f32) / (sigma * sigma);
        if p > bp {
            bp = p;
            bk = k;
            bm = moro;
        }
    }
    (bk, bp, bm)
}

fn sinal_audio(n: usize) -> Field {
    let nf = n as f32;
    let mut f = Field::from_fn_1d(n, |i| {
        let t = i as f32;
        let low = 0.35 * (std::f32::consts::TAU * 0.004 * t).sin();
        let mid = if t > nf * 0.3 && t < nf * 0.55 {
            0.30 * (std::f32::consts::TAU * 0.035 * t).sin()
        } else {
            0.0
        };
        let high = if t > nf * 0.75 { 0.20 * (std::f32::consts::TAU * 0.18 * t).sin() } else { 0.0 };
        low + mid + high
    });
    let m = f.mean();
    for v in f.data.iter_mut() {
        *v -= m;
    }
    f
}

fn sinal_imagem(size: usize) -> Field {
    let s = size as f32;
    let mut f = Field::from_fn_2d(size, size, |x, y| {
        let (xf, yf) = (x as f32, y as f32);
        let blob = 0.5
            * (-((xf - s * 0.35).powi(2) + (yf - s * 0.35).powi(2)) / (2.0 * (s * 0.13).powi(2))).exp();
        let stripes = if xf > s * 0.55 { 0.25 * (std::f32::consts::TAU * 0.15 * xf).sin() } else { 0.0 };
        blob + stripes
    });
    let m = f.mean();
    for v in f.data.iter_mut() {
        *v -= m;
    }
    f
}

pub fn run(seed: u64) {
    let b2 = AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
    let b1 = AtomBounds { freq_min: 0.002, freq_max: 0.45, sigma_min: 2.0, sigma_max: 400.0 };
    let d2 = [64usize, 64];
    let d1 = [8192usize];

    println!("=== γ É UNIVERSAL? predição cross-modal do k do denoise ===");
    println!("met_ótimo = γ·met_c = γ·c·σ²  ⇒  k = γ·c  ⇒  k_áudio = k_img·(c_áudio/c_img)\n");

    // c medido em RUÍDO PURO, config do denoise (não vê nenhum sinal limpo)
    let c_img = c_denoise::<Gabor2D>(&d2, &b2, seed, 800);
    let c_aud = c_denoise::<Gabor1D>(&d1, &b1, seed, 800);
    println!("c medido em ruído puro (config do denoise, 800 ticks):");
    println!("   imagem (Gabor2D): c = {c_img:.2}");
    println!("   áudio  (Gabor1D): c = {c_aud:.2}");
    println!("   razão c_áudio/c_img = {:.3}\n", c_aud / c_img);

    // grade FINA (a anterior era grossa demais pra cravar "γ não é universal")
    let ks: Vec<f32> = vec![8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 75.0, 90.0];

    // ── imagem: reproduz o k=40 e mede γ ────────────────────────────────
    println!("--- IMAGEM: varrendo k (σ=0.15) ---");
    let img = sinal_imagem(64);
    let (k_img, p_img, m_img) = melhor_k::<Gabor2D>(&img, &b2, 0.15, seed, 800, &ks);
    let gamma_img = k_img / c_img;
    println!("   k ótimo = {k_img}   PSNR = {p_img:.2} dB   Morozov = {m_img:.2}");
    println!("   γ = k/c = {gamma_img:.2}\n");

    // ── predição para áudio, ANTES de medir ─────────────────────────────
    let k_pred = k_img * (c_aud / c_img);
    println!("--- PREDIÇÃO: k_áudio = {k_img} × {:.3} = {k_pred:.1} ---", c_aud / c_img);

    println!("--- ÁUDIO: varrendo k (σ=0.15) ---");
    let aud = sinal_audio(8192);
    let (k_aud, p_aud, m_aud) = melhor_k::<Gabor1D>(&aud, &b1, 0.15, seed, 800, &ks);
    let gamma_aud = k_aud / c_aud;
    println!("   k ótimo = {k_aud}   PSNR = {p_aud:.2} dB   Morozov = {m_aud:.2}");
    println!("   γ = k/c = {gamma_aud:.2}\n");

    println!("=== VEREDITO ===");
    println!("k_áudio previsto {k_pred:.1}  |  medido {k_aud}  →  erro {:.0}%", 100.0 * (k_aud - k_pred).abs() / k_pred);
    println!("γ imagem {gamma_img:.2}  |  γ áudio {gamma_aud:.2}  →  {}",
        if (gamma_img - gamma_aud).abs() / gamma_img < 0.35 { "γ UNIVERSAL (dentro de 35%)" } else { "γ NÃO é universal" });
    println!("Morozov (‖r‖²/n ÷ σ² no ótimo): imagem {m_img:.2} | áudio {m_aud:.2}   (Lei 4 prevê ≈1)");
    println!("\ngrade de k: {ks:?}");

    // ── A REGRA QUE DISPENSA k E γ: escolher met por MOROZOV ─────────────
    // Se em ambas as modalidades o k ótimo cai onde ‖r‖²/n ≈ σ², então não é
    // preciso k nem γ: basta ajustar o metabolismo até o resíduo bater σ̂².
    // Nenhum hiperparâmetro, nenhuma calibração por modalidade.
    println!("\n=== A REGRA UNIVERSAL: escolher met por MOROZOV (sem k, sem γ) ===");
    println!("Bisseção no metabolismo até ‖r‖²/n = σ̂². Compara com o k-ótimo (oráculo).\n");
    println!("{:>10} {:>16} {:>16} {:>14} {:>10}", "modalidade", "PSNR k-ótimo", "PSNR Morozov", "met escolhido", "Δ (dB)");

    let p_img_mor = morozov::<Gabor2D>(&img, &b2, 0.15, seed, 800);
    println!("{:>10} {p_img:>16.2} {:>16.2} {:>14.4} {:>+10.2}", "imagem", p_img_mor.0, p_img_mor.1, p_img_mor.0 - p_img);
    let p_aud_mor = morozov::<Gabor1D>(&aud, &b1, 0.15, seed, 800);
    println!("{:>10} {p_aud:>16.2} {:>16.2} {:>14.4} {:>+10.2}", "áudio", p_aud_mor.0, p_aud_mor.1, p_aud_mor.0 - p_aud);
    println!("\nSe Δ ≈ 0, a regra de Morozov substitui a calibração de k — e é universal.\n");

    // ── de quem é a culpa: do estimador de σ ou da regra? ────────────────
    println!("=== ISOLANDO A CULPA: regra de Morozov com σ ORÁCULO ===");
    println!("Se o PSNR volta com σ verdadeiro, o culpado é o estimador (MAD viesado).");
    println!("Se não volta, a regra de Morozov é que não serve pra essa modalidade.\n");
    println!("{:>10} {:>14} {:>18} {:>18} {:>10}", "modalidade", "PSNR k-ótimo", "Morozov (σ̂ MAD)", "Morozov (σ oráculo)", "σ̂/σ");
    let im_or = morozov_sig::<Gabor2D>(&img, &b2, 0.15, seed, 800, Some(0.15));
    let au_or = morozov_sig::<Gabor1D>(&aud, &b1, 0.15, seed, 800, Some(0.15));
    let sh_i = estimar_sigma(&add_noise(&img, 0.15, seed ^ 0x0157)) / 0.15;
    let sh_a = estimar_sigma(&add_noise(&aud, 0.15, seed ^ 0x0157)) / 0.15;
    println!("{:>10} {p_img:>14.2} {:>18.2} {:>18.2} {sh_i:>10.3}", "imagem", p_img_mor.0, im_or.0);
    println!("{:>10} {p_aud:>14.2} {:>18.2} {:>18.2} {sh_a:>10.3}", "áudio", p_aud_mor.0, au_or.0);
}

/// bisseção no metabolismo até a discrepância bater σ̂² (σ̂ do próprio ruidoso).
/// Devolve (psnr_vs_limpo, met_escolhido).
fn morozov<A: Atom>(clean: &Field, b: &AtomBounds, sigma: f32, seed: u64, ticks: u64) -> (f32, f32) {
    morozov_sig::<A>(clean, b, sigma, seed, ticks, None)
}

/// `sigma_dado = Some(σ)` usa o σ VERDADEIRO (oráculo) — serve pra isolar
/// o erro do estimador do erro da própria regra.
fn morozov_sig<A: Atom>(
    clean: &Field,
    b: &AtomBounds,
    sigma: f32,
    seed: u64,
    ticks: u64,
    sigma_dado: Option<f32>,
) -> (f32, f32) {
    let noisy = add_noise(clean, sigma, seed ^ 0x0157);
    // σ̂ estimado do ruidoso (MAD do detalhe) — nunca olha o limpo
    let sig_hat = sigma_dado.unwrap_or_else(|| estimar_sigma(&noisy));
    let dc = noisy.mean();
    let mut noisy0 = noisy.clone();
    for v in noisy0.data.iter_mut() {
        *v -= dc;
    }
    let alvo = sig_hat * sig_hat;

    let rodar = |met: f32| -> (f32, f32) {
        let mut w: World<A> = World::new(noisy0.clone(), denoise_cfg(met), b.clone(), seed);
        for _ in 0..ticks {
            w.step();
        }
        w.rebuild();
        let disc = w.residual.energy() / w.residual.len() as f32;
        let mut rec = w.recon.clone();
        for v in rec.data.iter_mut() {
            *v += dc;
        }
        (disc, psnr(mse(&rec, clean), 1.0))
    };

    // met alto ⇒ poucos átomos ⇒ discrepância ALTA; met baixo ⇒ discrepância baixa
    let (mut lo, mut hi) = (0.01f32 * sigma * sigma, 400.0f32 * sigma * sigma);
    let mut melhor = (0.0f32, 0.0f32);
    for _ in 0..7 {
        let m = (lo * hi).sqrt();
        let (disc, p) = rodar(m);
        melhor = (p, m);
        if disc > alvo {
            hi = m; // ajustou de menos ⇒ baratear a vida
        } else {
            lo = m; // ajustou de mais (comeu ruído) ⇒ encarecer
        }
    }
    melhor
}

/// σ̂ por MAD do detalhe (Donoho–Johnstone), sem ver o limpo
fn estimar_sigma(f: &Field) -> f32 {
    let mut d: Vec<f32> = Vec::new();
    if f.dims.len() == 2 {
        let (h, w) = (f.dims[0], f.dims[1]);
        for y in 0..h - 1 {
            for x in 0..w - 1 {
                let a = f.data[y * w + x];
                let b = f.data[y * w + x + 1];
                let c = f.data[(y + 1) * w + x];
                let e = f.data[(y + 1) * w + x + 1];
                d.push(((a - b - c + e) * 0.5).abs());
            }
        }
    } else {
        for i in 0..f.len() - 1 {
            d.push(((f.data[i + 1] - f.data[i]) / std::f32::consts::SQRT_2).abs());
        }
    }
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    d[d.len() / 2] / 0.6745
}
