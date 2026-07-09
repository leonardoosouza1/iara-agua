//! `agua-lab law5` — UNIVERSALIDADE: a lei depende da FORMA do ruído?
//!
//! A lei met_c = c·σ² veio de teoria de valores extremos, e extremos dependem
//! da cauda. Ingenuamente: ruído Laplace (cauda pesada) deveria inflar c.
//!
//! MAS: o que o ser mede é ⟨r,ψ̂⟩ — uma soma ponderada de centenas de pixels
//! de ruído independentes. Pelo Teorema Central do Limite isso é gaussiano
//! SEJA QUAL FOR a distribuição do ruído, desde que a variância seja a mesma.
//!
//! PREDIÇÃO (universalidade): c é o MESMO para Gauss, Laplace e Uniforme com
//! variância σ². A forma da cauda não importa.
//!
//! PREDIÇÃO DE QUEBRA: sal-e-pimenta concentra toda a energia em pouquíssimos
//! pixels. Se o suporte do átomo pegar poucos impulsos, o CLT não tem amostras
//! e a universalidade DEVE falhar. Testamos com átomos de suporte pequeno.

use agua_atoms::{Atom, AtomBounds, Gabor2D};
use agua_core::{Field, Rng};
use agua_life::{LifeConfig, World};

/// Todos com média 0 e variância σ² — só a FORMA muda.
fn ruido(dims: &[usize], sigma: f32, seed: u64, tipo: &str) -> Field {
    let mut rng = Rng::new(seed);
    let n: usize = dims.iter().product();
    let data: Vec<f32> = (0..n)
        .map(|_| match tipo {
            "gauss" => rng.normal() * sigma,
            // Laplace(b): var = 2b² ⇒ b = σ/√2. Cauda exponencial (pesada).
            "laplace" => {
                let u = rng.f32() - 0.5;
                let b = sigma / std::f32::consts::SQRT_2;
                -b * u.signum() * (1.0 - 2.0 * u.abs()).max(1e-9).ln()
            }
            // Uniforme(-a,a): var = a²/3 ⇒ a = σ√3. Cauda LEVE (limitada).
            "uniforme" => {
                let a = sigma * 3.0f32.sqrt();
                rng.range(-a, a)
            }
            // Sal-e-pimenta: p pixels com ±A, resto 0. var = p·A².
            // Escolhemos p=1% ⇒ A = σ/√0.01 = 10σ. Energia ultra concentrada.
            t if t.starts_with("sal") => {
                let p: f32 = t[4..].parse().unwrap_or(0.01);
                if rng.f32() < p {
                    let a = sigma / p.sqrt();
                    if rng.f32() < 0.5 { a } else { -a }
                } else {
                    0.0
                }
            }
            _ => unreachable!(),
        })
        .collect();
    Field { dims: dims.to_vec(), data }
}

fn cfg(comp: usize) -> LifeConfig {
    let mut c = LifeConfig::default();
    c.cap = 400;
    c.init = 200;
    c.rain_every = 1;
    c.rain_count = 4;
    c.rain_competition = comp;
    c.tune_every = 1;
    c.gradient = true;
    c.metabolism = 1e-6;
    c.initial_energy = 50.0;
    c.repro_threshold = 1e9;
    c
}

fn c_of(dims: &[usize], sigma: f32, b: &AtomBounds, seed: u64, tipo: &str, ticks: u64) -> f32 {
    let campo = ruido(dims, sigma, seed, tipo);
    // variância REAL do campo (não a nominal) — é o que a lei usa
    let var = campo.data.iter().map(|v| v * v).sum::<f32>() / campo.len() as f32;
    let mut w: World<Gabor2D> = World::new(campo, cfg(8), b.clone(), seed);
    for _ in 0..ticks {
        w.step();
    }
    let e = w.explained_now();
    if e.is_empty() {
        return f32::NAN;
    }
    (e.iter().sum::<f32>() / e.len() as f32) / var
}

pub fn run(seed: u64) {
    let dims = [64usize, 64];
    // suporte GRANDE: centenas de pixels ⇒ o CLT deve funcionar
    let grande = AtomBounds { freq_min: 0.015, freq_max: 0.20, sigma_min: 4.0, sigma_max: 24.0 };
    // suporte PEQUENO: poucos pixels ⇒ o CLT não tem amostras
    let pequeno = AtomBounds { freq_min: 0.20, freq_max: 0.45, sigma_min: 1.0, sigma_max: 1.6 };

    println!("=== UNIVERSALIDADE: c depende da FORMA do ruído? ===");
    println!("Predição (CLT): NÃO, desde que a variância seja igual e o suporte grande.");
    println!("Predição de QUEBRA: sal-e-pimenta + suporte pequeno ⇒ c dispara.\n");

    for (nome, bounds, suporte) in [
        ("SUPORTE GRANDE (σ_atom 4–24 px, ~centenas de pixels)", &grande, "grande"),
        ("SUPORTE PEQUENO (σ_atom 1.0–1.6 px, ~dezenas de pixels)", &pequeno, "pequeno"),
    ] {
        println!("--- {nome} ---");
        println!("{:>14} {:>10} {:>10} {:>10} {:>12}", "ruído", "σ=0.05", "σ=0.10", "σ=0.20", "CV entre σ");
        let mut base = 0.0f32;
        for tipo in ["gauss", "laplace", "uniforme", "sal_pimenta"] {
            let cs: Vec<f32> = [0.05f32, 0.10, 0.20]
                .iter()
                .map(|&s| c_of(&dims, s, bounds, seed, tipo, 200))
                .collect();
            let m = cs.iter().sum::<f32>() / 3.0;
            let sd = (cs.iter().map(|c| (c - m).powi(2)).sum::<f32>() / 3.0).sqrt();
            if tipo == "gauss" {
                base = m;
            }
            let rel = 100.0 * (m - base) / base;
            println!(
                "{tipo:>14} {:>10.2} {:>10.2} {:>10.2} {:>11.1}%   média {m:.2}  ({rel:+.1}% vs gauss)",
                cs[0], cs[1], cs[2], 100.0 * sd / m
            );
        }
        let _ = suporte;
        println!();
    }
    println!("Leitura: se as 3 primeiras linhas empatam, o CLT domina e a lei é UNIVERSAL.");
    println!("Se sal-e-pimenta dispara (sobretudo no suporte pequeno), a quebra prevista ocorreu.\n");

    // ── CONVERGÊNCIA: o parâmetro é k = p·|suporte| (impulsos por átomo)? ──
    // Se sim, aumentar a densidade p faz o sal-e-pimenta convergir pro Gauss.
    println!("=== CONVERGÊNCIA CLT: c(sal-e-pimenta) → c(gauss) quando k = p·|S| cresce? ===");
    let c_gauss_g = c_of(&dims, 0.1, &grande, seed, "gauss", 200);
    let c_gauss_p = c_of(&dims, 0.1, &pequeno, seed, "gauss", 200);
    println!("referência gauss: suporte grande {c_gauss_g:.2}  |  suporte pequeno {c_gauss_p:.2}\n");
    println!("{:>10} {:>16} {:>18} {:>16} {:>18}", "p", "c (sup. grande)", "desvio vs gauss", "c (sup. pequeno)", "desvio vs gauss");
    for p in [0.01f32, 0.03, 0.10, 0.30, 0.60, 1.00] {
        let tipo = format!("sal_{p}");
        let cg = c_of(&dims, 0.1, &grande, seed, &tipo, 200);
        let cp = c_of(&dims, 0.1, &pequeno, seed, &tipo, 200);
        println!(
            "{p:>10.2} {cg:>16.2} {:>17.1}% {cp:>16.2} {:>17.1}%",
            100.0 * (cg - c_gauss_g) / c_gauss_g,
            100.0 * (cp - c_gauss_p) / c_gauss_p
        );
    }
    println!("\np=1.0 é ruído binário denso (±σ) — ainda não-gaussiano, mas denso.");
    println!("Se o desvio → 0 com p, o parâmetro de controle é k = p·|S|, não a cauda.");
}
