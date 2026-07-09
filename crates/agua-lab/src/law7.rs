//! `agua-lab law7` — dof EFETIVO: a lei usa o posto do Jacobiano, não o nº de botões.
//!
//! Lei 2 mediu c = 2·dof + 1 com R² 0.99 no Gabor2D (dof=6 ⇒ c=13.19, previsto 13).
//! Mas no Gabor1D (dof=4) deu 7.50, não 9 — 17% abaixo, sistematicamente.
//!
//! HIPÓTESE: os parâmetros não são independentes. Em 1D, mover o centro c é quase
//! equivalente a mudar a fase φ (a portadora é cos(k·dx+φ)). O ser tem 4 botões
//! mas menos direções REAIS de busca. O que conta é o posto do Jacobiano ∂ψ/∂p.
//!
//! MEDIDA INDEPENDENTE: monta a matriz de Gram G_ij = ⟨∂ψ/∂pᵢ, ∂ψ/∂pⱼ⟩ (por
//! diferenças finitas, sem tocar no código do gradiente) e calcula a RAZÃO DE
//! PARTICIPAÇÃO dos autovalores:
//!
//!     dof_eff = (Σλ)² / Σλ²
//!
//! (= nº de direções que realmente contribuem; vale d se todas iguais, 1 se uma domina)
//!
//! PREDIÇÃO: dof_eff medido pela Gram ≈ (c − 1)/2 medido pela sobrevivência.
//!   Gabor2D: (13.19−1)/2 = 6.10       Gabor1D: (7.50−1)/2 = 3.25
//! Se baterem, a lei fecha na forma  met_c = σ²·(1 + 2·dof_eff).

use agua_atoms::{Atom, AtomBounds, Gabor1D, Gabor2D};
use agua_core::Rng;

/// ∂ψ/∂pₖ por diferença central, avaliado num campo denso do tamanho de `dims`.
fn dpsi<A: Atom>(a: &A, dims: &[usize], b: &AtomBounds, k: usize, eps: f32) -> Vec<f32> {
    let n: usize = dims.iter().product();
    let mut dense = |atom: &A| -> Vec<f32> {
        let p = atom.render_unit(dims);
        let mut v = vec![0.0f32; n];
        let mut off = 0usize;
        for run in &p.runs {
            v[run.field_start..run.field_start + run.len]
                .copy_from_slice(&p.data[off..off + run.len]);
            off += run.len;
        }
        v
    };
    let mut d = vec![0.0f32; A::GEOM];
    d[k] = eps;
    let mut ap = a.clone();
    ap.apply_delta(&d, b, dims);
    d[k] = -eps;
    let mut am = a.clone();
    am.apply_delta(&d, b, dims);
    let (vp, vm) = (dense(&ap), dense(&am));
    vp.iter().zip(&vm).map(|(x, y)| (x - y) / (2.0 * eps)).collect()
}

/// Autovalores de uma matriz simétrica pequena via Jacobi cíclico.
fn eigvals(mut a: Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    for _ in 0..80 {
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += a[i][j] * a[i][j];
            }
        }
        if off < 1e-18 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() < 1e-15 {
                    continue;
                }
                let theta = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }
    (0..n).map(|i| a[i][i]).collect()
}

/// dof efetivo = (Σλ)²/Σλ² da Gram dos gradientes, com cada ∂ψ/∂pₖ NORMALIZADO
/// (senão a escala arbitrária de cada parâmetro domina a conta).
fn dof_efetivo<A: Atom>(dims: &[usize], b: &AtomBounds, seed: u64, amostras: usize, eps: &[f32]) -> f32 {
    let mut rng = Rng::new(seed);
    let mut soma = 0.0f64;
    let mut n_ok = 0usize;
    for _ in 0..amostras {
        let pos: Vec<f32> = dims.iter().map(|d| rng.range(0.2 * *d as f32, 0.8 * *d as f32)).collect();
        let atom = A::rain_at(&pos, b, &mut rng);
        let gs: Vec<Vec<f32>> = (0..A::GEOM).map(|k| dpsi(&atom, dims, b, k, eps[k])).collect();
        // normaliza cada direção
        let gs: Vec<Vec<f32>> = gs
            .into_iter()
            .map(|g| {
                let nrm = g.iter().map(|v| v * v).sum::<f32>().sqrt();
                if nrm < 1e-9 {
                    g
                } else {
                    g.iter().map(|v| v / nrm).collect()
                }
            })
            .collect();
        if gs.iter().any(|g| g.iter().all(|v| v.abs() < 1e-12)) {
            continue;
        }
        let d = A::GEOM;
        let mut gram = vec![vec![0.0f64; d]; d];
        for i in 0..d {
            for j in 0..d {
                gram[i][j] = gs[i].iter().zip(&gs[j]).map(|(x, y)| (*x as f64) * (*y as f64)).sum();
            }
        }
        let ev = eigvals(gram);
        let s1: f64 = ev.iter().map(|v| v.max(0.0)).sum();
        let s2: f64 = ev.iter().map(|v| v.max(0.0) * v.max(0.0)).sum();
        if s2 > 1e-12 {
            soma += s1 * s1 / s2;
            n_ok += 1;
        }
    }
    (soma / n_ok.max(1) as f64) as f32
}

pub fn run(seed: u64) {
    println!("=== dof EFETIVO: a lei conta DIREÇÕES, não botões ===");
    println!("dof_eff = (Σλ)²/Σλ² da Gram dos ∂ψ/∂p (razão de participação).");
    println!("PREDIÇÃO: dof_eff ≈ (c − 1)/2, com c medido pela SOBREVIVÊNCIA.\n");

    let b2 = AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
    let b1 = AtomBounds { freq_min: 0.002, freq_max: 0.45, sigma_min: 2.0, sigma_max: 400.0 };

    // eps por parâmetro (mesma ordem de apply_delta)
    let e2 = [0.05f32, 0.05, 5e-4, 5e-3, 5e-3, 0.05];
    let e1 = [0.05f32, 5e-4, 5e-3, 0.05];

    let d2 = dof_efetivo::<Gabor2D>(&[64, 64], &b2, seed, 60, &e2);
    let d1 = dof_efetivo::<Gabor1D>(&[4096], &b1, seed, 60, &e1);

    println!("{:>16} {:>10} {:>12} {:>16} {:>10}", "átomo", "dof nominal", "dof_eff", "(c−1)/2 medido", "erro");
    println!(
        "{:>16} {:>11} {:>12.2} {:>16.2} {:>9.1}%",
        "Gabor2D", 6, d2, 6.10, 100.0 * (d2 - 6.10) / 6.10
    );
    println!(
        "{:>16} {:>11} {:>12.2} {:>16.2} {:>9.1}%",
        "Gabor1D", 4, d1, 3.25, 100.0 * (d1 - 3.25) / 3.25
    );

    println!("\n(c medido na fase harmônica: Gabor2D 13.19, Gabor1D 7.50)");
    println!("Se dof_eff ≈ (c−1)/2, a lei fecha:  met_c = σ²·(1 + 2·dof_eff)");
    println!("\nPredição de c a partir de dof_eff:");
    println!("  Gabor2D: c = 1 + 2·{d2:.2} = {:.2}   (medido 13.19)", 1.0 + 2.0 * d2);
    println!("  Gabor1D: c = 1 + 2·{d1:.2} = {:.2}   (medido 7.50)", 1.0 + 2.0 * d1);

    println!("\n=== VEREDITO: hipótese da razão de participação REFUTADA ===");
    println!("Não vou trocar de métrica até alguma bater — isso seria p-hacking.");
    println!("Registro o diagnóstico bruto da Gram, sem alegar que explica c:\n");

    // diagnóstico: espectro da Gram (autovalores normalizados) — só descritivo
    for (nome, ev) in [
        ("Gabor2D", espectro::<Gabor2D>(&[64, 64], &b2, seed, 40, &e2)),
        ("Gabor1D", espectro::<Gabor1D>(&[4096], &b1, seed, 40, &e1)),
    ] {
        let lmax = ev[0].max(1e-12);
        let posto_1pct = ev.iter().filter(|l| **l > 0.01 * lmax).count();
        let posto_01pct = ev.iter().filter(|l| **l > 0.001 * lmax).count();
        let cond = lmax / ev.last().copied().unwrap_or(1e-12).max(1e-12);
        let espec: Vec<String> = ev.iter().map(|l| format!("{:.3}", l / lmax)).collect();
        println!("{nome}: espectro λ/λmax = [{}]", espec.join(", "));
        println!("   posto numérico: {posto_1pct} (tol 1%), {posto_01pct} (tol 0.1%); condicionamento {cond:.1e}");
    }
    println!("\nO 1D tem direções mais degeneradas que o 2D (mover c ≈ mudar φ),");
    println!("o que é CONSISTENTE com c_1D < c_2D por grau — mas a razão de");
    println!("participação não prevê o valor. Fica em aberto.");
}

/// autovalores médios (ordenados desc) da Gram das direções normalizadas
fn espectro<A: Atom>(dims: &[usize], b: &AtomBounds, seed: u64, amostras: usize, eps: &[f32]) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    let mut acc = vec![0.0f64; A::GEOM];
    let mut n_ok = 0usize;
    for _ in 0..amostras {
        let pos: Vec<f32> = dims.iter().map(|d| rng.range(0.2 * *d as f32, 0.8 * *d as f32)).collect();
        let atom = A::rain_at(&pos, b, &mut rng);
        let gs: Vec<Vec<f32>> = (0..A::GEOM)
            .map(|k| {
                let g = dpsi(&atom, dims, b, k, eps[k]);
                let n = g.iter().map(|v| v * v).sum::<f32>().sqrt();
                if n < 1e-9 { g } else { g.iter().map(|v| v / n).collect() }
            })
            .collect();
        let d = A::GEOM;
        let mut gram = vec![vec![0.0f64; d]; d];
        for i in 0..d {
            for j in 0..d {
                gram[i][j] = gs[i].iter().zip(&gs[j]).map(|(x, y)| (*x as f64) * (*y as f64)).sum();
            }
        }
        let mut ev = eigvals(gram);
        ev.sort_by(|a, b| b.partial_cmp(a).unwrap());
        if ev[0].is_finite() {
            for (a, e) in acc.iter_mut().zip(&ev) {
                *a += e.max(0.0);
            }
            n_ok += 1;
        }
    }
    acc.iter().map(|v| v / n_ok.max(1) as f64).collect()
}
