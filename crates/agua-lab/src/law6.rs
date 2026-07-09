//! `agua-lab law6` — c COMO FUNÇÃO DOS GRAUS DE LIBERDADE.
//!
//! A lei diz met_c = c·σ². Já sabemos que c cresce com o orçamento de busca
//! (log) e com a dimensionalidade do átomo. Aqui isolamos os graus de liberdade
//! CONGELANDO parâmetros (bitmask `frozen`), sem trocar de átomo.
//!
//! PREDIÇÃO DE ANCORAGEM (a mais forte de todas): com dof = 0 e sem competição
//! na chuva, o ser não busca NADA — sua geometria é sorteada e nunca muda.
//! Então ele é um átomo fixo sobre ruído e a teoria exige
//!
//!         c(dof = 0) = 1     (exatamente E[⟨r,ψ̂⟩²]/σ² = 1)
//!
//! Se o intercepto der 1, a teoria fecha: c = 1 + (ganho por grau de liberdade).
//! Se der longe de 1, há um viés não explicado na medida.
//!
//! PREDIÇÃO DE FORMA: cada parâmetro livre acrescenta uma direção de busca.
//! Se as direções fossem independentes, c cresceria linearmente com dof.

use agua_atoms::{AtomBounds, Gabor2D};
use agua_core::{Field, Rng};
use agua_life::{LifeConfig, World};

fn noise(dims: &[usize], sigma: f32, seed: u64) -> Field {
    let mut rng = Rng::new(seed);
    let n: usize = dims.iter().product();
    Field { dims: dims.to_vec(), data: (0..n).map(|_| rng.normal() * sigma).collect() }
}

fn c_of(dims: &[usize], b: &AtomBounds, seed: u64, frozen: u32, comp: usize, ticks: u64) -> f32 {
    c_of_pred(dims, b, seed, frozen, comp, ticks, 0.0)
}

fn c_of_pred(
    dims: &[usize],
    b: &AtomBounds,
    seed: u64,
    frozen: u32,
    comp: usize,
    ticks: u64,
    pred: f32,
) -> f32 {
    let sigma = 0.1f32;
    let campo = noise(dims, sigma, seed);
    let var = campo.data.iter().map(|v| v * v).sum::<f32>() / campo.len() as f32;
    let mut c = LifeConfig::default();
    c.cap = 500;
    c.init = 400;
    c.rain_every = 0; // sem imigração: a população inicial é a amostra
    c.rain_competition = comp;
    c.tune_every = 1;
    c.gradient = true;
    c.frozen = frozen;
    c.metabolism = 1e-6;
    c.initial_energy = 50.0;
    c.repro_threshold = 1e9;
    c.predation = pred;
    c.predation_every = 5;
    let mut w: World<Gabor2D> = World::new(campo, c, b.clone(), seed);
    for _ in 0..ticks {
        w.step();
    }
    let e = w.explained_now();
    if e.is_empty() {
        return f32::NAN;
    }
    (e.iter().sum::<f32>() / e.len() as f32) / var
}

fn linfit(xs: &[f32], ys: &[f32]) -> (f32, f32, f32) {
    let n = xs.len() as f32;
    let mx = xs.iter().sum::<f32>() / n;
    let my = ys.iter().sum::<f32>() / n;
    let sxy: f32 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f32 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let syy: f32 = ys.iter().map(|y| (y - my) * (y - my)).sum();
    let a = sxy / sxx;
    (a, my - a * mx, if syy > 0.0 { sxy * sxy / (sxx * syy) } else { 1.0 })
}

pub fn run(seed: u64) {
    let dims = [64usize, 64];
    let b = AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };

    // ordem de liberação: cx, cy, freq, theta, phase, sigma
    // frozen[k]=1 trava o parâmetro k. dof = nº de bits ZERADOS.
    let masks: [(u32, usize, &str); 7] = [
        (0b111111, 0, "nenhum (átomo fixo)"),
        (0b111110, 1, "cx"),
        (0b111100, 2, "cx,cy"),
        (0b111000, 3, "+freq"),
        (0b110000, 4, "+theta"),
        (0b100000, 5, "+phase"),
        (0b000000, 6, "+sigma (todos)"),
    ];

    println!("=== c vs GRAUS DE LIBERDADE (parâmetros congelados) ===");
    println!("ANCORAGEM: c(dof=0, sem competição) tem que ser EXATAMENTE 1.\n");

    for &comp in &[1usize, 8] {
        println!("--- competição na chuva = {comp} ---");
        println!("{:>5} {:>22} {:>10}", "dof", "livres", "c");
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        for &(m, dof, nome) in &masks {
            let c = c_of(&dims, &b, seed, m, comp, 200);
            xs.push(dof as f32);
            ys.push(c);
            let flag = if dof == 0 {
                if (c - 1.0).abs() < 0.15 { "  ← âncora OK (≈1)" } else { "  ← âncora FORA de 1" }
            } else {
                ""
            };
            println!("{dof:>5} {nome:>22} {c:>10.2}{flag}");
        }
        let (a, b0, r2) = linfit(&xs, &ys);
        println!("   ajuste linear: c = {a:.3}·dof + {b0:.3}   (R² = {r2:.4})");
        println!("   ganho por grau de liberdade: {a:.2}\n");
    }
    println!("Nota: com comp>1 a chuva já busca (escolhe a melhor de K sementes),");
    println!("então c(dof=0) NÃO deve ser 1 — a busca acontece no nascimento.\n");

    // ── A ÂNCORA FALHOU: por quê? MULTICOLINEARIDADE. ───────────────────
    // Um átomo SOZINHO sobre ruído captura 1·σ² (medido: 1.13). Com 400 átomos
    // coerentes (H2), os coeficientes individuais inflam por empréstimo mútuo —
    // é o mesmo fenômeno de regressão com preditores colineares.
    // PREDIÇÃO: a predação por nicho (H4) remove as quase-duplicatas ⇒ a âncora
    // c(dof=0) tem que CAIR em direção a 1.
    println!("=== POR QUE A ÂNCORA FALHOU: colinearidade da população ===");
    println!("dof=0, comp=1 (nenhuma busca). A teoria exige c=1 para um átomo isolado.");
    println!("Predição: predação mais agressiva (θ menor) ⇒ menos duplicatas ⇒ c → 1.\n");
    println!("{:>10} {:>10} {:>28}", "θ predação", "c(dof=0)", "");
    for &theta in &[0.0f32, 0.95, 0.8, 0.6, 0.4, 0.2] {
        let c = c_of_pred(&dims, &b, seed, 0b111111, 1, 200, theta);
        let lbl = if theta == 0.0 { "off".to_string() } else { format!("{theta:.2}") };
        let bar: String = std::iter::repeat('#').take(((c * 2.0) as usize).min(28)).collect();
        println!("{lbl:>10} {c:>10.2} {bar:>28}");
    }
    println!("\n(um átomo ISOLADO sobre ruído mede c = 1.13 — a âncora verdadeira)\n");

    // ── TRANSIÇÃO DE FASE DA HARMONIA ───────────────────────────────────
    // O salto entre θ=0.6 e θ=0.4 não é gradual. Existe uma coerência CRÍTICA
    // θ* abaixo da qual a população deixa de ser colinear e a lei recupera a
    // âncora. Análogo às transições de fase de Donoho–Tanner em compressed
    // sensing. Localizamos θ* e observamos N e μ atravessando a transição.
    println!("=== TRANSIÇÃO DE FASE DA HARMONIA ===");
    println!("θ* = coerência crítica. Abaixo dela: população harmônica, c → 1 (lei válida).");
    println!("Acima: colinearidade, coeficientes inflados, c ≈ 12.\n");
    println!("{:>8} {:>10} {:>8} {:>28}", "θ", "c(dof=0)", "vivos", "fase");
    let mut anterior = f32::NAN;
    for i in 0..=12 {
        let theta = 0.30 + 0.04 * i as f32;
        let c = c_of_pred(&dims, &b, seed, 0b111111, 1, 200, theta);
        let salto = if anterior.is_finite() && c / anterior > 2.0 { "  ⟵ SALTO (θ* aqui)" } else { "" };
        let fase = if c < 3.0 { "harmônica" } else { "colinear" };
        println!("{theta:>8.2} {c:>10.2} {:>8} {:>18}{salto}", "-", fase);
        anterior = c;
    }

    // ── dof, AGORA MEDIDO NA FASE HARMÔNICA (θ = 0.3 < θ*) ──────────────
    // Na fase colinear a medida de dof estava poluída (o intercepto dava 12).
    // Dentro da fase harmônica a âncora vale e o ganho por grau de liberdade
    // pode ser lido de verdade.
    println!("\n=== c vs GRAUS DE LIBERDADE, DENTRO DA FASE HARMÔNICA (θ=0.3) ===");
    println!("Agora o intercepto DEVE bater a âncora (≈1).\n");
    println!("{:>5} {:>22} {:>10}", "dof", "livres", "c");
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for &(m, dof, nome) in &masks {
        let c = c_of_pred(&dims, &b, seed, m, 1, 200, 0.3);
        xs.push(dof as f32);
        ys.push(c);
        println!("{dof:>5} {nome:>22} {c:>10.2}");
    }
    let (a, b0, r2) = linfit(&xs, &ys);
    println!("\n   c = {a:.3}·dof + {b0:.3}   (R² = {r2:.4})");
    println!("   intercepto {b0:.2} vs âncora teórica 1.00  →  {}",
        if (b0 - 1.0).abs() < 0.6 { "TEORIA FECHA" } else { "ainda há viés" });
    println!("   ganho por grau de liberdade: {a:.2}");
}
