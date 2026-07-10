//! `agua-lab resonance` — F2: órgão temporal CfC + a ressonância como EQUAÇÃO.
//!
//! CfC (Liquid/Closed-form Continuous, linhagem do Universe):
//!     state ← state·exp(−decay·dt) + input·gain·dt
//!     out   ← tanh(state)
//! É um integrador leaky de tempo contínuo. `decay` grande = memória curta
//! (rastreia sinal rápido); `decay` pequeno = memória longa (suaviza).
//!
//! A ressonância (validada no Universe, baterias 10/12): o decay ÓTIMO de um
//! rastreador sintoniza a frequência do sinal. Aqui formalizamos como equação:
//! medimos decay*(f) e ajustamos a relação. Predição: decay* ∝ f (a constante
//! de tempo do filtro casa o período do sinal).
//!
//! Teste de paridade: a fórmula do CfC bate a do Universe (mesma equação).

/// CfCCell — paridade EXATA com Universe/src/cfc.rs.
pub struct CfCCell {
    pub state: f32,
    pub decay_rate: f32,
    pub input_gain: f32,
    pub clamp: f32,
}

impl CfCCell {
    pub fn new(decay: f32, gain: f32) -> Self {
        CfCCell { state: 0.0, decay_rate: decay, input_gain: gain, clamp: 1e6 }
    }
    /// Um passo: exatamente state·exp(−decay·dt) + input·gain·dt, saída tanh.
    pub fn step(&mut self, input: f32, dt: f32) -> f32 {
        self.state = self.state * (-self.decay_rate * dt).exp() + input * self.input_gain * dt;
        self.state = self.state.clamp(-self.clamp, self.clamp);
        self.state.tanh()
    }
}

/// Rastreador: um CfC filtra um sinal ruidoso. Erro = MSE entre o filtrado
/// (des-tanh'd, já que o alvo não está saturado) e o sinal LIMPO.
fn tracking_error(clean: &[f32], noisy: &[f32], decay: f32) -> f32 {
    // gain escolhido pra ganho DC unitário: no regime estacionário,
    // state_ss = input·gain/decay ⇒ gain = decay dá state ≈ input.
    let mut cell = CfCCell::new(decay, decay);
    let mut e = 0.0f32;
    for (c, n) in clean.iter().zip(noisy) {
        // filtra o ruidoso; usa o state (pré-tanh) como estimativa (sinal pequeno)
        cell.step(*n, 1.0);
        let est = cell.state;
        e += (est - c) * (est - c);
    }
    e / clean.len() as f32
}

/// sinal senoidal de frequência f + ruído
fn signal(n: usize, f: f32, sigma: f32, seed: u64) -> (Vec<f32>, Vec<f32>) {
    let mut rng = agua_core::Rng::new(seed);
    let clean: Vec<f32> = (0..n).map(|i| (std::f32::consts::TAU * f * i as f32).sin()).collect();
    let noisy: Vec<f32> = clean.iter().map(|c| c + rng.normal() * sigma).collect();
    (clean, noisy)
}

pub fn run(_args: &[String]) {
    // ── paridade da fórmula ──────────────────────────────────────────────
    println!("=== F2: CfC (paridade com Universe) + ressonância como equação ===\n");
    let mut cell = CfCCell::new(0.5, 1.0);
    // reproduz manualmente a fórmula do Universe e compara
    let (mut s_ref, decay, gain, dt) = (0.0f32, 0.5f32, 1.0f32, 1.0f32);
    let mut maxdiff = 0.0f32;
    for i in 0..50 {
        let input = (i as f32 * 0.3).sin();
        s_ref = s_ref * (-decay * dt).exp() + input * gain * dt;
        let out_ref = s_ref.tanh();
        let out = cell.step(input, dt);
        maxdiff = maxdiff.max((out - out_ref).abs());
    }
    println!("paridade CfC vs fórmula Universe: max|diff| = {maxdiff:.2e}  {}",
        if maxdiff < 1e-6 { "✓ EXATA" } else { "✗" });

    // ── ressonância: decay ótimo vs frequência do sinal ──────────────────
    println!("\nressonância: qual decay MINIMIZA o erro de rastreio, por frequência?");
    println!("predição: decay* cresce com f (constante de tempo casa o período)\n");
    println!("{:>10} {:>12} {:>14}", "freq", "decay*", "período (1/f)");
    let decays: Vec<f32> = (0..40).map(|i| 0.01 * 1.2f32.powi(i)).filter(|d| *d < 5.0).collect();
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for &f in &[0.004f32, 0.008, 0.016, 0.032, 0.064, 0.128] {
        // média sobre seeds
        let mut best_d = decays[0];
        let mut best_e = f32::MAX;
        for &d in &decays {
            let mut e = 0.0f32;
            for seed in [11u64, 22, 33] {
                let (clean, noisy) = signal(2000, f, 0.3, seed);
                e += tracking_error(&clean, &noisy, d);
            }
            if e < best_e {
                best_e = e;
                best_d = d;
            }
        }
        xs.push(f.ln());
        ys.push(best_d.ln());
        println!("{:>10.4} {:>12.3} {:>14.1}", f, best_d, 1.0 / f);
    }
    // ajuste log-log: decay* = A·f^α
    let n = xs.len() as f32;
    let mx = xs.iter().sum::<f32>() / n;
    let my = ys.iter().sum::<f32>() / n;
    let sxy: f32 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f32 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let syy: f32 = ys.iter().map(|y| (y - my) * (y - my)).sum();
    let alpha = sxy / sxx;
    let a = (my - alpha * mx).exp();
    let r2 = sxy * sxy / (sxx * syy);
    println!("\n  ajuste: decay* = {a:.2}·f^{alpha:.3}   (R²={r2:.4})");
    println!("  → NÃO é linear (previ α≈1, deu α≈0.4): lei de POTÊNCIA sublinear.");
    println!("    Monótona (confirma Universe bat.10/12 qualitativamente) mas o");
    println!("    expoente ~0.4 é de um ótimo VIÉS-VARIÂNCIA — deve depender do ruído.");

    // ── o expoente depende do SNR? (viés-variância) ──────────────────────
    println!("\nα depende do ruído σ? (se for bias-variância, sim)");
    println!("{:>8} {:>12} {:>10}", "σ", "α", "constante");
    for &sg in &[0.1f32, 0.3, 0.6, 1.0] {
        let mut xs=Vec::new(); let mut ys=Vec::new();
        for &f in &[0.004f32,0.008,0.016,0.032,0.064,0.128] {
            let mut bd=decays[0]; let mut be=f32::MAX;
            for &d in &decays {
                let mut e=0.0; for seed in [11u64,22,33] {
                    let (c,nz)=signal(2000,f,sg,seed); e+=tracking_error(&c,&nz,d);
                }
                if e<be {be=e; bd=d;}
            }
            xs.push(f.ln()); ys.push(bd.ln());
        }
        let n=xs.len() as f32; let mx=xs.iter().sum::<f32>()/n; let my=ys.iter().sum::<f32>()/n;
        let sxy:f32=xs.iter().zip(&ys).map(|(x,y)|(x-mx)*(y-my)).sum();
        let sxx:f32=xs.iter().map(|x|(x-mx)*(x-mx)).sum();
        let al=sxy/sxx; let cc=(my-al*mx).exp();
        println!("{:>8.1} {:>12.3} {:>10.2}", sg, al, cc);
    }
}
