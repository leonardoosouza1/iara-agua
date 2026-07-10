//! `agua-lab synapse` — as leis da SINAPSE: acoplamento, alcance, avalanche.
//!
//! Pergunta do Leonardo: existe lei de sinapse? disparo? raio de distância? como
//! um neurônio que dispara atinge outro? por que o outro dispara?
//!
//! A resposta já está no modelo, implícita. A reconstrução ótima resolve as
//! equações normais G·a = b, onde G_ij = ⟨ψ_i, ψ_j⟩ é a matriz de Gram.
//! **G É a matriz de pesos sinápticos.** Quando o átomo j "dispara" (muda amp em
//! Δ), o resíduo que i sente muda em −Δ·G_ij; a amplitude ótima de i se desloca
//! em −Δ·G_ij/G_ii. Essa é a TRANSMISSÃO SINÁPTICA. G_ij decai com a distância
//! porque o overlap de gaussianas decai — é o RAIO DE ALCANCE.
//!
//! LEI 1 (alcance): G_ij(d) decai gaussiano com a distância d entre centros;
//!   comprimento característico ~ σ (o tamanho do átomo). Medimos e ajustamos.
//!
//! LEI 2 (avalanche/criticalidade): perturbar um átomo (um disparo) propaga pela
//!   rede de overlap. Predição forte, ligada à Lei 6 (transição de fase da
//!   harmonia): no regime HARMÔNICO (baixa coerência) o disparo DAMPA rápido; no
//!   regime COLINEAR (alta coerência) ele PROPAGA longe. O alcance da avalanche
//!   DIVERGE perto do ponto crítico θ* — assinatura de criticalidade, como as
//!   avalanches neurais reais (Beggs & Plenz 2003, self-organized criticality).

use agua_atoms::{Atom, AtomBounds, Gabor2D};
use agua_core::{Field, Rng};

fn bounds() -> AtomBounds {
    AtomBounds { freq_min: 0.02, freq_max: 0.3, sigma_min: 3.0, sigma_max: 8.0 }
}

/// LEI 1 — acoplamento sináptico vs distância.
fn coupling_vs_distance(seed: u64) {
    println!("=== LEI DA SINAPSE 1 — acoplamento G_ij vs distância (o raio de alcance) ===");
    println!("G_ij = ⟨ψ̂_i, ψ̂_j⟩ entre átomos IDÊNTICOS (mesma freq/θ/σ) separados por d.\n");
    let dims = [64usize, 64];
    let b = bounds();
    let mut rng = Rng::new(seed);
    // átomo de referência no centro
    let sigma = 5.0f32;
    let base = Gabor2D { cx: 32.0, cy: 32.0, freq: 0.12, theta: 0.6, phase: 0.3, sigma };
    let pu = base.render_unit(&dims);
    let nb = pu.norm2().sqrt();
    println!("desloco AO LONGO de θ (a portadora faz o sinal oscilar = chapéu-mexicano):");
    println!("{:>8} {:>14} {:>12}", "d (px)", "G_ij (SINAL)", "sinal");
    // envelope: desloco PERPENDICULAR a θ (só o envelope gaussiano, sem oscilar)
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    let (st, ct) = base.theta.sin_cos();
    for d in [0.0f32, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0] {
        // ao longo de θ (oscila)
        let along = Gabor2D { cx: 32.0 + d*ct, cy: 32.0 + d*st, ..base.clone() };
        let pa = along.render_unit(&dims); let na = pa.norm2().sqrt();
        let g_along = if nb*na>1e-9 { pu.dot_patch(&pa)/(nb*na) } else {0.0};
        // perpendicular a θ (só envelope)
        let perp = Gabor2D { cx: 32.0 - d*st, cy: 32.0 + d*ct, ..base.clone() };
        let pp = perp.render_unit(&dims); let np_ = pp.norm2().sqrt();
        let g_perp = if nb*np_>1e-9 { pu.dot_patch(&pp).abs()/(nb*np_) } else {0.0};
        let sinal = if g_along>0.05 {"EXCITA (+)"} else if g_along < -0.05 {"INIBE (−)"} else {"~0"};
        println!("{:>8.1} {:>14.4} {:>12}", d, g_along, sinal);
        if d>0.0 && g_perp>1e-3 { xs.push(d*d); ys.push(g_perp.ln()); }
    }
    // ajuste ln G = a·d² + b ⇒ G = exp(b)·exp(a·d²), comprimento λ = 1/√(-a)
    let n = xs.len() as f32;
    let mx = xs.iter().sum::<f32>() / n;
    let my = ys.iter().sum::<f32>() / n;
    let sxy: f32 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f32 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let a = sxy / sxx;
    let lambda = if a < 0.0 { (-1.0 / a).sqrt() } else { f32::NAN };
    println!("\n  envelope (perpendicular a θ): G(d) ∝ exp(−d²/λ²), λ = {lambda:.2} px");
    println!("  → o ENVELOPE do acoplamento é gaussiano (alcance ~σ√2={:.1}); a PORTADORA", sigma*std::f32::consts::SQRT_2);
    println!("    faz o sinal ALTERNAR excita(+)/inibe(−) com a distância = chapéu-mexicano,");
    println!("    a inibição lateral do córtex visual. A sinapse = a própria forma do átomo.\n");
    let _ = (rng.f32(), &b);
}

/// mundo pequeno pra medir avalanche: ajusta uma população a um alvo, com
/// coerência controlada pela predação (ligação com a Lei 6).
fn build_population(theta_pred: f32, seed: u64) -> (Vec<Gabor2D>, Vec<f32>, Field, [usize; 2]) {
    use agua_life::{LifeConfig, World};
    let size = 48usize;
    let dims = [size, size];
    let s = size as f32;
    let mut t = Field::from_fn_2d(size, size, |x, y| {
        let (xf, yf) = (x as f32, y as f32);
        0.4 + 0.3 * (std::f32::consts::TAU * 0.1 * xf).sin() * (std::f32::consts::TAU * 0.08 * yf).cos()
            + 0.2 * (-((xf - s * 0.5).powi(2) + (yf - s * 0.5).powi(2)) / (2.0 * (s * 0.2).powi(2))).exp()
    });
    let m = t.mean();
    for v in t.data.iter_mut() {
        *v -= m;
    }
    let mut cfg = LifeConfig::default();
    cfg.cap = 300;
    cfg.metabolism = 0.05;
    cfg.rain_every = 1;
    cfg.rain_count = 4;
    cfg.rain_competition = 8;
    cfg.tune_every = 1;
    cfg.gradient = false; // mutação → população mais COERENTE (aglomera nos hotspots)
    cfg.metabolism = 0.01; // baixo → over-completa → alta coerência
    cfg.predation = theta_pred;
    cfg.predation_every = 3;
    let mut w: World<Gabor2D> = World::new(t, cfg, bounds(), seed);
    for _ in 0..400 {
        w.step();
    }
    w.rebuild();
    let atoms: Vec<Gabor2D> = w.beings.iter().map(|b| b.atom.clone()).collect();
    let amps: Vec<f32> = w.beings.iter().map(|b| b.amp).collect();
    (atoms, amps, w.residual.clone(), dims)
}

/// LEI 2 — avalanche: perturba um átomo, relaxa (Gauss-Seidel), mede quantos
/// átomos se moveram (tamanho da avalanche) e quão longe (alcance espacial).
#[allow(dead_code)]
fn avalanche(seed: u64) {
    println!("=== LEI DA SINAPSE 2 — avalanche: o disparo propaga? diverge no crítico? ===");
    println!("perturba a amp de 1 átomo (disparo), relaxa a rede, mede o tamanho da");
    println!("cascata. Predição (liga Lei 6): DAMPA no harmônico, PROPAGA no colinear.\n");
    println!("{:>10} {:>12} {:>16} {:>16}", "θ predação", "coerência μ", "tamanho aval.", "alcance (px)");

    // predação ALTA baixa μ (harmônico); SEM predação = alta μ (colinear)
    for &theta in &[0.0f32, 0.85, 0.7, 0.55, 0.45, 0.35] {
        let (atoms, amps0, _resid, dims) = build_population(theta, seed);
        let n = atoms.len();
        if n < 5 {
            continue;
        }
        // patches unitários e centros
        let units: Vec<_> = atoms.iter().map(|a| a.render_unit(&dims)).collect();
        let n2: Vec<f32> = units.iter().map(|u| u.norm2()).collect();
        let centers: Vec<(f32, f32)> = atoms.iter().map(|a| (a.cx, a.cy)).collect();
        // coerência média (só pares que se tocam)
        let mut mu = 0.0f32;
        let mut np = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let ni = n2[i].sqrt();
                let nj = n2[j].sqrt();
                if ni * nj < 1e-9 {
                    continue;
                }
                let g = units[i].dot_patch(&units[j]).abs() / (ni * nj);
                if g > 1e-4 {
                    mu += g;
                    np += 1;
                }
            }
        }
        mu /= np.max(1) as f32;

        // reconstrói o campo com as amps de equilíbrio
        let mut recon = Field::zeros(&dims);
        for (u, a) in units.iter().zip(&amps0) {
            recon.add_scaled_patch(u, *a);
        }
        // DISPARO: perturba o átomo 0 em +50% da sua amp (ou +0.2 se amp~0)
        let src = 0;
        let delta = if amps0[src].abs() > 1e-3 { 0.5 * amps0[src] } else { 0.2 };
        let mut amps = amps0.clone();
        amps[src] += delta;
        recon.add_scaled_patch(&units[src], delta);

        // RELAXA: Gauss-Seidel — cada átomo reajusta amp contra o resíduo atual.
        // (o resíduo do alvo = alvo − recon; sem alvo aqui, relaxamos a recon de
        //  volta ao equilíbrio: cada átomo minimiza o desvio local causado pelo
        //  disparo — é a resposta pós-sináptica em cadeia.)
        // resíduo = (recon perturbada) − (recon equilíbrio)  → queremos zerar
        let mut perturb = Field::zeros(&dims);
        perturb.add_scaled_patch(&units[src], delta); // o que o disparo injetou
        // relaxa: cada átomo j (≠src) absorve parte da perturbação via overlap
        let mut moved = vec![false; n];
        for _sweep in 0..4 {
            for j in 0..n {
                if j == src {
                    continue;
                }
                // amp que cancela a projeção da perturbação em ψ_j
                let d = units[j].dot_field(&perturb);
                let adj = d / n2[j].max(1e-9);
                if adj.abs() > 0.02 * delta.abs().max(1e-6) {
                    perturb.add_scaled_patch(&units[j], -adj);
                    amps[j] -= adj;
                    moved[j] = true;
                }
            }
        }
        let size_av = moved.iter().filter(|&&m| m).count();
        // alcance = distância máxima de um átomo que se moveu ao src
        let (sx, sy) = centers[src];
        let reach = (0..n)
            .filter(|&j| moved[j])
            .map(|j| {
                let (x, y) = centers[j];
                ((x - sx).powi(2) + (y - sy).powi(2)).sqrt()
            })
            .fold(0.0f32, f32::max);
        let lbl = if theta == 0.0 { "off".to_string() } else { format!("{theta:.2}") };
        let flag = if mu > 0.35 { " ◄ colinear" } else { " harmônico" };
        println!("{lbl:>10} {mu:>12.3} {size_av:>10}/{n:<4} {reach:>16.1}{flag}");
    }
    println!("\n  Se o alcance/tamanho da avalanche cresce quando μ passa de ~0.35 (θ*),");
    println!("  o disparo propaga longe no regime colinear — criticalidade, como as");
    println!("  avalanches neurais. A sinapse (overlap) é forte demais lá: um dispara,");
    println!("  muitos respondem. No harmônico, cada um cuida do seu (dampa).");
}

/// LEI 2b — CRITICALIDADE isolada da dinâmica: constrói populações de coerência
/// CONTROLADA (átomos empilhados = alta μ; espalhados = baixa μ), perturba um,
/// mede a avalanche. Testa a física do acoplamento sem a sobrevivência interferir.
fn criticality_direct(seed: u64) {
    println!("=== LEI DA SINAPSE 2b — criticalidade isolada (coerência CONSTRUÍDA) ===");
    println!("N átomos idênticos num raio R: R pequeno = empilhados (alta μ) = colinear.");
    println!("perturba 1, relaxa, mede a cascata. Predição: avalanche DIVERGE em alta μ.\n");
    let dims = [48usize, 48];
    let n = 40usize;
    println!("{:>10} {:>12} {:>16} {:>16}", "raio R", "coerência μ", "tamanho aval.", "alcance (px)");
    for &r in &[2.0f32, 4.0, 6.0, 9.0, 13.0, 18.0] {
        let mut rng = Rng::new(seed);
        // N átomos idênticos espalhados num disco de raio R em torno do centro
        let atoms: Vec<Gabor2D> = (0..n).map(|_| {
            let ang = rng.range(0.0, std::f32::consts::TAU);
            let rr = r * rng.f32().sqrt();
            Gabor2D { cx: 24.0 + rr*ang.cos(), cy: 24.0 + rr*ang.sin(),
                      freq: 0.12, theta: 0.6, phase: 0.3, sigma: 4.0 }
        }).collect();
        let units: Vec<_> = atoms.iter().map(|a| a.render_unit(&dims)).collect();
        let n2: Vec<f32> = units.iter().map(|u| u.norm2()).collect();
        let centers: Vec<(f32,f32)> = atoms.iter().map(|a|(a.cx,a.cy)).collect();
        let mut mu=0.0; let mut npair=0;
        for i in 0..n { for j in (i+1)..n {
            let ni=n2[i].sqrt(); let nj=n2[j].sqrt();
            if ni*nj<1e-9 {continue;}
            let g=units[i].dot_patch(&units[j]).abs()/(ni*nj);
            if g>1e-3 { mu+=g; npair+=1; }
        }}
        mu/= npair.max(1) as f32;
        // disparo no átomo central (0), relaxa
        let delta=0.2f32;
        let mut perturb=Field::zeros(&dims);
        perturb.add_scaled_patch(&units[0], delta);
        let mut moved=vec![false;n];
        for _ in 0..6 { for j in 1..n {
            let d=units[j].dot_field(&perturb); let adj=d/n2[j].max(1e-9);
            if adj.abs()>0.02*delta { perturb.add_scaled_patch(&units[j], -adj); moved[j]=true; }
        }}
        let sz=moved.iter().filter(|&&m|m).count();
        let (sx,sy)=centers[0];
        let reach=(0..n).filter(|&j|moved[j]).map(|j|{let(x,y)=centers[j];((x-sx).powi(2)+(y-sy).powi(2)).sqrt()}).fold(0.0f32,f32::max);
        let fase = if mu>0.35 {" ◄ COLINEAR (crítico)"} else {" harmônico"};
        println!("{r:>10.1} {mu:>12.3} {sz:>10}/{n:<4} {reach:>16.1}{fase}");
    }
    // A MÉTRICA CERTA: susceptibilidade = ‖Δamp‖ / ‖disparo‖ ao resolver as
    // equações normais G·a = b. Colinear ⇒ G quase-singular ⇒ resposta explode.
    // É a MESMA colinearidade que inflou o c na Lei 6 (uma história só).
    println!("\n  SUSCEPTIBILIDADE — quanto as amplitudes reagem a um disparo (G·a=b):");
    println!("  {:>10} {:>12} {:>16} {:>16}", "raio R", "coerência μ", "cond(G)", "‖Δa‖/‖disparo‖");
    for &r in &[2.0f32, 4.0, 6.0, 9.0, 13.0, 18.0] {
        let mut rng = Rng::new(seed);
        let atoms: Vec<Gabor2D> = (0..n).map(|_| {
            let ang=rng.range(0.0,std::f32::consts::TAU); let rr=r*rng.f32().sqrt();
            Gabor2D{cx:24.0+rr*ang.cos(),cy:24.0+rr*ang.sin(),freq:0.12,theta:0.6,phase:0.3,sigma:4.0}
        }).collect();
        let units: Vec<_> = atoms.iter().map(|a|a.render_unit(&dims)).collect();
        // matriz de Gram normalizada G_ij = ⟨ψ̂_i,ψ̂_j⟩
        let nrm: Vec<f32> = units.iter().map(|u|u.norm2().sqrt()).collect();
        let mut g = vec![vec![0.0f32;n];n];
        let mut mu=0.0; let mut npair=0;
        for i in 0..n { g[i][i]=1.0;
            for j in (i+1)..n {
                let v = if nrm[i]*nrm[j]>1e-9 { units[i].dot_patch(&units[j])/(nrm[i]*nrm[j]) } else {0.0};
                g[i][j]=v; g[j][i]=v;
                if v.abs()>1e-3 { mu+=v.abs(); npair+=1; }
            }
        }
        mu/=npair.max(1) as f32;
        // autovalores de G por iteração de potência (maior) + deflação p/ menor
        let eig_extremes=|m:&Vec<Vec<f32>>|->(f32,f32){
            let d=m.len();
            let mut v=vec![1.0f32/ (d as f32).sqrt();d];
            let mut lmax=0.0;
            for _ in 0..200 { let mut w=vec![0.0f32;d];
                for i in 0..d { for j in 0..d { w[i]+=m[i][j]*v[j]; } }
                let nn=(w.iter().map(|x|x*x).sum::<f32>()).sqrt().max(1e-12);
                lmax=nn; for i in 0..d { v[i]=w[i]/nn; }
            }
            // menor: itera em (lmax·I − G)
            let mut u=vec![1.0f32/(d as f32).sqrt();d];
            let mut shift=0.0;
            for _ in 0..200 { let mut w=vec![0.0f32;d];
                for i in 0..d { for j in 0..d { w[i]+= if i==j {lmax-m[i][j]} else {-m[i][j]}*u[j]; } }
                let nn=(w.iter().map(|x|x*x).sum::<f32>()).sqrt().max(1e-12);
                shift=nn; for i in 0..d { u[i]=w[i]/nn; }
            }
            (lmax, (lmax-shift).max(1e-6))
        };
        let (lmax,lmin)=eig_extremes(&g);
        let cond=lmax/lmin;
        // susceptibilidade: resolve G·Δa = e0 (disparo unitário no átomo 0) via
        // gradiente conjugado simples; ‖Δa‖ = resposta
        let mut da=vec![0.0f32;n]; let mut rr=vec![0.0f32;n]; rr[0]=1.0;
        let mut pp=rr.clone(); let mut rs:f32=rr.iter().map(|x|x*x).sum();
        for _ in 0..n.min(80) {
            let mut ap=vec![0.0f32;n];
            for i in 0..n { for j in 0..n { ap[i]+=g[i][j]*pp[j]; } }
            let alpha=rs/pp.iter().zip(&ap).map(|(a,b)|a*b).sum::<f32>().max(1e-12);
            for i in 0..n { da[i]+=alpha*pp[i]; rr[i]-=alpha*ap[i]; }
            let rs2:f32=rr.iter().map(|x|x*x).sum();
            if rs2<1e-10 {break;}
            let beta=rs2/rs; for i in 0..n { pp[i]=rr[i]+beta*pp[i]; } rs=rs2;
        }
        let suscep=(da.iter().map(|x|x*x).sum::<f32>()).sqrt();
        println!("  {r:>10.1} {mu:>12.3} {cond:>16.1} {suscep:>16.2}");
    }
    println!("\n  cond(G) e a susceptibilidade DISPARAM quando μ cresce (átomos colineares");
    println!("  ⇒ G quase-singular). É a MESMA colinearidade que inflou o c na Lei 6.");
    println!("  UMA HISTÓRIA: alta coerência = sinapses fortes demais = matriz mal-condicionada");
    println!("  = disparo causa resposta gigante instável. A sobrevivência foge disso");
    println!("  (auto-organiza no harmônico), mantendo a rede ESTÁVEL — subcrítica por design.");
}

pub fn run(args: &[String]) {
    let seed = args
        .iter()
        .position(|a| a == "--seed")
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(11u64);
    coupling_vs_distance(seed);
    criticality_direct(seed);
}
