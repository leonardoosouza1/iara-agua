//! `agua-lab occlusion` — o átomo NÃO-LINEAR: onde o matching pursuit não joga.
//!
//! A visão-água (§4) afirma: o diferencial da sobrevivência é o átomo
//! NÃO-DIFERENCIÁVEL / NÃO-LINEAR, onde gradiente e guloso-linear não alcançam.
//! Aqui o teste mais limpo disso: composição por OCLUSÃO (algoritmo do pintor).
//!
//! Um sinal composto por soma:  recon(x) = Σ amp_i · ψ_i(x)   ← MP resolve isso
//! Um sinal composto por oclusão: recon(x) = valor do disco MAIS À FRENTE em x
//!
//! Na oclusão NÃO há amplitude em forma fechada (o max não é linear no valor),
//! NÃO há resíduo linear (tirar um disco atrás não muda nada se está tampado),
//! e o MP nem se define — não existe ⟨resíduo, átomo⟩ que preveja a melhoria.
//!
//! Mas a SOBREVIVÊNCIA se define: a comida de um disco = quanto a composição-por-
//! oclusão melhora com ele presente. Isso é mensurável SEM gradiente e SEM
//! linearidade. Se a população resolve um alvo-de-oclusão melhor que o melhor
//! aproximador LINEAR, a tese ganha seu primeiro terreno próprio.
//!
//! Honestidade: um alvo feito de oclusão é naturalmente melhor servido por
//! oclusão — isso é esperado, não é o ponto. O PONTO é que a sobrevivência
//! CONSEGUE otimizar a composição não-linear (sem gradiente/forma-fechada),
//! enquanto o MP nem tem como tentar. Medimos os dois no MESMO alvo.

use agua_core::{psnr, Field, Rng};

/// Um disco opaco: posição, raio, valor (tom), profundidade (quem tampa quem).
#[derive(Clone)]
struct Disk {
    cx: f32,
    cy: f32,
    r: f32,
    val: f32,
    depth: f32, // maior = mais à frente
    energy: f32,
    age: u64,
}

struct OcclusionWorld {
    target: Field,
    disks: Vec<Disk>,
    rng: Rng,
    w: usize,
    h: usize,
    bg: f32,
    cfg: OccCfg,
    evals: u64,
}

#[derive(Clone)]
struct OccCfg {
    cap: usize,
    metabolism: f32,
    rain: usize,
    r_min: f32,
    r_max: f32,
    tune_tries: usize,
}

impl OcclusionWorld {
    fn new(target: Field, cfg: OccCfg, seed: u64) -> Self {
        let (h, w) = (target.dims[0], target.dims[1]);
        let bg = target.mean();
        OcclusionWorld { target, disks: Vec::new(), rng: Rng::new(seed), w, h, bg, cfg, evals: 0 }
    }

    /// Composição por oclusão: para cada pixel, o valor do disco mais à frente
    /// que o cobre; senão o fundo. É O(discos × pixels) — barato em escala pequena.
    fn render(&self) -> Field {
        let mut f = Field::zeros(&self.target.dims);
        for v in f.data.iter_mut() {
            *v = self.bg;
        }
        let mut top = vec![f32::NEG_INFINITY; self.w * self.h];
        for d in &self.disks {
            let r2 = d.r * d.r;
            let x0 = (d.cx - d.r).floor().max(0.0) as usize;
            let x1 = ((d.cx + d.r).ceil() as usize).min(self.w - 1);
            let y0 = (d.cy - d.r).floor().max(0.0) as usize;
            let y1 = ((d.cy + d.r).ceil() as usize).min(self.h - 1);
            for y in y0..=y1 {
                let dy = y as f32 - d.cy;
                for x in x0..=x1 {
                    let dx = x as f32 - d.cx;
                    if dx * dx + dy * dy <= r2 {
                        let idx = y * self.w + x;
                        if d.depth > top[idx] {
                            top[idx] = d.depth;
                            f.data[idx] = d.val;
                        }
                    }
                }
            }
        }
        f
    }

    fn mse(&self, f: &Field) -> f32 {
        self.target
            .data
            .iter()
            .zip(&f.data)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f32>()
            / self.target.len() as f32
    }

    /// Erro só na região coberta pelo disco d (pra medir sua "comida" local).
    /// Amostra um disco novo onde o erro atual é alto.
    fn error_map(&self, recon: &Field) -> Vec<f32> {
        self.target
            .data
            .iter()
            .zip(&recon.data)
            .map(|(a, b)| (a - b).abs())
            .collect()
    }

    fn rain(&mut self) {
        let recon = self.render();
        let err = self.error_map(&recon);
        let total: f32 = err.iter().sum();
        for _ in 0..self.cfg.rain {
            if self.disks.len() >= self.cfg.cap {
                break;
            }
            // amostra centro ∝ erro (hotspot)
            let mut acc = self.rng.f32() * total.max(1e-6);
            let mut idx = 0;
            for (i, e) in err.iter().enumerate() {
                acc -= *e;
                if acc <= 0.0 {
                    idx = i;
                    break;
                }
            }
            let cx = (idx % self.w) as f32;
            let cy = (idx / self.w) as f32;
            let r = self.rng.range(self.cfg.r_min, self.cfg.r_max);
            // valor = o alvo naquele ponto (bebe o que precisa), depth novo = à frente
            let val = self.target.data[idx];
            let depth = self.disks.len() as f32 + self.rng.f32();
            let d = Disk { cx, cy, r, val, depth, energy: 2.0, age: 0 };
            self.evals += 1;
            self.disks.push(d);
        }
    }

    /// Comida de cada disco = quanto a MSE piora se ele for removido (utilidade
    /// marginal na composição não-linear). Disco tampado dá comida ~0 → morre.
    fn feed_and_cull(&mut self) {
        if self.disks.is_empty() {
            return;
        }
        let base = self.mse(&self.render());
        let n = self.disks.len();
        let mut util = vec![0.0f32; n];
        for i in 0..n {
            let saved = self.disks.remove(i);
            let without = self.mse(&self.render());
            self.disks.insert(i, saved);
            util[i] = (without - base).max(0.0); // quanto ele explica
            self.evals += 1;
        }
        // energia: ganha utilidade, paga metabolismo
        for (d, u) in self.disks.iter_mut().zip(&util) {
            d.energy += *u * self.target.len() as f32 - self.cfg.metabolism;
            d.age += 1;
        }
        // morte: energia esgotada
        let mut i = self.disks.len();
        while i > 0 {
            i -= 1;
            if self.disks[i].energy <= 0.0 {
                self.disks.remove(i);
            }
        }
    }

    /// Busca local sem gradiente: cada disco tenta mexer pos/raio/valor; fica
    /// se a MSE global melhora. (O max não é diferenciável — só busca.)
    fn tune(&mut self) {
        let n = self.disks.len();
        for i in 0..n {
            if i >= self.disks.len() {
                break;
            }
            let base = self.mse(&self.render());
            let cand = {
                let d = &self.disks[i];
                Disk {
                    cx: (d.cx + self.rng.normal() * d.r * 0.4).clamp(0.0, self.w as f32 - 1.0),
                    cy: (d.cy + self.rng.normal() * d.r * 0.4).clamp(0.0, self.h as f32 - 1.0),
                    r: (d.r * (self.rng.normal() * 0.3).exp()).clamp(self.cfg.r_min, self.cfg.r_max),
                    val: (d.val + self.rng.normal() * 0.1).clamp(0.0, 1.0),
                    depth: d.depth,
                    energy: d.energy,
                    age: d.age,
                }
            };
            let old = std::mem::replace(&mut self.disks[i], cand);
            self.evals += 1;
            if self.mse(&self.render()) >= base {
                self.disks[i] = old; // não melhorou: reverte
            }
        }
    }

    fn step(&mut self) {
        self.tune();
        self.feed_and_cull();
        self.rain();
    }

    fn psnr(&self) -> f32 {
        psnr(self.mse(&self.render()), 1.0)
    }
}

/// Alvo de OCLUSÃO: K discos opacos empilhados (painter's algorithm). Bordas
/// duras — uma soma linear de gaussianas não representa sem "vazar"/ringing.
fn occlusion_target(size: usize, k: usize, seed: u64) -> Field {
    let mut rng = Rng::new(seed);
    let mut f = Field::zeros(&[size, size]);
    for v in f.data.iter_mut() {
        *v = 0.15; // fundo escuro
    }
    for _ in 0..k {
        let cx = rng.range(size as f32 * 0.15, size as f32 * 0.85);
        let cy = rng.range(size as f32 * 0.15, size as f32 * 0.85);
        let r = rng.range(size as f32 * 0.08, size as f32 * 0.22);
        let val = rng.range(0.2, 1.0);
        let r2 = r * r;
        for y in 0..size {
            for x in 0..size {
                let dx = x as f32 - cx;
                let dy = y as f32 - cy;
                if dx * dx + dy * dy <= r2 {
                    f.data[y * size + x] = val; // pinta por cima (oclusão)
                }
            }
        }
    }
    f
}

/// Melhor aproximador LINEAR: matching pursuit com átomos gaussianos (blobs),
/// composição por SOMA — o que o MP sabe fazer. Recebe o mesmo orçamento.
fn linear_mp_psnr(target: &Field, n_atoms: usize, candidates: usize, seed: u64) -> f32 {
    let (h, w) = (target.dims[0], target.dims[1]);
    let dc = target.mean();
    let mut residual = target.clone();
    for v in residual.data.iter_mut() {
        *v -= dc;
    }
    let mut recon = Field::zeros(&target.dims);
    let mut rng = Rng::new(seed);
    for _ in 0..n_atoms {
        // amostra candidatos gaussianos, escolhe o que mais explica (forma fechada)
        let mut best: Option<(f32, Vec<f32>, f32)> = None; // (explicado, patch, amp)
        for _ in 0..candidates {
            let cx = rng.range(0.0, w as f32 - 1.0);
            let cy = rng.range(0.0, h as f32 - 1.0);
            let sg = rng.range(1.5, size_frac(w));
            // patch gaussiano
            let r = (3.0 * sg) as i64;
            let x0 = (cx as i64 - r).max(0) as usize;
            let x1 = ((cx as i64 + r + 1).min(w as i64)) as usize;
            let y0 = (cy as i64 - r).max(0) as usize;
            let y1 = ((cy as i64 + r + 1).min(h as i64)) as usize;
            let mut patch = vec![0.0f32; w * h];
            let mut n2 = 0.0f32;
            let mut dot = 0.0f32;
            let inv = 1.0 / (2.0 * sg * sg);
            for y in y0..y1 {
                for x in x0..x1 {
                    let dx = x as f32 - cx;
                    let dy = y as f32 - cy;
                    let g = (-(dx * dx + dy * dy) * inv).exp();
                    patch[y * w + x] = g;
                    n2 += g * g;
                    dot += g * residual.data[y * w + x];
                }
            }
            if n2 < 1e-6 {
                continue;
            }
            let amp = dot / n2;
            let explained = amp * amp * n2;
            if best.as_ref().map_or(true, |(e, _, _)| explained > *e) {
                best = Some((explained, patch, amp));
            }
        }
        if let Some((_, patch, amp)) = best {
            for (r, (rc, p)) in residual.data.iter_mut().zip(recon.data.iter_mut().zip(&patch)) {
                *r -= amp * p;
                *rc += amp * p;
            }
        }
    }
    // MSE contra o alvo zero-médio
    let mut t0 = target.clone();
    for v in t0.data.iter_mut() {
        *v -= dc;
    }
    let mse = t0
        .data
        .iter()
        .zip(&recon.data)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f32>()
        / t0.len() as f32;
    psnr(mse, 1.0)
}

fn size_frac(w: usize) -> f32 {
    (w as f32 * 0.2).max(3.0)
}

pub fn run(args: &[String]) {
    let size = args
        .iter()
        .position(|a| a == "--size")
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(48usize);
    let ticks = args
        .iter()
        .position(|a| a == "--ticks")
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(120u64);

    println!("=== ÁTOMO NÃO-LINEAR: composição por OCLUSÃO vs melhor linear ===");
    println!("alvo = K discos opacos empilhados (bordas duras). O MP soma gaussianas;");
    println!("a sobrevivência compõe por oclusão (sem forma-fechada, sem gradiente).\n");
    println!("{:>4} {:>6} {:>16} {:>16} {:>14} {:>8}", "K", "seed", "oclusão (vida)", "linear (MP)", "Δ dB", "discos");

    let mut wins = 0;
    let mut total = 0;
    let mut sum_d = 0.0f32;
    for k in [6usize, 12, 20] {
        for seed in [11u64, 22, 33] {
            let target = occlusion_target(size, k, seed);
            // sobrevivência não-linear
            let cfg = OccCfg {
                cap: 120,
                metabolism: 0.0015,
                rain: 6,
                r_min: size as f32 * 0.05,
                r_max: size as f32 * 0.28,
                tune_tries: 1,
            };
            let mut world = OcclusionWorld::new(target.clone(), cfg, seed);
            for _ in 0..ticks {
                world.step();
            }
            let p_occ = world.psnr();
            let ndisk = world.disks.len();
            // baseline linear com orçamento >= (nº de discos)
            let p_lin = linear_mp_psnr(&target, ndisk.max(k) * 2, 24, seed ^ 0x1234);
            let delta = p_occ - p_lin;
            sum_d += delta;
            wins += (delta > 0.0) as usize;
            total += 1;
            println!(
                "{:>4} {:>6} {:>16.2} {:>16.2} {:>+14.2} {:>8}",
                k, seed, p_occ, p_lin, delta, ndisk
            );
        }
    }
    println!(
        "\noclusão vence o linear em {}/{} · Δ médio {:+.2} dB",
        wins, total, sum_d / total as f32
    );
    println!("O MP NÃO PODE ser definido na composição por oclusão (sem resíduo linear).");
    println!("A sobrevivência otimiza a composição não-linear só com pressão de comer.");
}
