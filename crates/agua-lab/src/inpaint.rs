//! `agua-lab inpaint` — COMPLETAR por sobrevivência: a SELEÇÃO vira GERAÇÃO.
//!
//! A tese-mãe (seletor-não-gerador) diz: a IARA não pinta a imagem pixel a pixel;
//! ela SELECIONA quais seres-átomo sobrevivem, e a soma dos vivos É a imagem.
//! Aqui a seleção vira geração NA SUA FRENTE: damos uma imagem com um BURACO.
//! Os seres sobrevivem comendo só o resíduo VISÍVEL (fit em forma-fechada
//! mascarado — o buraco não conta na comida), mas o CORPO Gabor deles se estende
//! pra dentro do buraco. A soma dos sobreviventes PREENCHE o vazio: ela gerou a
//! parte que faltava puramente selecionando quem harmoniza com o que se vê.
//! É a água escorrendo pro vazio.
//!
//! Nada de alvo dentro do buraco: o fit e a chuva só leem pixels visíveis. O
//! valor no buraco é 100% soma-de-corpos. Só medimos contra a verdade no fim.
//!
//! Por que é diferente de um gerador: a difusão de Laplace (inpainting clássico)
//! só SUAVIZA — borra pro buraco, não recria estrutura. A seleção estende os
//! seres OSCILATÓRIOS (ondas Gabor continuam a listra pra dentro do buraco).
//! Por isso separamos dois buracos: um sobre região LISA (blob/gradiente, onde a
//! difusão é forte) e um sobre TEXTURA (listras, o terreno da seleção).
//!
//! Baselines honestos, medidos SÓ no buraco:
//!  - média: preenche com a média do visível (piso trivial)
//!  - difusão de Laplace: inpainting clássico por PDE (competidor real, suaviza)

use crate::targets::synthetic_image;
use agua_atoms::{Atom as _, AtomBounds, Gabor2D};
use agua_core::{psnr, unflatten, Field, HotspotSampler, Patch, Rng};

/// Um ser-átomo vivo: patch (corpo Gabor recortado) + amplitude forma-fechada + comida.
struct Being {
    amp: f32,
    unit: Patch,
    food: f32,
}

/// amp forma-fechada MASCARADA: amp* = ⟨ψ, r⟩/‖ψ‖² só nos pixels VISÍVEIS.
/// comida = energia explicada no visível (o buraco não alimenta ninguém).
fn fit_amp_masked(unit: &Patch, resid: &Field, mask: &[bool]) -> Option<(f32, f32)> {
    let (mut d, mut n2) = (0.0f32, 0.0f32);
    let mut off = 0usize;
    for run in &unit.runs {
        for k in 0..run.len {
            let idx = run.field_start + k;
            if mask[idx] {
                let s = unit.data[off + k];
                d += s * resid.data[idx];
                n2 += s * s;
            }
        }
        off += run.len;
    }
    if n2 < 1e-6 {
        return None;
    }
    let amp = d / n2;
    if !amp.is_finite() {
        return None;
    }
    Some((amp, amp * amp * n2))
}

/// MSE só nos pixels selecionados (ex.: dentro de um buraco).
fn mse_sel(a: &Field, b: &Field, sel: &[bool]) -> f32 {
    let (mut s, mut c) = (0.0f32, 0usize);
    for i in 0..a.len() {
        if sel[i] {
            let d = a.data[i] - b.data[i];
            s += d * d;
            c += 1;
        }
    }
    s / c.max(1) as f32
}

/// Inpainting clássico por difusão de Laplace: pixel do buraco = média dos 4
/// vizinhos, iterado até relaxar. Competidor honesto (suaviza, não recria textura).
fn diffusion_inpaint(target: &Field, hole: &[bool], size: usize, iters: usize, seed_val: f32) -> Field {
    let n = target.len();
    let mut f = target.clone();
    for i in 0..n {
        if hole[i] {
            f.data[i] = seed_val;
        }
    }
    for _ in 0..iters {
        let cur = f.data.clone();
        for y in 0..size {
            for x in 0..size {
                let i = y * size + x;
                if !hole[i] {
                    continue;
                }
                let (mut s, mut c) = (0.0f32, 0.0f32);
                if x > 0 { s += cur[i - 1]; c += 1.0; }
                if x + 1 < size { s += cur[i + 1]; c += 1.0; }
                if y > 0 { s += cur[i - size]; c += 1.0; }
                if y + 1 < size { s += cur[i + size]; c += 1.0; }
                f.data[i] = s / c;
            }
        }
    }
    f
}

/// Salva [alvo | mascarado | reconstruído] lado a lado (PNG, cinza).
fn save_png(target: &Field, recon: &Field, mask: &[bool], size: usize, path: &str) {
    let gap = 4usize;
    let w = size * 3 + gap * 2;
    let mut img = image::GrayImage::new(w as u32, size as u32);
    let px = |v: f32| -> u8 { (v.clamp(0.0, 1.0) * 255.0) as u8 };
    for y in 0..size {
        for x in 0..size {
            let i = y * size + x;
            img.put_pixel(x as u32, y as u32, image::Luma([px(target.data[i])]));
            let m = if mask[i] { target.data[i] } else { 0.0 };
            img.put_pixel((size + gap + x) as u32, y as u32, image::Luma([px(m)]));
            img.put_pixel((2 * (size + gap) + x) as u32, y as u32, image::Luma([px(recon.data[i])]));
        }
    }
    let _ = img.save(path);
}

/// Roda o solver-por-sobrevivência com fit MASCARADO. Retorna a recon.
fn survive(target: &Field, mask: &[bool], hole: &[bool], bounds: &AtomBounds, met: f32, cap: usize, epochs: usize, seed: u64) -> (Field, usize) {
    let dims = target.dims.clone();
    let n = target.len();
    let mut rng = Rng::new(seed);
    let mut recon = Field::zeros(&dims);
    let mut resid = target.clone(); // resid = target - recon (só lido no visível)
    let mut beings: Vec<Being> = Vec::new();

    for _ep in 0..epochs {
        // GAUSS-SEIDEL: cada ser refita (mascarado) o resíduo DEPOIS dos outros
        for b in beings.iter_mut() {
            recon.add_scaled_patch(&b.unit, -b.amp);
            resid.add_scaled_patch(&b.unit, b.amp);
            let (amp, food) = fit_amp_masked(&b.unit, &resid, mask).unwrap_or((0.0, 0.0));
            b.amp = amp;
            b.food = food;
            recon.add_scaled_patch(&b.unit, b.amp);
            resid.add_scaled_patch(&b.unit, -b.amp);
        }
        // METABOLISMO: quem explica < met (no visível) morre; devolve o corpo ao resíduo
        let mut kept: Vec<Being> = Vec::with_capacity(beings.len());
        for b in beings.into_iter() {
            if b.food >= met {
                kept.push(b);
            } else {
                recon.add_scaled_patch(&b.unit, -b.amp);
                resid.add_scaled_patch(&b.unit, b.amp);
            }
        }
        beings = kept;
        // CHUVA ∝ resíduo² VISÍVEL (nunca amostra dentro do buraco)
        let mut tries = 0;
        while beings.len() < cap && tries < 6 {
            tries += 1;
            let mut resid_vis = resid.clone();
            for i in 0..n {
                if hole[i] {
                    resid_vis.data[i] = 0.0;
                }
            }
            let sampler = HotspotSampler::new(&resid_vis);
            let idx = sampler.sample(&mut rng);
            let pos: Vec<f32> = unflatten(&dims, idx).into_iter().map(|v| v as f32).collect();
            let g = Gabor2D::rain_at(&pos, bounds, &mut rng);
            let unit = g.render_unit(&dims);
            if let Some((amp, food)) = fit_amp_masked(&unit, &resid, mask) {
                if food >= met {
                    recon.add_scaled_patch(&unit, amp);
                    resid.add_scaled_patch(&unit, -amp);
                    beings.push(Being { amp, unit, food });
                }
            }
        }
    }
    (recon, beings.len())
}

// ═══════════════ PRIOR IMPORTADO: cada ladrilho é um ser que SELECIONA do codebook ═══════════════
// A visão do Leonardo: não um ser só — um SWARM de ladrilhos dividindo o "conhecimento de imagem"
// (o codebook, aprendido de um modelo/foto real) e selecionando dele. No buraco não há dado; só a
// costura com o vizinho manda → a textura da borda se PROPAGA pra dentro. É seleção+prior = geração.
const P: usize = 8; // lado do patch (ser-ladrilho)
const STRIDE: usize = 4; // sobreposição de meio-patch

/// Aprende o codebook VQ (k-means) de patches PxP — o PRIOR (conhecimento de imagem importado).
/// Se `mask` dado, só usa patches 100% VISÍVEIS (não espia o buraco); senão, usa a imagem toda.
fn learn_cb(img: &Field, mask: Option<&[bool]>, k: usize, seed: u64) -> Vec<Vec<f32>> {
    let (h, w) = (img.dims[0], img.dims[1]);
    let mut patches: Vec<Vec<f32>> = Vec::new();
    let (mut y, ok) = (0usize, |gx: usize, gy: usize| mask.map_or(true, |m| m[gy * w + gx]));
    while y + P <= h {
        let mut x = 0;
        while x + P <= w {
            if (0..P).all(|py| (0..P).all(|px| ok(x + px, y + py))) {
                let mut pt = vec![0.0f32; P * P];
                for py in 0..P {
                    for px in 0..P {
                        pt[py * P + px] = img.data[(y + py) * w + (x + px)];
                    }
                }
                patches.push(pt);
            }
            x += STRIDE;
        }
        y += STRIDE;
    }
    let mut rng = Rng::new(seed);
    let mut cb: Vec<Vec<f32>> = (0..k).map(|_| patches[rng.usize(patches.len())].clone()).collect();
    for _ in 0..12 {
        let mut sums = vec![vec![0.0f32; P * P]; k];
        let mut cnt = vec![0usize; k];
        for pt in &patches {
            let mut best = 0;
            let mut bd = f32::MAX;
            for (i, c) in cb.iter().enumerate() {
                let d: f32 = pt.iter().zip(c).map(|(a, b)| (a - b) * (a - b)).sum();
                if d < bd { bd = d; best = i; }
            }
            for (s, v) in sums[best].iter_mut().zip(pt) { *s += v; }
            cnt[best] += 1;
        }
        for i in 0..k {
            if cnt[i] > 0 {
                for (c, s) in cb[i].iter_mut().zip(&sums[i]) { *c = s / cnt[i] as f32; }
            } else {
                cb[i] = patches[rng.usize(patches.len())].clone();
            }
        }
    }
    cb
}

/// Grade de ladrilhos-seres: cada célula seleciona um índice do codebook. Comida = casar o dado
/// VISÍVEL (fixa quem vê a verdade) + costurar com os vizinhos (propaga a textura pro buraco).
struct Quilt<'a> {
    gw: usize,
    gh: usize,
    idx: Vec<usize>,
    cb: &'a [Vec<f32>],
    w: usize,
    h: usize,
    target: &'a Field,
    mask: &'a [bool],
}
impl<'a> Quilt<'a> {
    fn new(w: usize, h: usize, cb: &'a [Vec<f32>], target: &'a Field, mask: &'a [bool]) -> Self {
        let gw = (w - P) / STRIDE + 1;
        let gh = (h - P) / STRIDE + 1;
        Quilt { gw, gh, idx: vec![0; gw * gh], cb, w, h, target, mask }
    }
    fn seam(&self, ci: usize, cj: usize, cand: usize) -> f32 {
        let pc = &self.cb[cand];
        let (mut e, ov) = (0.0f32, P - STRIDE);
        if ci > 0 {
            let left = &self.cb[self.idx[cj * self.gw + (ci - 1)]];
            for py in 0..P { for k in 0..ov { let d = pc[py * P + k] - left[py * P + (STRIDE + k)]; e += d * d; } }
        }
        if ci + 1 < self.gw {
            let right = &self.cb[self.idx[cj * self.gw + (ci + 1)]];
            for py in 0..P { for k in 0..ov { let d = pc[py * P + (STRIDE + k)] - right[py * P + k]; e += d * d; } }
        }
        if cj > 0 {
            let up = &self.cb[self.idx[(cj - 1) * self.gw + ci]];
            for k in 0..ov { for px in 0..P { let d = pc[k * P + px] - up[(STRIDE + k) * P + px]; e += d * d; } }
        }
        if cj + 1 < self.gh {
            let down = &self.cb[self.idx[(cj + 1) * self.gw + ci]];
            for k in 0..ov { for px in 0..P { let d = pc[(STRIDE + k) * P + px] - down[k * P + px]; e += d * d; } }
        }
        e
    }
    /// erro contra o dado só nos pixels VISÍVEIS que a célula cobre (0 = célula 100% no buraco).
    fn data(&self, ci: usize, cj: usize, cand: usize) -> f32 {
        let pc = &self.cb[cand];
        let (ox, oy) = (ci * STRIDE, cj * STRIDE);
        let mut e = 0.0f32;
        for py in 0..P {
            for px in 0..P {
                let gi = (oy + py) * self.w + (ox + px);
                if self.mask[gi] {
                    let d = pc[py * P + px] - self.target.data[gi];
                    e += d * d;
                }
            }
        }
        e
    }
    fn survive(&mut self, sweeps: usize, seed: u64) {
        let mut rng = Rng::new(seed);
        for v in self.idx.iter_mut() { *v = rng.usize(self.cb.len()); }
        let ncells = self.gw * self.gh;
        const LAM: f32 = 1.0; // peso da costura vs dado
        for _ in 0..sweeps {
            for _ in 0..ncells {
                let ci = rng.usize(self.gw);
                let cj = rng.usize(self.gh);
                let cur = self.idx[cj * self.gw + ci];
                let mut best = cur;
                let mut bc = self.data(ci, cj, cur) + LAM * self.seam(ci, cj, cur);
                for _ in 0..10 {
                    let cand = rng.usize(self.cb.len());
                    let c = self.data(ci, cj, cand) + LAM * self.seam(ci, cj, cand);
                    if c < bc { bc = c; best = cand; }
                }
                self.idx[cj * self.gw + ci] = best;
            }
        }
    }
    fn render(&self) -> Field {
        let mut acc = vec![0.0f32; self.w * self.h];
        let mut cnt = vec![0.0f32; self.w * self.h];
        for cj in 0..self.gh {
            for ci in 0..self.gw {
                let pc = &self.cb[self.idx[cj * self.gw + ci]];
                let (ox, oy) = (ci * STRIDE, cj * STRIDE);
                for py in 0..P {
                    for px in 0..P {
                        let gi = (oy + py) * self.w + (ox + px);
                        acc[gi] += pc[py * P + px];
                        cnt[gi] += 1.0;
                    }
                }
            }
        }
        for (a, c) in acc.iter_mut().zip(&cnt) { if *c > 0.0 { *a /= *c; } }
        Field { dims: vec![self.h, self.w], data: acc }
    }
    /// custo total = dado(visível) + costura, somado em todas as células. Menor = quilt mais coerente e
    /// mais fiel ao visível → critério pra escolher o MELHOR de N realizações (mata a variância do sorteio).
    fn total_cost(&self) -> f32 {
        let mut c = 0.0f32;
        for cj in 0..self.gh {
            for ci in 0..self.gw {
                let idx = self.idx[cj * self.gw + ci];
                c += self.data(ci, cj, idx) + self.seam(ci, cj, idx);
            }
        }
        c
    }
}

/// Cena compartilhada: alvo sintético + 2 buracos (LISO sobre blob, TEXTURA sobre listras).
struct Scene {
    size: usize,
    target: Field,
    mask: Vec<bool>,
    hole: Vec<bool>,
    hole_sel: Vec<Vec<bool>>,
    names: Vec<&'static str>,
}
fn build_scene(size: usize) -> Scene {
    build_scene_img(size, None)
}
fn build_scene_img(size: usize, img: Option<&str>) -> Scene {
    let target = match img {
        Some(p) => {
            let im = image::open(p).expect("--img")
                .resize_exact(size as u32, size as u32, image::imageops::FilterType::CatmullRom)
                .to_luma8();
            Field { dims: vec![size, size], data: im.as_raw().iter().map(|v| *v as f32 / 255.0).collect() }
        }
        None => synthetic_image(size),
    };
    let n = target.len();
    let defs: [(usize, usize, usize, usize, &'static str); 2] = [
        (size * 22 / 100, size * 46 / 100, size * 20 / 100, size * 44 / 100, "LISO (blob)"),
        (size * 60 / 100, size * 82 / 100, size * 62 / 100, size * 90 / 100, "TEXTURA (listras)"),
    ];
    let mut mask = vec![true; n];
    let mut hole_sel = vec![vec![false; n]; defs.len()];
    let mut names = Vec::new();
    for (h, &(y0, y1, x0, x1, nm)) in defs.iter().enumerate() {
        names.push(nm);
        for y in y0..y1 {
            for x in x0..x1 {
                mask[y * size + x] = false;
                hole_sel[h][y * size + x] = true;
            }
        }
    }
    let hole: Vec<bool> = mask.iter().map(|&v| !v).collect();
    Scene { size, target, mask, hole, hole_sel, names }
}
fn baselines(s: &Scene) -> (Field, Field) {
    let n = s.target.len();
    let vis_mean = {
        let (mut a, mut c) = (0.0f32, 0usize);
        for i in 0..n { if s.mask[i] { a += s.target.data[i]; c += 1; } }
        a / c.max(1) as f32
    };
    let mut mean_field = s.target.clone();
    for i in 0..n { if s.hole[i] { mean_field.data[i] = vis_mean; } }
    let diff_field = diffusion_inpaint(&s.target, &s.hole, s.size, 800, vis_mean);
    (mean_field, diff_field)
}

/// Bateria COM prior: mede media · difusão · seleção-Gabor(sem prior) · SELEÇÃO+PRIOR.
fn run_prior(args: &[String]) {
    let prior_src: Option<String> = args.iter().position(|a| a == "--prior-src").and_then(|i| args.get(i + 1).cloned());
    let seeds: Vec<u64> = match crate::flag_value(args, "--seed") {
        Some(v) => vec![v.parse().unwrap()],
        None => vec![11, 22, 33],
    };
    println!("=== INPAINT + PRIOR: swarm de ladrilhos SELECIONA do codebook (conhecimento de imagem importado) ===");
    match &prior_src {
        Some(p) => println!("prior = {p} (conhecimento EXTERNO importado) → codebook k=64 dividido entre os seres-ladrilho.\n"),
        None => println!("prior = região visível do próprio alvo (auto-similaridade) → codebook k=64.\n"),
    }
    let mut agg = vec![[0.0f32; 4]; 2]; // por buraco: [media, difusao, gabor, prior]
    for (si, &seed) in seeds.iter().enumerate() {
        let r = run_prior_one(seed, prior_src.as_deref(), si == 0);
        for h in 0..2 { for j in 0..4 { agg[h][j] += r[h][j]; } }
    }
    let ns = seeds.len() as f32;
    println!("\n=== MÉDIA {} sementes — PSNR no buraco (parte gerada) ===", seeds.len());
    println!("  {:<20} {:>9} {:>9} {:>12} {:>13}", "buraco", "media", "difusao", "Gabor(s/prior)", "SELEC+PRIOR");
    for (h, nome) in ["LISO (blob)", "TEXTURA (listras)"].iter().enumerate() {
        let a = agg[h].map(|v| v / ns);
        let mark = if a[3] >= a[1] && a[3] >= a[0] && a[3] >= a[2] { " ← prior vence" } else if a[1] > a[3] { " (difusão vence)" } else { "" };
        println!("  {:<20} {:>9.2} {:>9.2} {:>12.2} {:>13.2}{}", nome, a[0], a[1], a[2], a[3], mark);
    }
    println!("→ o prior importado dá à seleção o que morder no vazio: os seres estendem a textura da borda.");
    println!("  Onde a difusão só borra (TEXTURA), a seleção+prior recria a estrutura — terreno da seleção.");
}

fn run_prior_one(seed: u64, prior_src: Option<&str>, save: bool) -> Vec<[f32; 4]> {
    let s = build_scene(96);
    // codebook = o conhecimento de imagem. Externo (importado) ou auto (região visível do alvo).
    let cb = match prior_src {
        Some(p) => {
            let img = image::open(p).expect("prior-src")
                .resize_exact(s.size as u32, s.size as u32, image::imageops::FilterType::CatmullRom)
                .to_luma8();
            let src = Field { dims: vec![s.size, s.size], data: img.as_raw().iter().map(|v| *v as f32 / 255.0).collect() };
            learn_cb(&src, None, 64, 7)
        }
        None => learn_cb(&s.target, Some(&s.mask), 64, 7),
    };
    // SWARM de ladrilhos seleciona do codebook (dado visível fixa; costura propaga pro buraco)
    let mut q = Quilt::new(s.size, s.size, &cb, &s.target, &s.mask);
    q.survive(50, seed);
    let quilt = q.render();
    let mut prior_fill = s.target.clone(); // visível = verdade; buraco = quilt
    for i in 0..s.target.len() { if s.hole[i] { prior_fill.data[i] = quilt.data[i]; } }
    // Gabor sem prior (a v1) pra comparar
    let bounds = AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
    let (gabor, _) = survive(&s.target, &s.mask, &s.hole, &bounds, 0.25, 300, 160, seed);
    let (mean_field, diff_field) = baselines(&s);

    println!("-- semente {seed}: codebook {} patches", cb.len());
    println!("  {:<20} {:>9} {:>9} {:>12} {:>13}", "buraco", "media", "difusao", "Gabor(s/prior)", "SELEC+PRIOR");
    let mut out = Vec::new();
    for h in 0..s.names.len() {
        let pm = psnr(mse_sel(&mean_field, &s.target, &s.hole_sel[h]), 1.0);
        let pd = psnr(mse_sel(&diff_field, &s.target, &s.hole_sel[h]), 1.0);
        let pg = psnr(mse_sel(&gabor, &s.target, &s.hole_sel[h]), 1.0);
        let pp = psnr(mse_sel(&prior_fill, &s.target, &s.hole_sel[h]), 1.0);
        let mark = if pp >= pd && pp >= pm && pp >= pg { " ← prior vence" } else if pd > pp { " (difusão vence)" } else { "" };
        println!("  {:<20} {:>9.2} {:>9.2} {:>12.2} {:>13.2}{}", s.names[h], pm, pd, pg, pp, mark);
        out.push([pm, pd, pg, pp]);
    }
    if save {
        save_png(&s.target, &prior_fill, &s.mask, s.size, &format!("inpaint_prior_seed{seed}.png"));
        println!("  [inpaint_prior_seed{seed}.png: alvo | com buraco | preenchido pelo swarm+prior]");
    }
    out
}

// ═══════════════ POLIGLOTA: cada ser-ladrilho ROTEIA pro expert certo (swarm+IARA) ═══════════════
// O gap do prior simples: perdia no LISO (Laplace é ~ótimo lá). Solução IARA = ROTEAMENTO por região:
// dois experts — SUAVE (difusão, ganha no liso) e TEXTURA (quilt do codebook, ganha na listra) — e cada
// pixel do buraco escolhe o expert pelo caráter LOCAL do que se VÊ na borda (variância visível difundida
// pro buraco). É o dicionário poliglota do JOURNAL aplicado: o valor é LOCAL, o roteamento emerge do dado.

/// Textura como campo de peso (mantido como alternativa de router por heurística; o padrão agora é o
/// roteamento VALIDADO por leave-out, mais robusto em foto real). Peso ∈ [0,1] (0=suave, 1=textura).
#[allow(dead_code)]
fn route_weights(s: &Scene, r: i32) -> Vec<f32> {
    let (h, w) = (s.size, s.size);
    let mut score = vec![f32::NAN; h * w];
    for y in 0..h {
        for x in 0..w {
            let (mut sum, mut sq, mut c) = (0.0f32, 0.0f32, 0usize);
            for dy in -r..=r {
                for dx in -r..=r {
                    let (yy, xx) = (y as i32 + dy, x as i32 + dx);
                    if yy < 0 || yy >= h as i32 || xx < 0 || xx >= w as i32 { continue; }
                    let j = (yy as usize) * w + (xx as usize);
                    if s.mask[j] {
                        let v = s.target.data[j];
                        sum += v;
                        sq += v * v;
                        c += 1;
                    }
                }
            }
            if c >= 3 {
                let m = sum / c as f32;
                score[y * w + x] = (sq / c as f32 - m * m).max(0.0).sqrt();
            }
        }
    }
    // difunde o score (Laplace) pros pixels sem visível suficiente (buraco fundo)
    for _ in 0..500 {
        let cur = score.clone();
        for y in 0..h {
            for x in 0..w {
                let i = y * w + x;
                if !cur[i].is_nan() { continue; }
                let (mut sm, mut cc) = (0.0f32, 0.0f32);
                for (yy, xx) in [(y as i32 - 1, x as i32), (y as i32 + 1, x as i32), (y as i32, x as i32 - 1), (y as i32, x as i32 + 1)] {
                    if yy < 0 || yy >= h as i32 || xx < 0 || xx >= w as i32 { continue; }
                    let j = (yy as usize) * w + (xx as usize);
                    if !cur[j].is_nan() { sm += cur[j]; cc += 1.0; }
                }
                if cc > 0.0 { score[i] = sm / cc; }
            }
        }
    }
    for v in score.iter_mut() { if v.is_nan() { *v = 0.0; } }
    // score → peso: faixa ESTREITA em torno de média+0.4σ do score VISÍVEL → roteamento quase-binário
    // (blend parcial de dois fills que discordam é DESTRUTIVO; melhor escolher o expert, não misturar).
    let vis: Vec<f32> = (0..h * w).filter(|&i| s.mask[i]).map(|i| score[i]).collect();
    let m = vis.iter().sum::<f32>() / vis.len() as f32;
    let sd = (vis.iter().map(|v| (v - m) * (v - m)).sum::<f32>() / vis.len() as f32).sqrt();
    let (t0, t1) = (m + 0.30 * sd, m + 0.55 * sd + 1e-6);
    score.iter().map(|&v| { let u = ((v - t0) / (t1 - t0)).clamp(0.0, 1.0); u * u * (3.0 - 2.0 * u) }).collect()
}

/// Carrega um codebook JSON (conhecimento de imagem destilado de um modelo treinado, ex.: SD-VAE).
fn load_cb_json(path: &str) -> Vec<Vec<f32>> {
    let s = std::fs::read_to_string(path).expect("--prior-json não encontrado");
    let v: serde_json::Value = serde_json::from_str(&s).expect("prior-json inválido");
    v["codebook"].as_array().expect("codebook").iter()
        .map(|row| row.as_array().unwrap().iter().map(|x| x.as_f64().unwrap() as f32).collect())
        .collect()
}

/// Expert TEXTURA por EXEMPLAR (Efros-Leung): preenche em casca-de-cebola, e pra cada pixel do buraco
/// COPIA o centro do trecho VISÍVEL cuja vizinhança (já preenchida) melhor casa. Seleção pura — copia
/// conteúdo REAL, não borra nem media → recupera textura e padrões periódicos (as listras). Auto-similar.
fn expert_exemplar(target: &Field, mask: &[bool], hole: &[bool], size: usize) -> Field {
    const W: i32 = 3; // raio da janela (patch 7×7)
    let sz = size as i32;
    let n = size * size;
    let mut f = target.clone();
    let mut filled = mask.to_vec();
    // candidatos = centros de patch 100% VISÍVEIS (fontes de onde copiar)
    let mut cand: Vec<usize> = Vec::new();
    for cy in W..sz - W {
        for cx in W..sz - W {
            let mut ok = true;
            'w: for dy in -W..=W {
                for dx in -W..=W {
                    if !mask[((cy + dy) as usize) * size + (cx + dx) as usize] { ok = false; break 'w; }
                }
            }
            if ok { cand.push((cy as usize) * size + cx as usize); }
        }
    }
    if cand.is_empty() { return expert_diffusion(target, hole, size); }
    let mut remaining = hole.iter().filter(|&&h| h).count();
    let neigh8 = [(-1i32, 0i32), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)];
    while remaining > 0 {
        // fronteira: buraco não-preenchido com vizinho preenchido; mais restritos (mais vizinhos) primeiro
        let mut border: Vec<(u32, usize)> = Vec::new();
        for i in 0..n {
            if filled[i] || !hole[i] { continue; }
            let (x, y) = ((i % size) as i32, (i / size) as i32);
            let cntf = neigh8.iter().filter(|(dx, dy)| {
                let (xx, yy) = (x + dx, y + dy);
                xx >= 0 && yy >= 0 && xx < sz && yy < sz && filled[(yy as usize) * size + xx as usize]
            }).count() as u32;
            if cntf > 0 { border.push((cntf, i)); }
        }
        if border.is_empty() { break; }
        border.sort_by(|a, b| b.0.cmp(&a.0));
        for &(_, i) in &border {
            if filled[i] { continue; }
            let (cx0, cy0) = ((i % size) as i32, (i / size) as i32);
            let (mut best, mut bestc) = (f32::MAX, cand[0]);
            for &c in &cand {
                let (ccx, ccy) = ((c % size) as i32, (c / size) as i32);
                let mut ssd = 0.0f32;
                for dy in -W..=W {
                    for dx in -W..=W {
                        let (nx, ny) = (cx0 + dx, cy0 + dy);
                        if nx < 0 || ny < 0 || nx >= sz || ny >= sz { continue; }
                        let ni = (ny as usize) * size + nx as usize;
                        if filled[ni] {
                            let d = f.data[ni] - f.data[((ccy + dy) as usize) * size + (ccx + dx) as usize];
                            ssd += d * d;
                        }
                    }
                }
                if ssd < best { best = ssd; bestc = c; }
            }
            f.data[i] = target.data[bestc];
            filled[i] = true;
            remaining -= 1;
        }
    }
    f
}

/// Expert TEXTURA: EXEMPLAR (auto-similar, alta qualidade) quando sem prior; QUILT de codebook quando
/// `cb_ext` = conhecimento de imagem importado (modelo treinado).
fn expert_texture(target: &Field, mask: &[bool], hole: &[bool], size: usize, seed: u64, cb_ext: Option<&Vec<Vec<f32>>>) -> Field {
    if cb_ext.is_none() { return expert_exemplar(target, mask, hole, size); }
    let cb_owned = if cb_ext.is_none() { Some(learn_cb(target, Some(mask), 64, 7)) } else { None };
    let cb: &[Vec<f32>] = cb_ext.map(|v| v.as_slice()).unwrap_or_else(|| cb_owned.as_ref().unwrap().as_slice());
    // BEST-OF-N: o quilt estocástico converge a ótimos-locais diferentes por sorteio (variância 15.6→20.0).
    // Roda N realizações e fica com a de MENOR custo total (dado+costura) = mais coerente e mais fiel ao
    // visível. Mata a variância → o verificador passa a ver a real qualidade do expert.
    const NBEST: usize = 4;
    let mut best: Option<(f32, Field)> = None;
    for n in 0..NBEST {
        let mut q = Quilt::new(size, size, cb, target, mask);
        q.survive(50, seed.wrapping_add(n as u64 * 0x9E37));
        let cost = q.total_cost();
        if best.as_ref().map_or(true, |(bc, _)| cost < *bc) { best = Some((cost, q.render())); }
    }
    let quilt = best.unwrap().1;
    let mut f = target.clone();
    for i in 0..target.len() { if hole[i] { f.data[i] = quilt.data[i]; } }
    f
}

/// Expert SUAVE: difusão de Laplace (preenche o buraco com a média do visível como semente).
fn expert_diffusion(target: &Field, hole: &[bool], size: usize) -> Field {
    let n = size * size;
    let (mut a, mut c) = (0.0f32, 0usize);
    for i in 0..n { if !hole[i] { a += target.data[i]; c += 1; } }
    diffusion_inpaint(target, hole, size, 800, a / c.max(1) as f32)
}

/// Componentes conexos do buraco (4-conexo) — cada um roteia pro seu expert independentemente.
fn hole_components(hole: &[bool], size: usize) -> Vec<Vec<usize>> {
    let n = size * size;
    let mut comp = vec![usize::MAX; n];
    let mut comps: Vec<Vec<usize>> = Vec::new();
    for start in 0..n {
        if !hole[start] || comp[start] != usize::MAX { continue; }
        let id = comps.len();
        comp[start] = id;
        let mut stack = vec![start];
        let mut cells = Vec::new();
        while let Some(i) = stack.pop() {
            cells.push(i);
            let (x, y) = ((i % size) as i32, (i / size) as i32);
            for (dx, dy) in [(-1, 0), (1, 0), (0, -1), (0, 1)] {
                let (xx, yy) = (x + dx, y + dy);
                if xx < 0 || yy < 0 || xx >= size as i32 || yy >= size as i32 { continue; }
                let j = (yy as usize) * size + (xx as usize);
                if hole[j] && comp[j] == usize::MAX { comp[j] = id; stack.push(j); }
            }
        }
        comps.push(cells);
    }
    comps
}

/// ROTEAMENTO VALIDADO (verificador, não heurística): reposiciona uma RÉPLICA do buraco (mesma forma)
/// numa região VISÍVEL vizinha, esconde-a, e faz cada expert PREENCHÊ-LA a partir do resto — MESMA
/// tarefa/dificuldade do buraco real (preencher de dentro, não um anel externo fácil). Escolhe quem
/// acerta mais a réplica (verdade conhecida). Só usa dado visível — honesto. (textura_vence?, ed, et).
fn validated_route(target: &Field, mask: &[bool], comp: &[usize], size: usize, seed: u64, cb_ext: Option<&Vec<Vec<f32>>>) -> (bool, f32, f32) {
    let n = size * size;
    let sz = size as i32;
    let (xs, ys): (Vec<i32>, Vec<i32>) = comp.iter().map(|&i| ((i % size) as i32, (i / size) as i32)).unzip();
    let (xmin, xmax) = (*xs.iter().min().unwrap(), *xs.iter().max().unwrap());
    let (ymin, ymax) = (*ys.iter().min().unwrap(), *ys.iter().max().unwrap());
    let (bw, bh) = (xmax - xmin + 1, ymax - ymin + 1);
    let rel: Vec<(i32, i32)> = comp.iter().map(|&i| ((i % size) as i32 - xmin, (i / size) as i32 - ymin)).collect();
    let (cx, cy) = ((xmin + xmax) / 2, (ymin + ymax) / 2);
    // acha a colocação VISÍVEL mais próxima do buraco (estatística local parecida) que não toca o buraco
    let try_place = |ox: i32, oy: i32| -> Option<Vec<usize>> {
        let mut cells = Vec::with_capacity(rel.len());
        for &(rx, ry) in &rel {
            let (x, y) = (ox + rx, oy + ry);
            if x < 0 || y < 0 || x >= sz || y >= sz { return None; }
            let j = (y as usize) * size + (x as usize);
            if !mask[j] { return None; } // precisa ser 100% visível (exclui o próprio buraco)
            cells.push(j);
        }
        Some(cells)
    };
    // coleta as réplicas VÁLIDAS mais próximas e usa as K melhores NÃO-sobrepostas (média do erro →
    // robusto: uma vizinhança sozinha pode não representar o miolo do buraco).
    let mut cands: Vec<(i32, Vec<usize>)> = Vec::new();
    let mut oy = 0;
    while oy + bh <= sz {
        let mut ox = 0;
        while ox + bw <= sz {
            if let Some(cells) = try_place(ox, oy) {
                let dist = (ox + bw / 2 - cx).abs() + (oy + bh / 2 - cy).abs();
                cands.push((dist, cells));
            }
            ox += 2;
        }
        oy += 2;
    }
    // referência de TEXTURA da borda do buraco: variância dos pixels visíveis vizinhos ao componente.
    let var_of = |vals: &[f32]| -> f32 {
        if vals.is_empty() { return 0.0; }
        let m = vals.iter().sum::<f32>() / vals.len() as f32;
        vals.iter().map(|v| (v - m) * (v - m)).sum::<f32>() / vals.len() as f32
    };
    let mut bvals: Vec<f32> = Vec::new();
    for &i in comp {
        let (x, y) = ((i % size) as i32, (i / size) as i32);
        for dy in -3..=3 { for dx in -3..=3 {
            let (xx, yy) = (x + dx, y + dy);
            if xx < 0 || yy < 0 || xx >= sz || yy >= sz { continue; }
            let j = (yy as usize) * size + (xx as usize);
            if mask[j] { bvals.push(target.data[j]); }
        } }
    }
    let bref = var_of(&bvals);
    // escolhe a réplica cuja TEXTURA (variância) CASA com a borda do buraco → representa o buraco de
    // verdade (a mais próxima podia cair numa região suave e enganar o roteamento). Desempata por perto.
    let mut best: Option<(f32, Vec<usize>)> = None;
    for (dist, cells) in &cands {
        let vals: Vec<f32> = cells.iter().map(|&i| target.data[i]).collect();
        let score = (var_of(&vals) - bref).abs() + 3e-5 * *dist as f32;
        if best.as_ref().map_or(true, |(bs, _)| score < *bs) { best = Some((score, cells.clone())); }
    }
    let Some((_, test_var)) = best else { return (false, 0.0, 0.0) }; // sem lugar → default SUAVE (seguro)
    // avalia um conjunto-teste: esconde-o (+ o buraco), cada expert prevê, retorna (erro_dif, erro_tex).
    let eval = |test: &[usize]| -> (f32, f32) {
        let mut hole2 = vec![false; n];
        for i in 0..n { if !mask[i] { hole2[i] = true; } }
        let mut mask2 = mask.to_vec();
        for &i in test { hole2[i] = true; mask2[i] = false; }
        let diff2 = expert_diffusion(target, &hole2, size);
        let tex2 = expert_texture(target, &mask2, &hole2, size, seed ^ 0x99, cb_ext);
        let (mut ed, mut et) = (0.0f32, 0.0f32);
        for &i in test {
            let t = target.data[i];
            ed += (diff2.data[i] - t) * (diff2.data[i] - t);
            et += (tex2.data[i] - t) * (tex2.data[i] - t);
        }
        let d = test.len().max(1) as f32;
        (ed / d, et / d)
    };
    // MARGEM proporcional ao TAMANHO do buraco: textura-quilt (patch 8px) precisa de ÁREA pra valer.
    let area = comp.len() as f32;
    let margin = 0.60 + 0.32 * ((area - 64.0) / 512.0).clamp(0.0, 1.0);
    let (edv, etv) = eval(&test_var);
    // Roteia pra EXEMPLAR só quando AMBOS: (1) a BORDA é texturada (bref alto) — barra o falso-positivo
    // no liso, onde a difusão é imbatível; (2) o exemplar BATE a difusão na réplica com margem — barra o
    // exemplar quando falta contexto (buraco grande sobre textura: sem trecho pra copiar, ele perde).
    let tex_wins = bref > 0.013 && etv < edv * margin;
    if std::env::var("ROUTE_DBG").is_ok() {
        eprintln!("  [route] area={:.0} bref={:.4} edv={:.5} etv={:.5} m={:.2} → tex={}", area, bref, edv, etv, margin, tex_wins);
    }
    (tex_wins, edv, etv)
}

/// Cena com corrupção ESPALHADA: ~40 buracos quadrados aleatórios (inpainting realista). Máscara fixa
/// por `mask_seed` (a corrupção não muda com a semente do solver). Um único grupo de avaliação.
fn build_scene_scatter(size: usize, img: Option<&str>, mask_seed: u64) -> Scene {
    let base = build_scene_img(size, img);
    let target = base.target;
    let n = target.len();
    let mut rng = Rng::new(mask_seed);
    let mut mask = vec![true; n];
    let mut sel = vec![false; n];
    for _ in 0..40 {
        let hw = 5 + rng.usize(5);
        let hh = 5 + rng.usize(5);
        let x0 = 1 + rng.usize(size.saturating_sub(hw + 2).max(1));
        let y0 = 1 + rng.usize(size.saturating_sub(hh + 2).max(1));
        for y in y0..(y0 + hh).min(size) {
            for x in x0..(x0 + hw).min(size) {
                mask[y * size + x] = false;
                sel[y * size + x] = true;
            }
        }
    }
    let hole: Vec<bool> = mask.iter().map(|&v| !v).collect();
    Scene { size, target, mask, hole, hole_sel: vec![sel], names: vec!["espalhado (~40 buracos)"] }
}

/// Cena com UM buraco GRANDE central (~60% da área) — contexto escasso: a auto-similaridade FALHA
/// (pouco visível pra aprender), o cenário onde o conhecimento de um MODELO TREINADO deve pagar.
fn build_scene_bighole(size: usize, img: Option<&str>) -> Scene {
    let base = build_scene_img(size, img);
    let target = base.target;
    let n = target.len();
    let (mut mask, mut sel) = (vec![true; n], vec![false; n]);
    let (lo, hi) = (size * 12 / 100, size * 88 / 100); // buraco central grande
    for y in lo..hi {
        for x in lo..hi {
            mask[y * size + x] = false;
            sel[y * size + x] = true;
        }
    }
    let hole: Vec<bool> = mask.iter().map(|&v| !v).collect();
    Scene { size, target, mask, hole, hole_sel: vec![sel], names: vec!["buraco GRANDE (contexto escasso)"] }
}

fn run_polyglot(args: &[String]) {
    let seeds: Vec<u64> = match crate::flag_value(args, "--seed") {
        Some(v) => vec![v.parse().unwrap()],
        None => vec![11, 22, 33],
    };
    let img: Option<String> = args.iter().position(|a| a == "--img").and_then(|i| args.get(i + 1).cloned());
    let scatter = args.iter().any(|a| a == "--scatter");
    let bighole = args.iter().any(|a| a == "--bighole");
    let prior_json: Option<String> = args.iter().position(|a| a == "--prior-json").and_then(|i| args.get(i + 1).cloned());
    let cb_ext: Option<Vec<Vec<f32>>> = prior_json.as_deref().map(load_cb_json);
    let names: Vec<&str> = if scatter { vec!["espalhado (~40 buracos)"] } else if bighole { vec!["buraco GRANDE (contexto escasso)"] } else { vec!["LISO (blob)", "TEXTURA (listras)"] };
    println!("=== INPAINT POLIGLOTA: cada ser-ladrilho ROTEIA pro expert certo (swarm+IARA) ===");
    let cenario = if scatter { " · corrupção ESPALHADA (buracos aleatórios)" } else { "" };
    let prior_lbl = match &cb_ext { Some(c) => format!("codebook de MODELO TREINADO ({} patches, {})", c.len(), prior_json.as_deref().unwrap()), None => "codebook auto (k-means do visível)".into() };
    match &img {
        Some(p) => println!("alvo = {p} (FOTO REAL){cenario} · prior TEXTURA = {prior_lbl}\n"),
        None => println!("alvo sintético{cenario} · prior TEXTURA = {prior_lbl}\n"),
    }
    let mut agg: Vec<[f32; 4]> = Vec::new();
    for (si, &seed) in seeds.iter().enumerate() {
        let r = run_polyglot_one(seed, img.as_deref(), scatter, bighole, cb_ext.as_ref(), si == 0);
        if agg.is_empty() { agg = vec![[0.0f32; 4]; r.len()]; }
        for h in 0..r.len() { for j in 0..4 { agg[h][j] += r[h][j]; } }
    }
    let ns = seeds.len() as f32;
    let nh = agg.len() as f32;
    println!("\n=== MÉDIA {} sementes — PSNR no buraco (parte gerada) ===", seeds.len());
    println!("  {:<22} {:>9} {:>9} {:>10} {:>12} {:>8}", "buraco", "media", "difusao", "textura", "POLIGLOTA", "oráculo");
    let (mut sd, mut st, mut sp, mut so) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    for (h, nome) in names.iter().enumerate() {
        let a = agg[h].map(|v| v / ns);
        let oracle = a[1].max(a[2]); // melhor expert SE soubéssemos a verdade do buraco (teto inalcançável)
        sd += a[1]; st += a[2]; sp += a[3]; so += oracle;
        let mark = if a[3] >= oracle - 0.05 { " ✓rota ótima" } else { " (perde só p/ oráculo)" };
        println!("  {:<22} {:>9.2} {:>9.2} {:>10.2} {:>12.2} {:>8.2}{}", nome, a[0], a[1], a[2], a[3], oracle, mark);
    }
    // critério HONESTO: bater a melhor ESTRATÉGIA FIXA implantável (sempre-difusão ou sempre-textura)
    let (ad, at, ap, ao) = (sd / nh, st / nh, sp / nh, so / nh);
    let best_fixed = ad.max(at);
    println!("\n  média imagem → sempre-difusão {ad:.2} · sempre-textura {at:.2} · POLIGLOTA {ap:.2} · oráculo {ao:.2}");
    if ap >= best_fixed - 0.10 {
        println!("  ✅ a poliglota ROTEADA ≥ melhor estratégia FIXA implantável ({best_fixed:.2}) e nunca cai abaixo da difusão.");
        if ao - ap > 0.1 { println!("     (gap p/ o oráculo {:.2} dB = onde a vizinhança visível não revela o expert do buraco — limite do roteamento CEGO.)", ao - ap); }
    } else {
        println!("  ⚠️ abaixo da melhor estratégia fixa — investigar.");
    }
}

/// Comparação COLORIDA [original | com buraco (vermelho) | difusão | POLIGLOTA], upscale nearest.
fn save_compare(target: &Field, mask: &[bool], diff: &Field, poly: &Field, size: usize, scale: usize, path: &str) {
    let gap = 4usize;
    let pw = size * scale;
    let w = pw * 4 + gap * 3;
    let mut img = image::RgbImage::from_pixel(w as u32, pw as u32, image::Rgb([255, 255, 255]));
    let g = |v: f32| -> u8 { (v.clamp(0.0, 1.0) * 255.0) as u8 };
    for oy in 0..pw {
        let sy = oy / scale;
        for oxp in 0..pw {
            let i = sy * size + (oxp / scale);
            let (gv, dv, pv) = (g(target.data[i]), g(diff.data[i]), g(poly.data[i]));
            img.put_pixel(oxp as u32, oy as u32, image::Rgb([gv, gv, gv]));
            let p1 = image::Rgb(if !mask[i] { [225, 35, 35] } else { [gv, gv, gv] });
            img.put_pixel((pw + gap + oxp) as u32, oy as u32, p1);
            img.put_pixel((2 * (pw + gap) + oxp) as u32, oy as u32, image::Rgb([dv, dv, dv]));
            img.put_pixel((3 * (pw + gap) + oxp) as u32, oy as u32, image::Rgb([pv, pv, pv]));
        }
    }
    let _ = img.save(path);
}

fn run_polyglot_one(seed: u64, img: Option<&str>, scatter: bool, bighole: bool, cb_ext: Option<&Vec<Vec<f32>>>, save: bool) -> Vec<[f32; 4]> {
    let s = if scatter { build_scene_scatter(96, img, 999) } else if bighole { build_scene_bighole(96, img) } else { build_scene_img(96, img) };
    let (mean_field, diff_field) = baselines(&s); // expert SUAVE = difusão
    let tex_field = expert_texture(&s.target, &s.mask, &s.hole, s.size, seed, cb_ext); // expert TEXTURA
    // ROTEAMENTO VALIDADO: por componente do buraco, cada expert prevê um anel visível escondido
    // (leave-out) e vence quem reconstrói melhor o observável. Hard-fill (sem blend destrutivo).
    let comps = hole_components(&s.hole, s.size);
    let mut poly = s.target.clone();
    let mut ntex = 0usize;
    for comp in &comps {
        let (tex_wins, _ed, _et) = validated_route(&s.target, &s.mask, comp, s.size, seed, cb_ext);
        let src = if tex_wins { &tex_field } else { &diff_field };
        if tex_wins { ntex += 1; }
        for &i in comp { poly.data[i] = src.data[i]; }
    }
    println!("-- semente {seed}: {} componentes · {ntex}→TEXTURA {}→SUAVE (validado por anel)", comps.len(), comps.len() - ntex);
    println!("  {:<20} {:>9} {:>9} {:>10} {:>12}", "buraco", "media", "difusao", "textura", "POLIGLOTA");
    let mut out = Vec::new();
    for h in 0..s.names.len() {
        let pm = psnr(mse_sel(&mean_field, &s.target, &s.hole_sel[h]), 1.0);
        let pd = psnr(mse_sel(&diff_field, &s.target, &s.hole_sel[h]), 1.0);
        let pt = psnr(mse_sel(&tex_field, &s.target, &s.hole_sel[h]), 1.0);
        let pp = psnr(mse_sel(&poly, &s.target, &s.hole_sel[h]), 1.0);
        println!("  {:<20} {:>9.2} {:>9.2} {:>10.2} {:>12.2}", s.names[h], pm, pd, pt, pp);
        out.push([pm, pd, pt, pp]);
    }
    if save {
        let p = format!("inpaint_poliglota_seed{seed}.png");
        save_compare(&s.target, &s.mask, &diff_field, &poly, s.size, 4, &p);
        println!("  [{p}: original | com buraco | difusão | POLIGLOTA]");
    }
    out
}

pub fn run(args: &[String]) {
    if args.iter().any(|a| a == "--poliglota") {
        run_polyglot(args);
        return;
    }
    if args.iter().any(|a| a == "--prior" || a == "--prior-src") {
        run_prior(args);
        return;
    }
    println!("=== INPAINT: a SELEÇÃO vira GERAÇÃO — seres preenchem o buraco só selecionando quem harmoniza ===");
    println!("alvo sintético 96×96; fit MASCARADO (o buraco não alimenta); corpo Gabor se estende pro vazio.\n");
    let seeds: Vec<u64> = match crate::flag_value(args, "--seed") {
        Some(v) => vec![v.parse().unwrap()],
        None => vec![11, 22, 33],
    };
    let mut agg: Vec<(f32, f32, f32)> = vec![(0.0, 0.0, 0.0); 2]; // por buraco: (media, difusao, selecao)
    for (si, &seed) in seeds.iter().enumerate() {
        let r = run_one(seed, si == 0);
        for h in 0..2 {
            agg[h].0 += r[h].0;
            agg[h].1 += r[h].1;
            agg[h].2 += r[h].2;
        }
    }
    let ns = seeds.len() as f32;
    println!("\n=== MÉDIA {} sementes — PSNR no buraco (parte gerada) ===", seeds.len());
    println!("  {:<20} {:>10} {:>10} {:>12}", "buraco", "media", "difusao", "SELECAO");
    for (h, nome) in ["LISO (blob)", "TEXTURA (listras)"].iter().enumerate() {
        let (m, d, s) = (agg[h].0 / ns, agg[h].1 / ns, agg[h].2 / ns);
        let mark = if s >= d && s >= m { " ← seleção vence" } else if d > s { " (difusão vence)" } else { "" };
        println!("  {:<20} {:>10.2} {:>10.2} {:>12.2}{}", nome, m, d, s, mark);
    }
    println!("→ a seleção-por-sobrevivência GEROU o que faltava sem nunca ver o buraco. Onde há TEXTURA,");
    println!("  estender os seres oscilatórios recria a estrutura que a difusão só borra.");
}

/// Uma semente: roda, imprime a tabela, salva PNG (se `save`). Retorna [(media,difusao,selecao); 2 buracos].
fn run_one(seed: u64, save: bool) -> Vec<(f32, f32, f32)> {
    let size = 96usize;
    let target = synthetic_image(size);
    let n = target.len();

    // dois buracos: H0 sobre região LISA (blob/gradiente), H1 sobre TEXTURA (listras)
    let holes_def: [(usize, usize, usize, usize, &str); 2] = [
        (size * 22 / 100, size * 46 / 100, size * 20 / 100, size * 44 / 100, "LISO (blob)"),
        (size * 60 / 100, size * 82 / 100, size * 62 / 100, size * 90 / 100, "TEXTURA (listras)"),
    ];
    let mut mask = vec![true; n]; // true = visível
    let mut hole_sel: Vec<Vec<bool>> = vec![vec![false; n]; holes_def.len()];
    for (h, &(y0, y1, x0, x1, _)) in holes_def.iter().enumerate() {
        for y in y0..y1 {
            for x in x0..x1 {
                mask[y * size + x] = false;
                hole_sel[h][y * size + x] = true;
            }
        }
    }
    let hole: Vec<bool> = mask.iter().map(|&v| !v).collect();

    let bounds = AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
    let met = 0.25f32;
    let (recon, pop) = survive(&target, &mask, &hole, &bounds, met, 300, 160, seed);

    // baselines
    let vis_mean = {
        let (mut s, mut c) = (0.0f32, 0usize);
        for i in 0..n {
            if mask[i] {
                s += target.data[i];
                c += 1;
            }
        }
        s / c.max(1) as f32
    };
    let mut mean_field = target.clone();
    for i in 0..n {
        if hole[i] {
            mean_field.data[i] = vis_mean;
        }
    }
    let diff_field = diffusion_inpaint(&target, &hole, size, 800, vis_mean);

    println!("-- semente {seed}: {pop} seres (emergiu do metabolismo) · PSNR visível {:.2} dB (sanidade)", psnr(mse_sel(&recon, &target, &mask), 1.0));
    println!("  {:<20} {:>10} {:>10} {:>12}", "buraco", "media", "difusao", "SELECAO(vida)");
    let mut out = Vec::with_capacity(holes_def.len());
    for (h, &(_, _, _, _, nome)) in holes_def.iter().enumerate() {
        let pm = psnr(mse_sel(&mean_field, &target, &hole_sel[h]), 1.0);
        let pd = psnr(mse_sel(&diff_field, &target, &hole_sel[h]), 1.0);
        let ps = psnr(mse_sel(&recon, &target, &hole_sel[h]), 1.0);
        let mark = if ps >= pd && ps >= pm { " ← seleção vence" } else if pd > ps { " (difusão vence)" } else { "" };
        println!("  {:<20} {:>10.2} {:>10.2} {:>12.2}{}", nome, pm, pd, ps, mark);
        out.push((pm, pd, ps));
    }
    if save {
        let path = format!("inpaint_seed{seed}.png");
        save_png(&target, &recon, &mask, size, &path);
        println!("  [{path}: alvo | com buraco | preenchido pelos seres]");
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    // INVARIANTE 1: o roteamento nunca deixa a poliglota abaixo da difusão (segurança do baseline forte).
    #[test]
    fn poliglota_nunca_pior_que_difusao() {
        let s = build_scene(96);
        let (_, diff) = baselines(&s);
        let tex = expert_texture(&s.target, &s.mask, &s.hole, s.size, 11, None);
        let mut poly = s.target.clone();
        for comp in &hole_components(&s.hole, s.size) {
            let (tw, _, _) = validated_route(&s.target, &s.mask, comp, s.size, 11, None);
            let src = if tw { &tex } else { &diff };
            for &i in comp { poly.data[i] = src.data[i]; }
        }
        let pd = psnr(mse_sel(&diff, &s.target, &s.hole), 1.0);
        let pp = psnr(mse_sel(&poly, &s.target, &s.hole), 1.0);
        assert!(pp >= pd - 0.5, "poliglota {pp:.2} caiu abaixo da difusão {pd:.2}");
    }
    // INVARIANTE 2: determinismo por seed (regra da casa).
    #[test]
    fn determinismo_por_seed() {
        let s = build_scene(96);
        let a = expert_texture(&s.target, &s.mask, &s.hole, s.size, 7, None);
        let b = expert_texture(&s.target, &s.mask, &s.hole, s.size, 7, None);
        assert_eq!(a.data, b.data);
    }
    // INVARIANTE 3: fit mascarado é a projeção forma-fechada (amp minimiza o erro no visível).
    #[test]
    fn fit_mascarado_e_otimo() {
        let s = build_scene(48);
        let g = Gabor2D::rain_at(&[24.0, 24.0], &AtomBounds { freq_min: 0.05, freq_max: 0.3, sigma_min: 3.0, sigma_max: 8.0 }, &mut Rng::new(1));
        let unit = g.render_unit(&s.target.dims);
        let (amp, _) = fit_amp_masked(&unit, &s.target, &s.mask).unwrap();
        let e = |a: f32| -> f32 {
            let mut acc = 0.0f32; let mut off = 0usize;
            for run in &unit.runs { for k in 0..run.len {
                let idx = run.field_start + k;
                if s.mask[idx] { let d = a * unit.data[off + k] - s.target.data[idx]; acc += d * d; }
            } off += run.len; }
            acc
        };
        assert!(e(amp) <= e(amp + 0.05) && e(amp) <= e(amp - 0.05), "amp forma-fechada não é o mínimo");
    }
}
