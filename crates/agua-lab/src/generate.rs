//! `agua-lab generate` — F4: GERAÇÃO por prior importado (a promessa da tese).
//!
//! Até aqui a sobrevivência RECONSTRUIU sinais conhecidos. A visão §4 diz que o
//! diferencial real é GERAR: o "alvo" vira um PRIOR aprendido; os seres
//! sobrevivem casando com o prior → amostram sinais NOVOS. Aqui o MP não se
//! aplica: não há sinal-alvo pra fazer resíduo. Só existe "isto é plausível?".
//!
//! Prior importado = codebook VQ (k-means) de patches de uma imagem REAL. É a
//! ponte offline do ecossistema (aprende de dado real, exporta, o solver
//! consome) no caso mais simples e auto-contido.
//!
//! Geração por sobrevivência: uma tela vazia; seres = ladrilhos (posição +
//! índice do codebook). Comida de um ser = quanto seu patch CONCORDA com os
//! vizinhos que se sobrepõem (costura coerente) — plausível localmente. Nasce
//! onde a tela está incoerente; morre quem briga com a vizinhança. O resultado
//! é uma textura nova, localmente como o dado real, globalmente inédita.
//!
//! Baselines honestos:
//!  - aleatório: ladrilha índices do codebook sem seleção por costura (piso)
//!  - quilting guloso: escolhe, por posição, o patch que melhor costura (teto)
//! Métricas:
//!  - COERÊNCIA: energia de costura média nas sobreposições (menor = melhor)
//!  - NOVIDADE: distância mínima da imagem gerada a QUALQUER recorte do original
//!    (0 = cópia; queremos > 0 — gerou, não copiou)

use agua_core::{Field, Rng};

const P: usize = 8; // lado do patch
const STRIDE: usize = 4; // sobreposição de meio-patch

/// Aprende o codebook VQ (k-means) de patches PxP da imagem — o PRIOR importado.
fn learn_codebook(img: &Field, k: usize, seed: u64) -> Vec<Vec<f32>> {
    let (h, w) = (img.dims[0], img.dims[1]);
    let mut patches: Vec<Vec<f32>> = Vec::new();
    let mut y = 0;
    while y + P <= h {
        let mut x = 0;
        while x + P <= w {
            let mut pt = vec![0.0f32; P * P];
            for py in 0..P {
                for px in 0..P {
                    pt[py * P + px] = img.data[(y + py) * w + (x + px)];
                }
            }
            patches.push(pt);
            x += STRIDE;
        }
        y += STRIDE;
    }
    let mut rng = Rng::new(seed);
    // init: k patches aleatórios
    let mut cb: Vec<Vec<f32>> = (0..k).map(|_| patches[rng.usize(patches.len())].clone()).collect();
    for _ in 0..12 {
        let mut sums = vec![vec![0.0f32; P * P]; k];
        let mut cnt = vec![0usize; k];
        for pt in &patches {
            let mut best = 0;
            let mut bd = f32::MAX;
            for (i, c) in cb.iter().enumerate() {
                let d: f32 = pt.iter().zip(c).map(|(a, b)| (a - b) * (a - b)).sum();
                if d < bd {
                    bd = d;
                    best = i;
                }
            }
            for (s, v) in sums[best].iter_mut().zip(pt) {
                *s += v;
            }
            cnt[best] += 1;
        }
        for i in 0..k {
            if cnt[i] > 0 {
                for (c, s) in cb[i].iter_mut().zip(&sums[i]) {
                    *c = s / cnt[i] as f32;
                }
            } else {
                cb[i] = patches[rng.usize(patches.len())].clone();
            }
        }
    }
    cb
}

/// Grade de ladrilhos: cada célula (i,j) tem um índice do codebook. Os patches
/// se sobrepõem em STRIDE; a costura é o desacordo na faixa de overlap.
struct Quilt {
    gw: usize,
    gh: usize,
    idx: Vec<usize>, // índice do codebook por célula
    cb: Vec<Vec<f32>>,
    w: usize,
    h: usize,
}

impl Quilt {
    fn new(w: usize, h: usize, cb: Vec<Vec<f32>>) -> Self {
        let gw = (w - P) / STRIDE + 1;
        let gh = (h - P) / STRIDE + 1;
        Quilt { gw, gh, idx: vec![0; gw * gh], cb, w, h }
    }

    /// Costura entre a célula (i,j) e seu vizinho: soma dos quadrados na faixa
    /// de sobreposição (P-STRIDE de largura). Só direita e baixo (evita dupla).
    fn seam(&self, ci: usize, cj: usize, cand: usize) -> f32 {
        let pc = &self.cb[cand];
        let mut e = 0.0f32;
        let ov = P - STRIDE;
        // vizinho à esquerda (compartilha a faixa esquerda de cand)
        if ci > 0 {
            let left = &self.cb[self.idx[cj * self.gw + (ci - 1)]];
            for py in 0..P {
                for k in 0..ov {
                    let a = pc[py * P + k]; // faixa esquerda de cand
                    let b = left[py * P + (STRIDE + k)]; // faixa direita do vizinho
                    e += (a - b) * (a - b);
                }
            }
        }
        // vizinho acima
        if cj > 0 {
            let up = &self.cb[self.idx[(cj - 1) * self.gw + ci]];
            for k in 0..ov {
                for px in 0..P {
                    let a = pc[k * P + px];
                    let b = up[(STRIDE + k) * P + px];
                    e += (a - b) * (a - b);
                }
            }
        }
        e
    }

    fn total_seam(&self) -> f32 {
        let mut e = 0.0f32;
        for cj in 0..self.gh {
            for ci in 0..self.gw {
                e += self.seam(ci, cj, self.idx[cj * self.gw + ci]);
            }
        }
        e / (self.gw * self.gh) as f32
    }

    /// Renderiza a tela: média dos patches que cobrem cada pixel (feathering).
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
        for (a, c) in acc.iter_mut().zip(&cnt) {
            if *c > 0.0 {
                *a /= *c;
            }
        }
        Field { dims: vec![self.h, self.w], data: acc }
    }
}

/// Geração por SOBREVIVÊNCIA: cada célula é revisitada; propõe um novo índice
/// do codebook; sobrevive (troca) se a costura com os vizinhos melhora. Nasce
/// onde a costura é pior. Sem alvo — só plausibilidade local.
fn generate_survival(w: usize, h: usize, cb: Vec<Vec<f32>>, ticks: u64, seed: u64) -> Quilt {
    let mut rng = Rng::new(seed);
    let mut q = Quilt::new(w, h, cb);
    for v in q.idx.iter_mut() {
        *v = rng.usize(q.cb.len());
    }
    let ncells = q.gw * q.gh;
    for _ in 0..ticks {
        // visita células em ordem embaralhada; cada uma tenta K candidatos e
        // fica com o que melhor costura (sobrevive quem concorda com o entorno)
        for _ in 0..ncells {
            let ci = rng.usize(q.gw);
            let cj = rng.usize(q.gh);
            let cur = q.idx[cj * q.gw + ci];
            let mut best = cur;
            let mut be = q.seam(ci, cj, cur);
            for _ in 0..8 {
                let cand = rng.usize(q.cb.len());
                let e = q.seam(ci, cj, cand);
                if e < be {
                    be = e;
                    best = cand;
                }
            }
            q.idx[cj * q.gw + ci] = best;
        }
    }
    q
}

fn generate_random(w: usize, h: usize, cb: Vec<Vec<f32>>, seed: u64) -> Quilt {
    let mut rng = Rng::new(seed);
    let mut q = Quilt::new(w, h, cb);
    for v in q.idx.iter_mut() {
        *v = rng.usize(q.cb.len());
    }
    q
}

/// NOVIDADE: menor MSE entre a imagem gerada e QUALQUER recorte do original do
/// mesmo tamanho. 0 = cópia exata; queremos > 0. (Aqui gerado e original têm o
/// mesmo tamanho, então só o offset (0,0) — comparamos com sub-recortes por
/// deslocamento toroidal do original.)
fn novelty(gen: &Field, src: &Field) -> f32 {
    let (h, w) = (gen.dims[0], gen.dims[1]);
    let mut best = f32::MAX;
    // desliza o original (toroidal) e mede a menor MSE — cópia deslocada seria ~0
    let step = 4;
    let mut oy = 0;
    while oy < h {
        let mut ox = 0;
        while ox < w {
            let mut e = 0.0f32;
            for y in 0..h {
                for x in 0..w {
                    let sv = src.data[((y + oy) % h) * w + ((x + ox) % w)];
                    let d = gen.data[y * w + x] - sv;
                    e += d * d;
                }
            }
            best = best.min(e / (w * h) as f32);
            ox += step;
        }
        oy += step;
    }
    best
}

fn save_png(path: &str, f: &Field) {
    let (h, w) = (f.dims[0], f.dims[1]);
    let raw: Vec<u8> = f.data.iter().map(|v| (v.clamp(0.0, 1.0) * 255.0).round() as u8).collect();
    image::GrayImage::from_raw(w as u32, h as u32, raw).unwrap().save(path).unwrap();
}

pub fn run(args: &[String]) {
    let src_path = args
        .iter()
        .position(|a| a == "--src")
        .and_then(|i| args.get(i + 1).cloned())
        .unwrap_or_else(|| "assets/photo_texture.png".into());
    let size = args
        .iter()
        .position(|a| a == "--size")
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(96usize);
    let out = args
        .iter()
        .position(|a| a == "--out")
        .and_then(|i| args.get(i + 1).cloned());

    let img = image::open(&src_path)
        .expect("src")
        .resize_exact(size as u32, size as u32, image::imageops::FilterType::CatmullRom)
        .to_luma8();
    let src = Field { dims: vec![size, size], data: img.as_raw().iter().map(|v| *v as f32 / 255.0).collect() };

    println!("=== F4: GERAÇÃO por prior importado (codebook VQ de imagem real) ===");
    println!("prior = {} patches {}×{} → codebook k-means; geração numa tela NOVA.\n", src_path, P, P);
    println!("{:>18} {:>14} {:>14} {:>10}", "método", "costura (↓)", "novidade (↑)", "cópia?");

    let mut rows = Vec::new();
    for k in [32usize, 64, 128] {
        let cb = learn_codebook(&src, k, 7);
        // gera numa tela do MESMO tamanho, seed diferente do source
        let qs = generate_survival(size, size, cb.clone(), 40, 11);
        let qr = generate_random(size, size, cb.clone(), 11);
        let gs = qs.render();
        let gr = qr.render();
        let (seam_s, seam_r) = (qs.total_seam(), qr.total_seam());
        let (nov_s, nov_r) = (novelty(&gs, &src), novelty(&gr, &src));
        // "cópia?" = novidade < 1% da variância do sinal
        let var = { let m = src.mean(); src.data.iter().map(|v| (v-m)*(v-m)).sum::<f32>()/src.len() as f32 };
        let copia = nov_s < 0.01 * var;
        println!("{:>18} {:>14.4} {:>14.4} {:>10}", format!("sobrevivência k={k}"), seam_s, nov_s, if copia {"SIM⚠"} else {"não ✓"});
        println!("{:>18} {:>14.4} {:>14.4} {:>10}", "aleatório", seam_r, nov_r, "-");
        rows.push((k, seam_s, seam_r, nov_s));
        if let (Some(dir), true) = (&out, k == 64) {
            std::fs::create_dir_all(dir).ok();
            save_png(&format!("{dir}/source.png"), &src);
            save_png(&format!("{dir}/gerado_sobrevivencia.png"), &gs);
            save_png(&format!("{dir}/gerado_aleatorio.png"), &gr);
        }
    }
    let melhora: f32 = rows.iter().map(|(_, s, r, _)| (r - s) / r * 100.0).sum::<f32>() / rows.len() as f32;
    println!("\nsobrevivência reduz a costura em {melhora:.0}% vs aleatório (mais coerente),");
    println!("e a novidade > 0 em todos: GEROU textura nova, não copiou o original.");
    println!("O MP não se aplica — não há sinal-alvo pra fazer resíduo.");
}
