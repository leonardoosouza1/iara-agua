//! `agua-lab polyglot` — a PONTE entre representações (tese IARA-Router).
//!
//! A pergunta do Leonardo: e se cada "língua" é boa em coisas diferentes?
//! bytes bons em X, tokens em Y, CfC em Z... Um tradutor que trabalha em várias.
//!
//! O seletor por sobrevivência JÁ É um roteador: átomos competem pelo resíduo.
//! Aqui os candidatos vêm de FAMÍLIAS diferentes (representações), e cada região
//! do sinal é reivindicada por quem a explica mais barato. O "tradutor" não é um
//! modelo que converte — é uma ARENA onde as línguas competem por região; o
//! roteamento EMERGE.
//!
//! 3 famílias (línguas), cada uma boa numa coisa:
//!   GABOR  — oscilação / textura periódica
//!   BLOB   — regiões lisas (gaussiana, freq 0)
//!   BORDA  — descontinuidade / meia-face com transição suave (o que Gabor
//!            aproxima mal: uma soma de senos ringa numa borda dura)
//!
//! Predição (com regime de validade honesto): num sinal HETEROGÊNEO (regiões de
//! cada tipo), o dicionário poliglota bate cada língua sozinha no mesmo orçamento
//! de átomos, E cada região é reivindicada pela língua certa (o roteador emerge).
//! Num sinal HOMOGÊNEO, uma língua vence e misturar só paga o imposto ln(M) da
//! Lei do Garimpo à toa. Liga com o M74 do resonance_forest (generalista > swarm
//! quando há estrutura comum) — aqui o oposto: línguas COMPLEMENTARES, não iguais.

use agua_core::{psnr, Field, Rng};

#[derive(Clone, Copy, PartialEq)]
enum Lang {
    Gabor,
    Blob,
    Edge,
}

fn lang_name(l: Lang) -> &'static str {
    match l {
        Lang::Gabor => "gabor",
        Lang::Blob => "blob",
        Lang::Edge => "borda",
    }
}

/// Renderiza um átomo unitário de uma língua num campo denso (patch completo).
struct Atom {
    lang: Lang,
    cx: f32,
    cy: f32,
    sigma: f32,
    // gabor
    freq: f32,
    theta: f32,
    phase: f32,
}

impl Atom {
    fn sample(lang: Lang, w: usize, h: usize, rng: &mut Rng) -> Atom {
        let cx = rng.range(0.0, w as f32 - 1.0);
        let cy = rng.range(0.0, h as f32 - 1.0);
        let sigma = rng.range(2.0, (w as f32 * 0.25).max(4.0));
        Atom {
            lang,
            cx,
            cy,
            sigma,
            freq: rng.range(0.03, 0.4),
            theta: rng.range(0.0, std::f32::consts::PI),
            phase: rng.range(0.0, std::f32::consts::TAU),
        }
    }

    /// campo denso do átomo unitário (amp=1)
    fn render(&self, w: usize, h: usize) -> Vec<f32> {
        let mut v = vec![0.0f32; w * h];
        let r = (3.0 * self.sigma) as i64;
        let x0 = (self.cx as i64 - r).max(0) as usize;
        let x1 = ((self.cx as i64 + r + 1).min(w as i64)) as usize;
        let y0 = (self.cy as i64 - r).max(0) as usize;
        let y1 = ((self.cy as i64 + r + 1).min(h as i64)) as usize;
        let inv = 1.0 / (2.0 * self.sigma * self.sigma);
        let (st, ct) = self.theta.sin_cos();
        let k = std::f32::consts::TAU * self.freq;
        for y in y0..y1 {
            let dy = y as f32 - self.cy;
            for x in x0..x1 {
                let dx = x as f32 - self.cx;
                let env = (-(dx * dx + dy * dy) * inv).exp();
                let u = dx * ct + dy * st;
                v[y * w + x] = match self.lang {
                    Lang::Gabor => env * (k * u + self.phase).cos(),
                    Lang::Blob => env, // gaussiana lisa
                    // borda: transição suave (tanh) projetada em θ, janelada
                    Lang::Edge => env * (u / self.sigma).tanh(),
                };
            }
        }
        v
    }
}

fn dot(a: &[f32], b: &[f32]) -> f32 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Matching pursuit POLIGLOTA: a cada passo amostra C candidatos de CADA língua
/// permitida, escolhe o globalmente melhor (forma fechada de amplitude), subtrai.
/// Devolve (psnr, uso por língua, mapa de língua por pixel dominante).
fn mp_polyglot(
    target: &Field,
    langs: &[Lang],
    n_atoms: usize,
    cands: usize,
    seed: u64,
) -> (f32, [usize; 3], Vec<i8>) {
    let (h, w) = (target.dims[0], target.dims[1]);
    let dc = target.mean();
    let mut resid: Vec<f32> = target.data.iter().map(|v| v - dc).collect();
    let mut rng = Rng::new(seed);
    let mut usage = [0usize; 3];
    let mut owner = vec![-1i8; w * h]; // qual língua "explicou" cada pixel

    for _ in 0..n_atoms {
        let mut best: Option<(f32, Vec<f32>, f32, Lang)> = None; // (explicado, patch, amp, lang)
        for &lang in langs {
            for _ in 0..cands {
                let a = Atom::sample(lang, w, h, &mut rng);
                let patch = a.render(w, h);
                let n2 = dot(&patch, &patch);
                if n2 < 1e-6 {
                    continue;
                }
                let amp = dot(&patch, &resid) / n2;
                let explained = amp * amp * n2;
                if best.as_ref().map_or(true, |(e, _, _, _)| explained > *e) {
                    best = Some((explained, patch, amp, lang));
                }
            }
        }
        if let Some((_, patch, amp, lang)) = best {
            let li = match lang {
                Lang::Gabor => 0,
                Lang::Blob => 1,
                Lang::Edge => 2,
            };
            usage[li] += 1;
            for (i, p) in patch.iter().enumerate() {
                resid[i] -= amp * p;
                if p.abs() * amp.abs() > 0.05 {
                    owner[i] = li as i8; // último a explicar forte esse pixel
                }
            }
        }
    }
    let mse = resid.iter().map(|r| r * r).sum::<f32>() / resid.len() as f32;
    (psnr(mse, 1.0), usage, owner)
}

/// Alvo HETEROGÊNEO: 3 faixas verticais — lisa (gradiente), textura (listras),
/// borda dura (dois patamares). Cada faixa tem a língua "certa".
fn heterogeneous(size: usize) -> Field {
    let s = size as f32;
    Field::from_fn_2d(size, size, |x, y| {
        let (xf, yf) = (x as f32, y as f32);
        if xf < s / 3.0 {
            // liso: gradiente suave → BLOB
            0.3 + 0.4 * (yf / s)
        } else if xf < 2.0 * s / 3.0 {
            // textura: listras → GABOR
            0.5 + 0.35 * (std::f32::consts::TAU * 0.15 * yf).sin()
        } else {
            // borda dura: dois patamares → BORDA
            if yf < s / 2.0 {
                0.15
            } else {
                0.85
            }
        }
    })
}

/// Alvo HOMOGÊNEO: só textura (para mostrar o regime onde misturar NÃO paga).
fn homogeneous(size: usize) -> Field {
    Field::from_fn_2d(size, size, |x, y| {
        0.5 + 0.35 * (std::f32::consts::TAU * 0.15 * (x + y) as f32).sin()
    })
}

pub fn run(args: &[String]) {
    let size = args
        .iter()
        .position(|a| a == "--size")
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(48usize);
    let n = 60usize;
    let cands = 20usize;

    println!("=== PONTE ENTRE LÍNGUAS — dicionário poliglota vs cada língua sozinha ===");
    println!("cada família é boa numa coisa: gabor=textura, blob=liso, borda=descontinuidade.");
    println!("o roteamento por região EMERGE da competição (Lei do Garimpo dá o custo).\n");

    for (nome, target) in [("HETEROGÊNEO (3 regiões)", heterogeneous(size)), ("HOMOGÊNEO (só textura)", homogeneous(size))] {
        println!("── alvo {nome} ──");
        println!("{:>22} {:>10} {:>28}", "dicionário", "PSNR", "uso por língua");
        // cada língua sozinha (mesmo orçamento TOTAL de candidatos por passo:
        // sozinha recebe 3×cands para igualar o poliglota que testa 3 famílias)
        for &l in &[Lang::Gabor, Lang::Blob, Lang::Edge] {
            let mut best = f32::MIN;
            for seed in [11u64, 22, 33] {
                let (p, _, _) = mp_polyglot(&target, &[l], n, cands * 3, seed);
                best = best.max(p);
            }
            println!("{:>22} {:>10.2} {:>28}", format!("só {}", lang_name(l)), best, "-");
        }
        // poliglota
        let mut pg = f32::MIN;
        let mut best_usage = [0; 3];
        for seed in [11u64, 22, 33] {
            let (p, u, _) = mp_polyglot(&target, &[Lang::Gabor, Lang::Blob, Lang::Edge], n, cands, seed);
            if p > pg {
                pg = p;
                best_usage = u;
            }
        }
        println!(
            "{:>22} {:>10.2} {:>28}",
            "POLIGLOTA (3 línguas)",
            pg,
            format!("g{} b{} e{}", best_usage[0], best_usage[1], best_usage[2])
        );

        // roteamento por região: no heterogêneo, qual língua domina cada faixa?
        if nome.starts_with("HETERO") {
            let (_, _, owner) = mp_polyglot(&target, &[Lang::Gabor, Lang::Blob, Lang::Edge], n, cands, 11);
            let (h, w) = (size, size);
            let mut faixa = [[0usize; 3]; 3]; // [faixa][língua]
            for y in 0..h {
                for x in 0..w {
                    let f = (x * 3 / w).min(2);
                    let o = owner[y * w + x];
                    if o >= 0 {
                        faixa[f][o as usize] += 1;
                    }
                }
            }
            // faixa 0=liso→blob(idx1), 1=textura→gabor(idx0), 2=degrau→borda(idx2)
            let esperado_idx = [1usize, 0, 2];
            let esperado_nome = ["blob (liso)", "gabor (textura)", "borda (degrau)"];
            println!("\n  roteamento EMERGENTE (língua dominante por faixa):");
            for fi in 0..3 {
                let dom = (0..3).max_by_key(|&l| faixa[fi][l]).unwrap();
                let ok = if dom == esperado_idx[fi] { "✓" } else { "✗" };
                println!(
                    "    faixa {}: dominante = {:<6} (esperado {}) {ok}   g{} b{} e{}",
                    fi + 1,
                    lang_name([Lang::Gabor, Lang::Blob, Lang::Edge][dom]),
                    esperado_nome[fi],
                    faixa[fi][0], faixa[fi][1], faixa[fi][2]
                );
            }

            // PSNR POR REGIÃO: poliglota vs gabor-sozinha (o valor é LOCAL)
            println!("\n  PSNR por faixa — poliglota vs melhor-única (gabor):");
            let region_psnr = |langs: &[Lang], seed: u64| -> [f32; 3] {
                let (h, w) = (size, size);
                let (_, _, _) = mp_polyglot(&target, langs, n, cands, seed);
                // recomputa resíduo p/ PSNR por faixa
                let dc = target.mean();
                let mut resid: Vec<f32> = target.data.iter().map(|v| v - dc).collect();
                let mut rng = Rng::new(seed);
                for _ in 0..n {
                    let mut best: Option<(f32, Vec<f32>, f32)> = None;
                    for &lang in langs {
                        for _ in 0..cands {
                            let a = Atom::sample(lang, w, h, &mut rng);
                            let patch = a.render(w, h);
                            let n2 = dot(&patch, &patch);
                            if n2 < 1e-6 { continue; }
                            let amp = dot(&patch, &resid)/n2;
                            let ex = amp*amp*n2;
                            if best.as_ref().map_or(true, |(e,_,_)| ex>*e) { best=Some((ex,patch,amp)); }
                        }
                    }
                    if let Some((_,patch,amp))=best { for (i,pp) in patch.iter().enumerate(){resid[i]-=amp*pp;} }
                }
                let mut se=[0.0f32;3]; let mut cnt=[0usize;3];
                for y in 0..h { for x in 0..w { let f=(x*3/w).min(2); se[f]+=resid[y*w+x].powi(2); cnt[f]+=1; } }
                [psnr(se[0]/cnt[0] as f32,1.0), psnr(se[1]/cnt[1] as f32,1.0), psnr(se[2]/cnt[2] as f32,1.0)]
            };
            let pg = region_psnr(&[Lang::Gabor, Lang::Blob, Lang::Edge], 11);
            let gb = region_psnr(&[Lang::Gabor], 11);
            for (fi, nm) in ["liso","textura","degrau"].iter().enumerate() {
                let d = pg[fi]-gb[fi];
                println!("    faixa {} ({}): poliglota {:.2} vs gabor {:.2}  → {:+.2} dB {}",
                    fi+1, nm, pg[fi], gb[fi], d, if d>0.3 {"◄ poliglota ganha"} else if d < -0.3 {"gabor ganha"} else {""});
            }
        }
        println!();
    }
    println!("Lei do Garimpo: o poliglota paga 2σ²·ln(3·C) de busca. Só vale se a");
    println!("diversidade das regiões comprar mais que o imposto — regime de validade.");
}
