//! `agua-lab compare` — gera os artefatos VISUAIS de um experimento de imagem:
//! alvo, reconstrução da sobrevivência, do MP iso-budget, do random; o resíduo;
//! a decomposição da PIRÂMIDE (recon por banda de frequência); e a evolução
//! temporal (recon a cada checkpoint). Tudo PNG no diretório de saída.

use agua_atoms::{AtomBounds, Gabor2D};
use agua_core::Field;
use agua_life::{LifeConfig, World};
use agua_metrics::{matching_pursuit_recon, random_atoms_recon};
use std::path::Path;

/// Fatias da pirâmide (ciclos/px): grosso, médio, fino.
const BANDS: [(f32, f32, &str); 3] =
    [(0.0, 0.06, "grosso"), (0.06, 0.12, "medio"), (0.12, 10.0, "fino")];

fn png(path: &Path, f: &Field, lo: f32, hi: f32) {
    let (h, w) = (f.dims[0], f.dims[1]);
    let span = (hi - lo).max(1e-6);
    let raw: Vec<u8> = f
        .data
        .iter()
        .map(|v| (((v - lo) / span).clamp(0.0, 1.0) * 255.0).round() as u8)
        .collect();
    image::GrayImage::from_raw(w as u32, h as u32, raw)
        .expect("dims")
        .save(path)
        .expect("png");
}

pub struct CompareArgs<'a> {
    pub target_dc: Field,
    pub life: LifeConfig,
    pub bounds: AtomBounds,
    pub ticks: u64,
    pub seed: u64,
    pub out: &'a Path,
    /// ticks em que salva a evolução
    pub checkpoints: &'a [u64],
}

pub fn run(args: CompareArgs) {
    std::fs::create_dir_all(args.out).expect("mkdir");
    let dc = args.target_dc.mean();
    let mut target = args.target_dc.clone();
    for v in target.data.iter_mut() {
        *v -= dc;
    }
    let with_dc = |f: &Field| {
        let mut o = f.clone();
        for v in o.data.iter_mut() {
            *v += dc;
        }
        o
    };

    png(&args.out.join("01_alvo.png"), &args.target_dc, 0.0, 1.0);

    let mut world: World<Gabor2D> =
        World::new(target.clone(), args.life, args.bounds.clone(), args.seed);

    for tick in 1..=args.ticks {
        world.step();
        if args.checkpoints.contains(&tick) {
            png(
                &args.out.join(format!("evo_t{tick:05}.png")),
                &with_dc(&world.recon),
                0.0,
                1.0,
            );
        }
    }
    world.rebuild();

    let n = world.beings.len().max(1);
    let psnr_life = world.psnr(1.0);
    png(&args.out.join("02_vida.png"), &with_dc(&world.recon), 0.0, 1.0);
    // resíduo: 0 = cinza, ±0.35 = preto/branco
    png(&args.out.join("03_residuo_vida.png"), &world.residual, -0.35, 0.35);

    // pirâmide: cada fatia soma amp·ψ só dos seres daquela banda
    for (lo, hi, name) in BANDS {
        let band = world.recon_band(lo, hi);
        let count = world
            .beings
            .iter()
            .filter(|b| {
                let f = b.atom.freq;
                f >= lo && f < hi
            })
            .count();
        png(&args.out.join(format!("pir_{name}.png")), &band, -0.35, 0.35);
        println!("  pirâmide[{name}]: {count} seres");
    }

    // baselines com o MESMO N e orçamento iso (evals/n candidatos)
    let cands = ((world.evals as usize) / n).max(1);
    let mp = matching_pursuit_recon::<Gabor2D>(&target, n, cands, &args.bounds, args.seed ^ 0x1234);
    let rnd = random_atoms_recon::<Gabor2D>(&target, n, &args.bounds, args.seed ^ 0xABCD);
    let mse = |r: &Field| {
        target
            .data
            .iter()
            .zip(&r.data)
            .map(|(t, c)| (t - c) * (t - c))
            .sum::<f32>()
            / target.len() as f32
    };
    png(&args.out.join("04_mp.png"), &with_dc(&mp), 0.0, 1.0);
    png(&args.out.join("05_random.png"), &with_dc(&rnd), 0.0, 1.0);

    println!(
        "seres={n} evals={} mp_cands={cands}\n  vida   {:.2} dB\n  MP iso {:.2} dB\n  random {:.2} dB",
        world.evals,
        psnr_life,
        agua_core::psnr(mse(&mp), 1.0),
        agua_core::psnr(mse(&rnd), 1.0)
    );
}
