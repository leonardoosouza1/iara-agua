//! agua-metrics — baselines HONESTOS pra comparação obrigatória.
//!
//! Regra do protocolo: nenhum resultado de sobrevivência é reportado sem
//! (a) random-atoms (piso: mesma família, geometria aleatória),
//! (b) matching pursuit (teto guloso: melhor-de-C candidatos por passo),
//! (c) JPEG (referência de mundo real em bytes×PSNR).
//! Sobrevivência que não bate (a) é ruído; distância pra (b) mede o custo
//! de não ter gradiente/ganância global.

use agua_atoms::{fit_amp, AtomBounds, Gabor2D};
use agua_core::{psnr, sample_hotspot, Field2D, Rng};
use serde::Serialize;

/// Ponto taxa-distorção de um método.
#[derive(Clone, Debug, Serialize)]
pub struct RdPoint {
    pub label: String,
    pub atoms: usize,
    pub bytes: usize,
    pub psnr_db: f32,
}

fn atom_bytes(n: usize) -> usize {
    n * Gabor2D::PARAMS * 4
}

/// Piso: N átomos com geometria aleatória (posição uniforme), amp em forma
/// fechada, adicionados sequencialmente. Curva em checkpoints.
pub fn random_atoms(
    target: &Field2D,
    n: usize,
    bounds: &AtomBounds,
    seed: u64,
    checkpoints: &[usize],
) -> Vec<RdPoint> {
    let mut rng = Rng::new(seed);
    let mut residual = target.clone();
    let mut out = Vec::new();
    let mut placed = 0usize;
    for k in 1..=n {
        let x = rng.range(0.0, (target.w - 1) as f32);
        let y = rng.range(0.0, (target.h - 1) as f32);
        let g = Gabor2D::rain_at(x, y, bounds, &mut rng);
        let unit = g.render_unit(target.w, target.h);
        if let Some((amp, _)) = fit_amp(&unit, &residual) {
            residual.add_scaled_patch(&unit, -amp);
            placed += 1;
        }
        if checkpoints.contains(&k) {
            let mse = residual.energy() / (target.w * target.h) as f32;
            out.push(RdPoint {
                label: "random".into(),
                atoms: placed,
                bytes: atom_bytes(placed),
                psnr_db: psnr(mse, 1.0),
            });
        }
    }
    out
}

/// Teto guloso (matching pursuit estocástico): a cada passo amostra C
/// candidatos em hotspots do resíduo e adiciona o que mais explica.
pub fn matching_pursuit(
    target: &Field2D,
    n: usize,
    candidates: usize,
    bounds: &AtomBounds,
    seed: u64,
    checkpoints: &[usize],
) -> Vec<RdPoint> {
    let mut rng = Rng::new(seed);
    let mut residual = target.clone();
    let mut out = Vec::new();
    for k in 1..=n {
        let mut best: Option<(f32, f32, agua_core::Patch)> = None; // (explicado, amp, unit)
        for _ in 0..candidates {
            let (x, y) = sample_hotspot(&residual, &mut rng, 8);
            let g = Gabor2D::rain_at(x as f32, y as f32, bounds, &mut rng);
            let unit = g.render_unit(target.w, target.h);
            if let Some((amp, explained)) = fit_amp(&unit, &residual) {
                if best.as_ref().map_or(true, |(e, _, _)| explained > *e) {
                    best = Some((explained, amp, unit));
                }
            }
        }
        if let Some((_, amp, unit)) = best {
            residual.add_scaled_patch(&unit, -amp);
        }
        if checkpoints.contains(&k) {
            let mse = residual.energy() / (target.w * target.h) as f32;
            out.push(RdPoint {
                label: "matching_pursuit".into(),
                atoms: k,
                bytes: atom_bytes(k),
                psnr_db: psnr(mse, 1.0),
            });
        }
    }
    out
}

/// Referência mundo-real: JPEG em várias qualidades (bytes reais do arquivo).
/// `target` em [0,1] (com DC), vira Luma8, roundtrip encode→decode→PSNR.
pub fn jpeg_curve(target: &Field2D, qualities: &[u8]) -> Vec<RdPoint> {
    let (w, h) = (target.w as u32, target.h as u32);
    let raw: Vec<u8> = target
        .data
        .iter()
        .map(|v| (v.clamp(0.0, 1.0) * 255.0).round() as u8)
        .collect();
    let img = image::GrayImage::from_raw(w, h, raw.clone()).expect("dims");
    let mut out = Vec::new();
    for &q in qualities {
        let mut buf = Vec::new();
        {
            let mut enc = image::codecs::jpeg::JpegEncoder::new_with_quality(&mut buf, q);
            if enc
                .encode(img.as_raw(), w, h, image::ExtendedColorType::L8)
                .is_err()
            {
                continue;
            }
        }
        let Ok(decoded) = image::load_from_memory(&buf) else {
            continue;
        };
        let dec = decoded.to_luma8();
        let mse: f32 = raw
            .iter()
            .zip(dec.as_raw())
            .map(|(a, b)| {
                let d = (*a as f32 - *b as f32) / 255.0;
                d * d
            })
            .sum::<f32>()
            / raw.len() as f32;
        out.push(RdPoint {
            label: format!("jpeg_q{q}"),
            atoms: 0,
            bytes: buf.len(),
            psnr_db: psnr(mse, 1.0),
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn target() -> Field2D {
        let mut f = Field2D::from_fn(48, 48, |x, y| {
            0.5 + 0.4 * (std::f32::consts::TAU * 0.05 * x as f32).sin()
                * (std::f32::consts::TAU * 0.03 * y as f32).cos()
        });
        let m = f.mean();
        for v in f.data.iter_mut() {
            *v -= m;
        }
        f
    }

    fn bounds() -> AtomBounds {
        AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 16.0 }
    }

    #[test]
    fn mp_monotonic_improvement() {
        let t = target();
        let pts = matching_pursuit(&t, 40, 12, &bounds(), 3, &[10, 20, 40]);
        assert_eq!(pts.len(), 3);
        assert!(pts[1].psnr_db >= pts[0].psnr_db - 0.01);
        assert!(pts[2].psnr_db >= pts[1].psnr_db - 0.01);
        assert!(pts[2].psnr_db > 15.0, "MP fraco demais: {:.1}", pts[2].psnr_db);
    }

    #[test]
    fn mp_beats_random() {
        let t = target();
        let mp = matching_pursuit(&t, 30, 12, &bounds(), 5, &[30]);
        let rnd = random_atoms(&t, 30, &bounds(), 5, &[30]);
        assert!(mp[0].psnr_db > rnd[0].psnr_db, "MP {:.1} <= random {:.1}", mp[0].psnr_db, rnd[0].psnr_db);
    }

    #[test]
    fn jpeg_roundtrip_sane() {
        // JPEG precisa do DC — usa alvo em [0,1]
        let f = Field2D::from_fn(48, 48, |x, y| {
            0.5 + 0.3 * (std::f32::consts::TAU * 0.04 * (x + y) as f32).sin()
        });
        let pts = jpeg_curve(&f, &[30, 80]);
        assert_eq!(pts.len(), 2);
        assert!(pts.iter().all(|p| p.psnr_db > 20.0 && p.bytes > 0));
        assert!(pts[1].psnr_db >= pts[0].psnr_db);
    }
}
