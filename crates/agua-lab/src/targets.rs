//! Alvos sintéticos multi-escala determinísticos por modalidade.
//! Todos têm componentes grossos E finos — testam a claim "a pirâmide emerge".

use agua_core::Field;
use std::f32::consts::TAU;

/// Imagem: gradiente suave (baixa freq) + blob gaussiano (média) + listras
/// finas num quadrante (alta freq) + quadrado de borda dura. Valores [0,1].
pub fn synthetic_image(size: usize) -> Field {
    let s = size as f32;
    Field::from_fn_2d(size, size, |x, y| {
        let (xf, yf) = (x as f32, y as f32);
        let grad = 0.25 * (TAU * xf / s * 0.5).sin() * (TAU * yf / s * 0.5).cos();
        let blob = 0.45
            * (-((xf - s * 0.3).powi(2) + (yf - s * 0.32).powi(2)) / (2.0 * (s * 0.12).powi(2)))
                .exp();
        let stripes = if xf > s * 0.55 && yf > s * 0.55 {
            0.22 * (TAU * 0.22 * xf).sin()
        } else {
            0.0
        };
        let square = if (xf - s * 0.7).abs() < s * 0.09 && (yf - s * 0.25).abs() < s * 0.09 {
            0.3
        } else {
            0.0
        };
        (0.45 + grad + blob + stripes + square).clamp(0.0, 1.0)
    })
}

/// Áudio: tom grave contínuo (baixa freq) + burst médio no meio + tom agudo
/// no último quarto + transiente (click). Valores [-1,1], média ~0.
pub fn synthetic_audio(n: usize) -> Field {
    let nf = n as f32;
    Field::from_fn_1d(n, |i| {
        let t = i as f32;
        let low = 0.35 * (TAU * 0.004 * t).sin();
        let mid = if t > nf * 0.28 && t < nf * 0.52 {
            let env = ((t - nf * 0.28) / (nf * 0.04)).min(1.0) * ((nf * 0.52 - t) / (nf * 0.04)).min(1.0);
            0.3 * env * (TAU * 0.035 * t).sin()
        } else {
            0.0
        };
        let high = if t > nf * 0.75 {
            0.2 * (TAU * 0.18 * t).sin()
        } else {
            0.0
        };
        let click = if t >= nf * 0.6 {
            0.5 * (-(t - nf * 0.6) / 6.0).exp() * (TAU * 0.3 * t).cos()
        } else {
            0.0
        };
        (low + mid + high + click).clamp(-1.0, 1.0)
    })
}

/// Vídeo: blob gaussiano que atravessa a cena (movimento = estrutura
/// temporal) + listras estáticas num canto + pulso global lento. [0,1].
pub fn synthetic_video(frames: usize, size: usize) -> Field {
    let s = size as f32;
    let tf = frames.max(2) as f32;
    Field::from_fn_3d(frames, size, size, |x, y, t| {
        let (xf, yf, ti) = (x as f32, y as f32, t as f32);
        let prog = ti / (tf - 1.0);
        let (bx, by) = (s * (0.2 + 0.6 * prog), s * (0.35 + 0.25 * (TAU * prog).sin() * 0.5));
        let blob = 0.5
            * (-((xf - bx).powi(2) + (yf - by).powi(2)) / (2.0 * (s * 0.1).powi(2))).exp();
        let stripes = if xf > s * 0.6 && yf > s * 0.6 {
            0.2 * (TAU * 0.2 * xf).sin()
        } else {
            0.0
        };
        let pulse = 0.1 * (TAU * 0.08 * ti).sin();
        (0.4 + blob + stripes + pulse).clamp(0.0, 1.0)
    })
}
