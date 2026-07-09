//! agua-atoms — os "seres/alfinetes": funções-base locais paramétricas.
//!
//! Modalidade vive AQUI, não na física de vida:
//!   `Gabor1D` — grão de áudio (envelope gaussiano × senoide)
//!   `Gabor2D` — átomo de imagem
//!   `Gabor3D` — átomo espaço-temporal de vídeo
//! Um Gabor com freq→0 degenera num blob gaussiano, então a mesma família
//! cobre do grosso (baixa freq) ao fino (alta freq) — a pirâmide pode emergir.
//!
//! A AMPLITUDE nunca é evoluída: é ajustada em forma fechada contra o
//! resíduo (amp* = ⟨r,ψ⟩/||ψ||²) — é o ser "bebendo" o resíduo. A GEOMETRIA
//! (posição/freq/orientação/fase/escala) é o que a sobrevivência esculpe.

use agua_core::{Field, Patch, Rng, Run};
use serde::{Deserialize, Serialize};

/// Limites do espaço de parâmetros (por experimento, vem do TOML).
/// Frequências em ciclos/amostra (ou ciclos/pixel, ciclos/frame).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AtomBounds {
    pub freq_min: f32,
    pub freq_max: f32,
    pub sigma_min: f32,
    pub sigma_max: f32,
}

/// Contrato do ser-átomo. A física (`agua-life`) só conhece isto —
/// implementar este trait = nova modalidade, física intocada.
pub trait Atom: Clone + Serialize {
    /// nº de parâmetros da codificação (geometria + 1 amplitude)
    const PARAMS: usize;
    /// Nasce da chuva numa posição (coords na ordem de `dims`).
    fn rain_at(pos: &[f32], b: &AtomBounds, rng: &mut Rng) -> Self;
    /// Filho/candidato mutado (clamps dentro de `dims` e `bounds`).
    fn mutate(&self, rng: &mut Rng, s: f32, b: &AtomBounds, dims: &[usize]) -> Self;
    /// Renderiza o patch UNITÁRIO (amp=1) recortado ao campo.
    fn render_unit(&self, dims: &[usize]) -> Patch;
    /// Frequência característica (pro histograma da pirâmide).
    fn freq(&self) -> f32;
}

/// Ajuste de amplitude em forma fechada: amp* = ⟨residual,ψ⟩/||ψ||².
/// Retorna (amp, energia_explicada = amp²·||ψ||²). `None` se suporte degenerado.
pub fn fit_amp(unit: &Patch, residual: &Field) -> Option<(f32, f32)> {
    let n2 = unit.norm2();
    if n2 < 1e-6 {
        return None;
    }
    let d = unit.dot_field(residual);
    let amp = d / n2;
    if !amp.is_finite() {
        return None;
    }
    Some((amp, amp * amp * n2))
}

#[inline]
fn sigma_from_freq(freq: f32, b: &AtomBounds, rng: &mut Rng) -> f32 {
    let cycles = rng.range(0.35, 0.9);
    (cycles / freq).clamp(b.sigma_min, b.sigma_max)
}

// ─────────────────────────── 1D — áudio ───────────────────────────

/// Grão de Gabor 1D: envelope gaussiano × senoide. dims = [n_amostras].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Gabor1D {
    pub c: f32,
    /// ciclos por amostra
    pub freq: f32,
    pub phase: f32,
    pub sigma: f32,
}

impl Atom for Gabor1D {
    const PARAMS: usize = 5;

    fn rain_at(pos: &[f32], b: &AtomBounds, rng: &mut Rng) -> Self {
        let freq = rng.log_range(b.freq_min, b.freq_max);
        Gabor1D {
            c: pos[0],
            freq,
            phase: rng.range(0.0, std::f32::consts::TAU),
            sigma: sigma_from_freq(freq, b, rng),
        }
    }

    fn mutate(&self, rng: &mut Rng, s: f32, b: &AtomBounds, dims: &[usize]) -> Self {
        let n = dims[0] as f32;
        Gabor1D {
            c: (self.c + rng.normal() * s * self.sigma).clamp(0.0, n - 1.0),
            freq: (self.freq * (rng.normal() * s * 0.5).exp()).clamp(b.freq_min, b.freq_max),
            phase: self.phase + rng.normal() * s * 1.2,
            sigma: (self.sigma * (rng.normal() * s * 0.5).exp()).clamp(b.sigma_min, b.sigma_max),
        }
    }

    fn render_unit(&self, dims: &[usize]) -> Patch {
        let n = dims[0];
        let r = (2.5 * self.sigma).ceil().max(1.0) as i64;
        let x0 = (self.c as i64 - r).max(0) as usize;
        let x1 = ((self.c as i64 + r + 1).min(n as i64)) as usize;
        let len = x1.saturating_sub(x0);
        let inv2s2 = 1.0 / (2.0 * self.sigma * self.sigma);
        let k = std::f32::consts::TAU * self.freq;
        let data: Vec<f32> = (0..len)
            .map(|i| {
                let dx = (x0 + i) as f32 - self.c;
                (-(dx * dx) * inv2s2).exp() * (k * dx + self.phase).cos()
            })
            .collect();
        Patch { runs: vec![Run { field_start: x0, len }], data }
    }

    fn freq(&self) -> f32 {
        self.freq
    }
}

// ─────────────────────────── 2D — imagem ───────────────────────────

/// Átomo Gabor 2D. dims = [h, w]; pos = [y, x]. `amp` vive fora (no ser).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Gabor2D {
    pub cx: f32,
    pub cy: f32,
    /// ciclos por pixel
    pub freq: f32,
    pub theta: f32,
    pub phase: f32,
    pub sigma: f32,
}

impl Atom for Gabor2D {
    const PARAMS: usize = 7;

    fn rain_at(pos: &[f32], b: &AtomBounds, rng: &mut Rng) -> Self {
        let freq = rng.log_range(b.freq_min, b.freq_max);
        Gabor2D {
            cy: pos[0],
            cx: pos[1],
            freq,
            theta: rng.range(0.0, std::f32::consts::PI),
            phase: rng.range(0.0, std::f32::consts::TAU),
            sigma: sigma_from_freq(freq, b, rng),
        }
    }

    fn mutate(&self, rng: &mut Rng, s: f32, b: &AtomBounds, dims: &[usize]) -> Self {
        let (h, w) = (dims[0] as f32, dims[1] as f32);
        Gabor2D {
            cx: (self.cx + rng.normal() * s * self.sigma).clamp(0.0, w - 1.0),
            cy: (self.cy + rng.normal() * s * self.sigma).clamp(0.0, h - 1.0),
            freq: (self.freq * (rng.normal() * s * 0.5).exp()).clamp(b.freq_min, b.freq_max),
            theta: self.theta + rng.normal() * s * 0.6,
            phase: self.phase + rng.normal() * s * 1.2,
            sigma: (self.sigma * (rng.normal() * s * 0.5).exp()).clamp(b.sigma_min, b.sigma_max),
        }
    }

    fn render_unit(&self, dims: &[usize]) -> Patch {
        let (h, w) = (dims[0], dims[1]);
        let r = (2.5 * self.sigma).ceil().max(1.0) as i64;
        let x0 = (self.cx as i64 - r).max(0) as usize;
        let y0 = (self.cy as i64 - r).max(0) as usize;
        let x1 = ((self.cx as i64 + r + 1).min(w as i64)) as usize;
        let y1 = ((self.cy as i64 + r + 1).min(h as i64)) as usize;
        let (pw, ph) = (x1.saturating_sub(x0), y1.saturating_sub(y0));
        let mut runs = Vec::with_capacity(ph);
        let mut data = Vec::with_capacity(pw * ph);
        let inv2s2 = 1.0 / (2.0 * self.sigma * self.sigma);
        let (st, ct) = self.theta.sin_cos();
        let k = std::f32::consts::TAU * self.freq;
        for y in y0..y1 {
            runs.push(Run { field_start: y * w + x0, len: pw });
            let dy = y as f32 - self.cy;
            for x in x0..x1 {
                let dx = x as f32 - self.cx;
                let env = (-(dx * dx + dy * dy) * inv2s2).exp();
                let carrier = (k * (dx * ct + dy * st) + self.phase).cos();
                data.push(env * carrier);
            }
        }
        Patch { runs, data }
    }

    fn freq(&self) -> f32 {
        self.freq
    }
}

// ─────────────────────────── 3D — vídeo ───────────────────────────

/// Átomo espaço-temporal separável: Gabor 2D espacial × envelope+portadora
/// temporal. dims = [t, h, w]; pos = [t, y, x].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Gabor3D {
    pub cx: f32,
    pub cy: f32,
    pub ct: f32,
    /// ciclos por pixel (espacial)
    pub freq_s: f32,
    pub theta: f32,
    pub phase: f32,
    pub sigma_s: f32,
    /// ciclos por frame (temporal; 0 ≈ estático)
    pub freq_t: f32,
    pub sigma_t: f32,
}

impl Atom for Gabor3D {
    const PARAMS: usize = 10;

    fn rain_at(pos: &[f32], b: &AtomBounds, rng: &mut Rng) -> Self {
        let freq_s = rng.log_range(b.freq_min, b.freq_max);
        let freq_t = rng.log_range(b.freq_min, b.freq_max);
        Gabor3D {
            ct: pos[0],
            cy: pos[1],
            cx: pos[2],
            freq_s,
            theta: rng.range(0.0, std::f32::consts::PI),
            phase: rng.range(0.0, std::f32::consts::TAU),
            sigma_s: sigma_from_freq(freq_s, b, rng),
            freq_t,
            sigma_t: sigma_from_freq(freq_t, b, rng),
        }
    }

    fn mutate(&self, rng: &mut Rng, s: f32, b: &AtomBounds, dims: &[usize]) -> Self {
        let (t, h, w) = (dims[0] as f32, dims[1] as f32, dims[2] as f32);
        Gabor3D {
            cx: (self.cx + rng.normal() * s * self.sigma_s).clamp(0.0, w - 1.0),
            cy: (self.cy + rng.normal() * s * self.sigma_s).clamp(0.0, h - 1.0),
            ct: (self.ct + rng.normal() * s * self.sigma_t).clamp(0.0, t - 1.0),
            freq_s: (self.freq_s * (rng.normal() * s * 0.5).exp()).clamp(b.freq_min, b.freq_max),
            theta: self.theta + rng.normal() * s * 0.6,
            phase: self.phase + rng.normal() * s * 1.2,
            sigma_s: (self.sigma_s * (rng.normal() * s * 0.5).exp()).clamp(b.sigma_min, b.sigma_max),
            freq_t: (self.freq_t * (rng.normal() * s * 0.5).exp()).clamp(b.freq_min, b.freq_max),
            sigma_t: (self.sigma_t * (rng.normal() * s * 0.5).exp()).clamp(b.sigma_min, b.sigma_max),
        }
    }

    fn render_unit(&self, dims: &[usize]) -> Patch {
        let (t, h, w) = (dims[0], dims[1], dims[2]);
        let rs = (2.5 * self.sigma_s).ceil().max(1.0) as i64;
        let rt = (2.5 * self.sigma_t).ceil().max(1.0) as i64;
        let x0 = (self.cx as i64 - rs).max(0) as usize;
        let y0 = (self.cy as i64 - rs).max(0) as usize;
        let t0 = (self.ct as i64 - rt).max(0) as usize;
        let x1 = ((self.cx as i64 + rs + 1).min(w as i64)) as usize;
        let y1 = ((self.cy as i64 + rs + 1).min(h as i64)) as usize;
        let t1 = ((self.ct as i64 + rt + 1).min(t as i64)) as usize;
        let pw = x1.saturating_sub(x0);
        let mut runs = Vec::new();
        let mut data = Vec::new();
        let inv2ss = 1.0 / (2.0 * self.sigma_s * self.sigma_s);
        let inv2st = 1.0 / (2.0 * self.sigma_t * self.sigma_t);
        let (st, ct_) = self.theta.sin_cos();
        let ks = std::f32::consts::TAU * self.freq_s;
        let kt = std::f32::consts::TAU * self.freq_t;
        for ti in t0..t1 {
            let dt = ti as f32 - self.ct;
            let env_t = (-(dt * dt) * inv2st).exp();
            for y in y0..y1 {
                runs.push(Run { field_start: (ti * h + y) * w + x0, len: pw });
                let dy = y as f32 - self.cy;
                for x in x0..x1 {
                    let dx = x as f32 - self.cx;
                    let env = (-(dx * dx + dy * dy) * inv2ss).exp() * env_t;
                    let carrier = (ks * (dx * ct_ + dy * st) + kt * dt + self.phase).cos();
                    data.push(env * carrier);
                }
            }
        }
        Patch { runs, data }
    }

    fn freq(&self) -> f32 {
        self.freq_s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bounds() -> AtomBounds {
        AtomBounds { freq_min: 0.02, freq_max: 0.45, sigma_min: 1.0, sigma_max: 20.0 }
    }

    #[test]
    fn planted_atom_recovered_2d() {
        let g = Gabor2D { cx: 32.0, cy: 32.0, freq: 0.1, theta: 0.6, phase: 0.3, sigma: 6.0 };
        let unit = g.render_unit(&[64, 64]);
        let mut field = Field::zeros(&[64, 64]);
        field.add_scaled_patch(&unit, 0.7);
        let (amp, explained) = fit_amp(&unit, &field).unwrap();
        assert!((amp - 0.7).abs() < 1e-4, "amp={amp}");
        let expect = 0.49 * unit.norm2();
        assert!((explained - expect).abs() / expect < 1e-3);
    }

    #[test]
    fn planted_atom_recovered_1d() {
        let g = Gabor1D { c: 500.0, freq: 0.05, phase: 1.0, sigma: 12.0 };
        let unit = g.render_unit(&[1024]);
        let mut field = Field::zeros(&[1024]);
        field.add_scaled_patch(&unit, -0.4);
        let (amp, _) = fit_amp(&unit, &field).unwrap();
        assert!((amp + 0.4).abs() < 1e-4, "amp={amp}");
    }

    #[test]
    fn planted_atom_recovered_3d() {
        let g = Gabor3D {
            cx: 16.0, cy: 16.0, ct: 6.0,
            freq_s: 0.12, theta: 0.4, phase: 0.2,
            sigma_s: 4.0, freq_t: 0.08, sigma_t: 3.0,
        };
        let dims = [12usize, 32, 32];
        let unit = g.render_unit(&dims);
        assert!(unit.norm2() > 0.0);
        let mut field = Field::zeros(&dims);
        field.add_scaled_patch(&unit, 0.55);
        let (amp, _) = fit_amp(&unit, &field).unwrap();
        assert!((amp - 0.55).abs() < 1e-4, "amp={amp}");
    }

    #[test]
    fn mutate_stays_in_bounds_all_modalities() {
        let b = bounds();
        let mut rng = Rng::new(9);
        let mut g2 = Gabor2D::rain_at(&[10.0, 10.0], &b, &mut rng);
        for _ in 0..300 {
            g2 = g2.mutate(&mut rng, 0.4, &b, &[64, 64]);
            assert!(g2.freq >= b.freq_min && g2.freq <= b.freq_max);
            assert!(g2.sigma >= b.sigma_min && g2.sigma <= b.sigma_max);
            assert!(g2.cx >= 0.0 && g2.cx <= 63.0 && g2.cy >= 0.0 && g2.cy <= 63.0);
        }
        let mut g1 = Gabor1D::rain_at(&[100.0], &b, &mut rng);
        for _ in 0..300 {
            g1 = g1.mutate(&mut rng, 0.4, &b, &[1024]);
            assert!(g1.c >= 0.0 && g1.c <= 1023.0);
            assert!(g1.freq >= b.freq_min && g1.freq <= b.freq_max);
        }
        let mut g3 = Gabor3D::rain_at(&[5.0, 10.0, 10.0], &b, &mut rng);
        for _ in 0..300 {
            g3 = g3.mutate(&mut rng, 0.4, &b, &[12, 32, 32]);
            assert!(g3.ct >= 0.0 && g3.ct <= 11.0);
            assert!(g3.cx >= 0.0 && g3.cx <= 31.0);
        }
    }

    #[test]
    fn render_clips_at_border() {
        let g = Gabor2D { cx: 1.0, cy: 1.0, freq: 0.1, theta: 0.0, phase: 0.0, sigma: 8.0 };
        let p = g.render_unit(&[32, 32]);
        assert!(p.norm2() > 0.0);
        for run in &p.runs {
            assert!(run.field_start + run.len <= 32 * 32);
        }
    }
}
