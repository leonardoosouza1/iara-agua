//! agua-atoms — os "seres/alfinetes": funções-base locais paramétricas.
//!
//! v0: Gabor 2D (envelope gaussiano × portadora senoidal). Um Gabor com
//! freq→0 degenera num blob gaussiano, então a mesma família cobre do
//! grosso (baixa freq) ao fino (alta freq) — a pirâmide pode emergir.
//!
//! A AMPLITUDE nunca é evoluída: é ajustada em forma fechada contra o
//! resíduo (amp* = ⟨r,ψ⟩/||ψ||²) — é o ser "bebendo" o resíduo. A GEOMETRIA
//! (posição/freq/orientação/fase/escala) é o que a sobrevivência esculpe.

use agua_core::{Field2D, Patch, Rng};
use serde::{Deserialize, Serialize};

/// Limites do espaço de parâmetros (por experimento, vem do TOML).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AtomBounds {
    pub freq_min: f32,
    pub freq_max: f32,
    pub sigma_min: f32,
    pub sigma_max: f32,
}

/// Átomo Gabor 2D. `amp` vive fora (no ser) — aqui só a geometria.
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

impl Gabor2D {
    /// Parâmetros por átomo (6 geometria + 1 amplitude) — pro custo em bytes.
    pub const PARAMS: usize = 7;

    /// Ser novo "chovendo" numa posição: geometria aleatória dentro dos limites.
    /// sigma acoplado à freq no nascimento (~meio ciclo dentro do envelope),
    /// depois a mutação desacopla.
    pub fn rain_at(x: f32, y: f32, b: &AtomBounds, rng: &mut Rng) -> Self {
        let freq = rng.log_range(b.freq_min, b.freq_max);
        let cycles = rng.range(0.35, 0.9);
        let sigma = (cycles / freq).clamp(b.sigma_min, b.sigma_max);
        Gabor2D {
            cx: x,
            cy: y,
            freq,
            theta: rng.range(0.0, std::f32::consts::PI),
            phase: rng.range(0.0, std::f32::consts::TAU),
            sigma,
        }
    }

    /// Filho/candidato: jitter multiplicativo em escalas, aditivo em pos/ângulos.
    pub fn mutate(&self, rng: &mut Rng, s: f32, b: &AtomBounds, w: f32, h: f32) -> Self {
        let mut g = self.clone();
        g.cx = (g.cx + rng.normal() * s * g.sigma).clamp(0.0, w - 1.0);
        g.cy = (g.cy + rng.normal() * s * g.sigma).clamp(0.0, h - 1.0);
        g.freq = (g.freq * (rng.normal() * s * 0.5).exp()).clamp(b.freq_min, b.freq_max);
        g.theta += rng.normal() * s * 0.6;
        g.phase += rng.normal() * s * 1.2;
        g.sigma = (g.sigma * (rng.normal() * s * 0.5).exp()).clamp(b.sigma_min, b.sigma_max);
        g
    }

    /// Renderiza o patch UNITÁRIO (amp=1) recortado ao campo w×h.
    /// Suporte = 2.5σ (fora disso o envelope é ~0).
    pub fn render_unit(&self, w: usize, h: usize) -> Patch {
        let r = (2.5 * self.sigma).ceil().max(1.0) as i64;
        let x0 = (self.cx as i64 - r).max(0) as usize;
        let y0 = (self.cy as i64 - r).max(0) as usize;
        let x1 = ((self.cx as i64 + r + 1).min(w as i64)) as usize;
        let y1 = ((self.cy as i64 + r + 1).min(h as i64)) as usize;
        let (pw, ph) = (x1.saturating_sub(x0), y1.saturating_sub(y0));
        let mut data = vec![0.0f32; pw * ph];
        let inv2s2 = 1.0 / (2.0 * self.sigma * self.sigma);
        let (st, ct) = self.theta.sin_cos();
        let k = std::f32::consts::TAU * self.freq;
        for py in 0..ph {
            let dy = (y0 + py) as f32 - self.cy;
            for px in 0..pw {
                let dx = (x0 + px) as f32 - self.cx;
                let env = (-(dx * dx + dy * dy) * inv2s2).exp();
                let carrier = (k * (dx * ct + dy * st) + self.phase).cos();
                data[py * pw + px] = env * carrier;
            }
        }
        Patch { x0, y0, w: pw, h: ph, data }
    }
}

/// Ajuste de amplitude em forma fechada: amp* = ⟨residual,ψ⟩/||ψ||².
/// Retorna (amp, energia_explicada = amp²·||ψ||²). `None` se suporte degenerado.
pub fn fit_amp(unit: &Patch, residual: &Field2D) -> Option<(f32, f32)> {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn bounds() -> AtomBounds {
        AtomBounds { freq_min: 0.02, freq_max: 0.45, sigma_min: 1.0, sigma_max: 20.0 }
    }

    #[test]
    fn planted_atom_recovered() {
        // planta um gabor com amp 0.7; fit com a MESMA geometria recupera amp
        let g = Gabor2D { cx: 32.0, cy: 32.0, freq: 0.1, theta: 0.6, phase: 0.3, sigma: 6.0 };
        let unit = g.render_unit(64, 64);
        let mut field = Field2D::zeros(64, 64);
        field.add_scaled_patch(&unit, 0.7);
        let (amp, explained) = fit_amp(&unit, &field).unwrap();
        assert!((amp - 0.7).abs() < 1e-4, "amp={amp}");
        let expect = 0.49 * unit.norm2();
        assert!((explained - expect).abs() / expect < 1e-3);
    }

    #[test]
    fn mutate_stays_in_bounds() {
        let b = bounds();
        let mut rng = Rng::new(9);
        let mut g = Gabor2D::rain_at(10.0, 10.0, &b, &mut rng);
        for _ in 0..500 {
            g = g.mutate(&mut rng, 0.4, &b, 64.0, 64.0);
            assert!(g.freq >= b.freq_min && g.freq <= b.freq_max);
            assert!(g.sigma >= b.sigma_min && g.sigma <= b.sigma_max);
            assert!(g.cx >= 0.0 && g.cx <= 63.0 && g.cy >= 0.0 && g.cy <= 63.0);
        }
    }

    #[test]
    fn render_clips_at_border() {
        let g = Gabor2D { cx: 1.0, cy: 1.0, freq: 0.1, theta: 0.0, phase: 0.0, sigma: 8.0 };
        let p = g.render_unit(32, 32);
        assert_eq!(p.x0, 0);
        assert_eq!(p.y0, 0);
        assert!(p.w <= 32 && p.h <= 32);
        assert!(p.norm2() > 0.0);
    }
}
