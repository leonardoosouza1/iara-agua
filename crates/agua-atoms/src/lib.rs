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

    /// nº de parâmetros de GEOMETRIA (os que o gradiente move; amp é forma fechada)
    const GEOM: usize;

    /// ∂E/∂p para cada parâmetro de geometria, onde E = ½‖r_sem_eu − amp·ψ‖².
    /// Como r_atual = r_sem_eu − amp·ψ, vale ∂E/∂p = −amp·⟨r_atual, ∂ψ/∂p⟩ —
    /// tudo local ao suporte do átomo (custo de um render).
    fn grad(&self, dims: &[usize], residual: &Field, amp: f32) -> Vec<f32>;

    /// Passo no espaço de parâmetros, com clamp nos limites.
    fn apply_delta(&mut self, delta: &[f32], b: &AtomBounds, dims: &[usize]);

    /// Escala típica de cada parâmetro (pro lr do Adam ser adimensional).
    fn param_scale(&self, dims: &[usize]) -> Vec<f32>;

    /// Centro no espaço do campo (ordem de `dims`) — pra vizinhança espacial.
    fn center(&self) -> Vec<f32>;
    /// Raio do suporte — dois seres só podem se predar se os suportes se tocam.
    fn radius(&self) -> f32;
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

    const GEOM: usize = 4; // c, freq, phase, sigma

    fn grad(&self, dims: &[usize], residual: &Field, amp: f32) -> Vec<f32> {
        let n = dims[0];
        let r = (2.5 * self.sigma).ceil().max(1.0) as i64;
        let x0 = (self.c as i64 - r).max(0) as usize;
        let x1 = ((self.c as i64 + r + 1).min(n as i64)) as usize;
        let inv2s2 = 1.0 / (2.0 * self.sigma * self.sigma);
        let invs2 = 1.0 / (self.sigma * self.sigma);
        let invs3 = invs2 / self.sigma;
        let k = std::f32::consts::TAU * self.freq;
        let mut g = vec![0.0f32; 4];
        for x in x0..x1 {
            let dx = x as f32 - self.c;
            let env = (-(dx * dx) * inv2s2).exp();
            let arg = k * dx + self.phase;
            let (sa, ca) = arg.sin_cos();
            let res = residual.data[x];
            let dpsi_c = env * (dx * invs2) * ca + env * k * sa;
            let dpsi_freq = -env * sa * std::f32::consts::TAU * dx;
            let dpsi_phase = -env * sa;
            let dpsi_sigma = env * (dx * dx * invs3) * ca;
            g[0] -= amp * res * dpsi_c;
            g[1] -= amp * res * dpsi_freq;
            g[2] -= amp * res * dpsi_phase;
            g[3] -= amp * res * dpsi_sigma;
        }
        g
    }

    fn apply_delta(&mut self, d: &[f32], b: &AtomBounds, dims: &[usize]) {
        self.c = (self.c + d[0]).clamp(0.0, dims[0] as f32 - 1.0);
        self.freq = (self.freq + d[1]).clamp(b.freq_min, b.freq_max);
        self.phase += d[2];
        self.sigma = (self.sigma + d[3]).clamp(b.sigma_min, b.sigma_max);
    }

    fn param_scale(&self, _dims: &[usize]) -> Vec<f32> {
        let s = self.sigma.max(1.0);
        vec![s, self.freq.max(1e-3), 1.0, s]
    }

    fn center(&self) -> Vec<f32> {
        vec![self.c]
    }
    fn radius(&self) -> f32 {
        2.5 * self.sigma
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

    const GEOM: usize = 6; // cx, cy, freq, theta, phase, sigma

    /// ψ = env·car, env = exp(−(dx²+dy²)/2σ²), car = cos(k·u + φ), u = dx·cosθ + dy·sinθ
    fn grad(&self, dims: &[usize], residual: &Field, amp: f32) -> Vec<f32> {
        let (h, w) = (dims[0], dims[1]);
        let r = (2.5 * self.sigma).ceil().max(1.0) as i64;
        let x0 = (self.cx as i64 - r).max(0) as usize;
        let y0 = (self.cy as i64 - r).max(0) as usize;
        let x1 = ((self.cx as i64 + r + 1).min(w as i64)) as usize;
        let y1 = ((self.cy as i64 + r + 1).min(h as i64)) as usize;
        let inv2s2 = 1.0 / (2.0 * self.sigma * self.sigma);
        let invs2 = 1.0 / (self.sigma * self.sigma);
        let invs3 = invs2 / self.sigma;
        let (st, ct) = self.theta.sin_cos();
        let k = std::f32::consts::TAU * self.freq;
        let mut g = vec![0.0f32; 6];
        for y in y0..y1 {
            let dy = y as f32 - self.cy;
            let row = y * w;
            for x in x0..x1 {
                let dx = x as f32 - self.cx;
                let rr = dx * dx + dy * dy;
                let env = (-rr * inv2s2).exp();
                let u = dx * ct + dy * st;
                let arg = k * u + self.phase;
                let (sa, ca) = arg.sin_cos();
                let res = residual.data[row + x];
                // ∂ψ/∂p
                let dpsi_cx = env * (dx * invs2) * ca + env * k * ct * sa;
                let dpsi_cy = env * (dy * invs2) * ca + env * k * st * sa;
                let dpsi_freq = -env * sa * std::f32::consts::TAU * u;
                let dpsi_theta = -env * sa * k * (-dx * st + dy * ct);
                let dpsi_phase = -env * sa;
                let dpsi_sigma = env * (rr * invs3) * ca;
                // ∂E/∂p = −amp·⟨r, ∂ψ/∂p⟩
                g[0] -= amp * res * dpsi_cx;
                g[1] -= amp * res * dpsi_cy;
                g[2] -= amp * res * dpsi_freq;
                g[3] -= amp * res * dpsi_theta;
                g[4] -= amp * res * dpsi_phase;
                g[5] -= amp * res * dpsi_sigma;
            }
        }
        g
    }

    fn apply_delta(&mut self, d: &[f32], b: &AtomBounds, dims: &[usize]) {
        let (h, w) = (dims[0] as f32, dims[1] as f32);
        self.cx = (self.cx + d[0]).clamp(0.0, w - 1.0);
        self.cy = (self.cy + d[1]).clamp(0.0, h - 1.0);
        self.freq = (self.freq + d[2]).clamp(b.freq_min, b.freq_max);
        self.theta += d[3];
        self.phase += d[4];
        self.sigma = (self.sigma + d[5]).clamp(b.sigma_min, b.sigma_max);
    }

    fn param_scale(&self, _dims: &[usize]) -> Vec<f32> {
        let s = self.sigma.max(1.0);
        // cada parâmetro anda em unidades da SUA escala natural
        vec![s, s, self.freq.max(1e-3), 1.0, 1.0, s]
    }

    fn center(&self) -> Vec<f32> {
        vec![self.cy, self.cx]
    }
    fn radius(&self) -> f32 {
        2.5 * self.sigma
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

    const GEOM: usize = 0; // gradiente 3D ainda não derivado — usa mutação

    fn grad(&self, _dims: &[usize], _residual: &Field, _amp: f32) -> Vec<f32> {
        Vec::new()
    }

    fn apply_delta(&mut self, _d: &[f32], _b: &AtomBounds, _dims: &[usize]) {}

    fn param_scale(&self, _dims: &[usize]) -> Vec<f32> {
        Vec::new()
    }

    fn center(&self) -> Vec<f32> {
        vec![self.ct, self.cy, self.cx]
    }
    fn radius(&self) -> f32 {
        2.5 * self.sigma_s.max(self.sigma_t)
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

#[cfg(test)]
mod grad_tests {
    use super::*;

    /// Energia RELATIVA em f64: E(p) = const − amp·⟨r_w,ψ⟩ + ½·amp²·‖ψ‖².
    /// O termo ½‖r_w‖² é constante em p e cancela na diferença central — omiti-lo
    /// evita cancelamento catastrófico em f32 (E≈100, ΔE≈1e-4).
    fn rel_energy(unit: &Patch, r_without: &Field, amp: f32) -> f64 {
        let dot = unit.dot_field(r_without) as f64;
        let n2 = unit.norm2() as f64;
        -(amp as f64) * dot + 0.5 * (amp as f64) * (amp as f64) * n2
    }

    fn check<A: Atom>(g: &A, dims: &[usize], b: &AtomBounds, amp: f32, eps: &[f32], seed: u64) {
        let mut rng = Rng::new(seed);
        let n: usize = dims.iter().product();
        let r_without =
            Field { dims: dims.to_vec(), data: (0..n).map(|_| rng.normal() * 0.3).collect() };
        // resíduo ATUAL = r_sem_eu − amp·ψ (é o que grad() recebe)
        let mut residual = r_without.clone();
        residual.add_scaled_patch(&g.render_unit(dims), -amp);

        let ana = g.grad(dims, &residual, amp);
        assert_eq!(ana.len(), A::GEOM);

        for k in 0..A::GEOM {
            let mut dp = vec![0.0; A::GEOM];
            dp[k] = eps[k];
            let mut gp = g.clone();
            gp.apply_delta(&dp, b, dims);
            dp[k] = -eps[k];
            let mut gm = g.clone();
            gm.apply_delta(&dp, b, dims);
            let ep = rel_energy(&gp.render_unit(dims), &r_without, amp);
            let em = rel_energy(&gm.render_unit(dims), &r_without, amp);
            let num = ((ep - em) / (2.0 * eps[k] as f64)) as f32;
            let denom = num.abs().max(ana[k].abs()).max(1e-2);
            let rel = (num - ana[k]).abs() / denom;
            assert!(rel < 0.02, "param {k}: analítico {} vs numérico {num} (rel {rel:.4})", ana[k]);
        }
    }

    #[test]
    fn gabor2d_grad_matches_finite_differences() {
        let dims = [48usize, 48];
        let b = AtomBounds { freq_min: 0.01, freq_max: 0.45, sigma_min: 1.0, sigma_max: 20.0 };
        // centros longe de inteiros: o suporte é recortado com `as i64` e saltaria
        // 2.5·σ = 12.7 (longe de inteiro) ⇒ o raio do suporte não salta com ±eps
        let g = Gabor2D { cx: 24.3, cy: 22.7, freq: 0.11, theta: 0.8, phase: 1.3, sigma: 5.08 };
        // eps grande o bastante pra ΔE sair do ruído de f32 do dot, pequeno pra 2ª ordem não pesar
        check(&g, &dims, &b, 0.6, &[2e-3, 2e-3, 2e-5, 1e-3, 1e-3, 2e-3], 5);
    }

    #[test]
    fn gabor1d_grad_matches_finite_differences() {
        let dims = [512usize];
        let b = AtomBounds { freq_min: 0.002, freq_max: 0.45, sigma_min: 2.0, sigma_max: 100.0 };
        // 2.5·σ = 49.75 (longe de inteiro); c longe de inteiro pelo mesmo motivo
        let g = Gabor1D { c: 250.4, freq: 0.05, phase: 0.9, sigma: 19.9 };
        check(&g, &dims, &b, -0.4, &[2e-3, 2e-6, 1e-3, 2e-3], 9);
    }

    /// O gradiente aponta ladeira abaixo: um passo pequeno reduz a energia.
    #[test]
    fn gradient_descends() {
        let dims = [48usize, 48];
        let b = AtomBounds { freq_min: 0.01, freq_max: 0.45, sigma_min: 1.0, sigma_max: 20.0 };
        let g = Gabor2D { cx: 24.3, cy: 22.7, freq: 0.11, theta: 0.8, phase: 1.3, sigma: 5.0 };
        let amp = 0.6f32;
        let mut rng = Rng::new(3);
        let r_without =
            Field { dims: dims.to_vec(), data: (0..48 * 48).map(|_| rng.normal() * 0.3).collect() };
        let mut residual = r_without.clone();
        residual.add_scaled_patch(&g.render_unit(&dims), -amp);
        let e0 = rel_energy(&g.render_unit(&dims), &r_without, amp);
        let grad = g.grad(&dims, &residual, amp);
        let scale = g.param_scale(&dims);
        let gn = grad.iter().map(|v| v * v).sum::<f32>().sqrt().max(1e-9);
        let step: Vec<f32> = grad
            .iter()
            .zip(&scale)
            .map(|(gr, sc)| -1e-3 * sc * gr / gn)
            .collect();
        let mut g2 = g.clone();
        g2.apply_delta(&step, &b, &dims);
        let e1 = rel_energy(&g2.render_unit(&dims), &r_without, amp);
        assert!(e1 < e0, "energia subiu: {e0} → {e1}");
    }
}
