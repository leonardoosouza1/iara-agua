//! agua-core — primitivas fundamentais do solver-por-sobrevivência.
//!
//! Zero dependências externas. Tudo determinístico: RNG próprio seedado,
//! sem `Instant`/`random` escondidos. Este crate NÃO sabe o que é um átomo
//! nem o que é sobrevivência — só campos 2D, patches esparsos e amostragem.

/// RNG xorshift64* — determinístico, barato, seedável.
#[derive(Clone, Debug)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut s = seed
            .wrapping_mul(0x9E37_79B9_7F4A_7C15)
            .wrapping_add(0x2545_F491_4F6C_DD1D);
        if s == 0 {
            s = 0xDEAD_BEEF_CAFE_F00D;
        }
        Rng(s)
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniforme em [0, 1).
    #[inline]
    pub fn f32(&mut self) -> f32 {
        (self.next_u64() >> 40) as f32 / (1u64 << 24) as f32
    }

    #[inline]
    pub fn range(&mut self, lo: f32, hi: f32) -> f32 {
        lo + (hi - lo) * self.f32()
    }

    /// Uniforme em escala log — certo pra frequências/escalas.
    #[inline]
    pub fn log_range(&mut self, lo: f32, hi: f32) -> f32 {
        self.range(lo.ln(), hi.ln()).exp()
    }

    #[inline]
    pub fn usize(&mut self, n: usize) -> usize {
        if n == 0 {
            0
        } else {
            (self.next_u64() % n as u64) as usize
        }
    }

    /// Normal(0,1) via Box-Muller.
    pub fn normal(&mut self) -> f32 {
        let u1 = self.f32().max(1e-7);
        let u2 = self.f32();
        (-2.0 * u1.ln()).sqrt() * (std::f32::consts::TAU * u2).cos()
    }

    /// Embaralha índices 0..n (Fisher-Yates) — ordem de atualização dos seres.
    pub fn shuffled_indices(&mut self, n: usize) -> Vec<usize> {
        let mut idx: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = self.usize(i + 1);
            idx.swap(i, j);
        }
        idx
    }
}

/// Campo escalar 2D denso (o sinal, a reconstrução, o resíduo).
#[derive(Clone, Debug)]
pub struct Field2D {
    pub w: usize,
    pub h: usize,
    pub data: Vec<f32>,
}

impl Field2D {
    pub fn zeros(w: usize, h: usize) -> Self {
        Field2D { w, h, data: vec![0.0; w * h] }
    }

    pub fn from_fn(w: usize, h: usize, f: impl Fn(usize, usize) -> f32) -> Self {
        let mut data = Vec::with_capacity(w * h);
        for y in 0..h {
            for x in 0..w {
                data.push(f(x, y));
            }
        }
        Field2D { w, h, data }
    }

    #[inline]
    pub fn at(&self, x: usize, y: usize) -> f32 {
        self.data[y * self.w + x]
    }

    #[inline]
    pub fn at_mut(&mut self, x: usize, y: usize) -> &mut f32 {
        &mut self.data[y * self.w + x]
    }

    /// Energia total (soma dos quadrados).
    pub fn energy(&self) -> f32 {
        self.data.iter().map(|v| v * v).sum()
    }

    pub fn mean(&self) -> f32 {
        self.data.iter().sum::<f32>() / self.data.len() as f32
    }

    pub fn mse_vs(&self, other: &Field2D) -> f32 {
        assert_eq!(self.data.len(), other.data.len());
        self.data
            .iter()
            .zip(&other.data)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f32>()
            / self.data.len() as f32
    }

    /// Soma `k * patch` no campo (k = +amp pra aplicar, -amp pra remover).
    pub fn add_scaled_patch(&mut self, p: &Patch, k: f32) {
        for py in 0..p.h {
            let row = (p.y0 + py) * self.w + p.x0;
            let prow = py * p.w;
            for px in 0..p.w {
                self.data[row + px] += k * p.data[prow + px];
            }
        }
    }
}

/// Patch esparso: contribuição local de um átomo (janela retangular).
#[derive(Clone, Debug)]
pub struct Patch {
    pub x0: usize,
    pub y0: usize,
    pub w: usize,
    pub h: usize,
    pub data: Vec<f32>,
}

impl Patch {
    /// ⟨patch, campo restrito ao suporte⟩.
    pub fn dot_field(&self, f: &Field2D) -> f32 {
        let mut acc = 0.0f32;
        for py in 0..self.h {
            let row = (self.y0 + py) * f.w + self.x0;
            let prow = py * self.w;
            for px in 0..self.w {
                acc += self.data[prow + px] * f.data[row + px];
            }
        }
        acc
    }

    /// ⟨self, other⟩ sobre a interseção dos suportes (0 se não sobrepõem).
    pub fn dot_patch(&self, other: &Patch) -> f32 {
        let x0 = self.x0.max(other.x0);
        let y0 = self.y0.max(other.y0);
        let x1 = (self.x0 + self.w).min(other.x0 + other.w);
        let y1 = (self.y0 + self.h).min(other.y0 + other.h);
        if x1 <= x0 || y1 <= y0 {
            return 0.0;
        }
        let mut acc = 0.0f32;
        for y in y0..y1 {
            let srow = (y - self.y0) * self.w + (x0 - self.x0);
            let orow = (y - other.y0) * other.w + (x0 - other.x0);
            for i in 0..(x1 - x0) {
                acc += self.data[srow + i] * other.data[orow + i];
            }
        }
        acc
    }

    /// ||patch||².
    pub fn norm2(&self) -> f32 {
        self.data.iter().map(|v| v * v).sum()
    }
}

/// Amostra uma posição ∝ energia local do campo (blocos grid×grid).
/// É a "comida sobrando": onde o resíduo é alto, chove ser novo.
pub fn sample_hotspot(f: &Field2D, rng: &mut Rng, grid: usize) -> (usize, usize) {
    let grid = grid.max(1).min(f.w.min(f.h));
    let bw = (f.w + grid - 1) / grid;
    let bh = (f.h + grid - 1) / grid;
    let mut weights = vec![0.0f32; grid * grid];
    for y in 0..f.h {
        for x in 0..f.w {
            let v = f.at(x, y);
            weights[(y / bh).min(grid - 1) * grid + (x / bw).min(grid - 1)] += v * v;
        }
    }
    let total: f32 = weights.iter().sum();
    let block = if total <= 1e-12 {
        rng.usize(grid * grid)
    } else {
        let mut r = rng.f32() * total;
        let mut chosen = grid * grid - 1;
        for (i, w) in weights.iter().enumerate() {
            if r < *w {
                chosen = i;
                break;
            }
            r -= w;
        }
        chosen
    };
    let bx = block % grid;
    let by = block / grid;
    let x = (bx * bw + rng.usize(bw)).min(f.w - 1);
    let y = (by * bh + rng.usize(bh)).min(f.h - 1);
    (x, y)
}

/// PSNR em dB dado MSE e pico do sinal.
pub fn psnr(mse: f32, peak: f32) -> f32 {
    if mse <= 1e-12 {
        return 99.0;
    }
    10.0 * (peak * peak / mse).log10()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rng_deterministic() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
        let mut c = Rng::new(43);
        assert_ne!(a.next_u64(), c.next_u64());
    }

    #[test]
    fn rng_f32_in_range() {
        let mut r = Rng::new(7);
        for _ in 0..1000 {
            let v = r.f32();
            assert!((0.0..1.0).contains(&v));
        }
    }

    #[test]
    fn patch_roundtrip_zeroes() {
        let mut f = Field2D::zeros(16, 16);
        let p = Patch {
            x0: 3,
            y0: 5,
            w: 4,
            h: 4,
            data: (0..16).map(|i| i as f32 * 0.1).collect(),
        };
        f.add_scaled_patch(&p, 2.5);
        f.add_scaled_patch(&p, -2.5);
        assert!(f.data.iter().all(|v| v.abs() < 1e-6));
    }

    #[test]
    fn dot_patch_overlap() {
        let a = Patch { x0: 0, y0: 0, w: 2, h: 2, data: vec![1.0; 4] };
        let b = Patch { x0: 1, y0: 1, w: 2, h: 2, data: vec![3.0; 4] };
        // interseção = 1 pixel → 1*3
        assert!((a.dot_patch(&b) - 3.0).abs() < 1e-6);
        let c = Patch { x0: 10, y0: 10, w: 2, h: 2, data: vec![1.0; 4] };
        assert_eq!(a.dot_patch(&c), 0.0);
    }

    #[test]
    fn hotspot_prefers_energy() {
        // energia concentrada no canto inferior direito
        let f = Field2D::from_fn(32, 32, |x, y| if x >= 24 && y >= 24 { 10.0 } else { 0.0 });
        let mut rng = Rng::new(1);
        let mut hits = 0;
        for _ in 0..200 {
            let (x, y) = sample_hotspot(&f, &mut rng, 8);
            if x >= 24 && y >= 24 {
                hits += 1;
            }
        }
        assert!(hits > 180, "hotspot ignorou a energia: {hits}/200");
    }
}
