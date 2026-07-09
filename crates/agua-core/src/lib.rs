//! agua-core — primitivas fundamentais do solver-por-sobrevivência.
//!
//! Zero dependências externas. Tudo determinístico: RNG próprio seedado,
//! sem `Instant`/`random` escondidos. Este crate NÃO sabe o que é um átomo,
//! sobrevivência ou modalidade — só campos N-dimensionais flat, patches
//! esparsos por runs e amostragem ponderada.
//!
//! AGNOSTICISMO: `Field` guarda `dims` (row-major, dims[0] mais lento) e
//! dados flat. `Patch` é uma lista de runs contíguos no índice flat — uma
//! linha de imagem, um intervalo de áudio, uma linha de um frame de vídeo.
//! A física de vida opera só nisso; a modalidade vive nos átomos.

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

/// Campo escalar N-dimensional denso, row-major (o sinal, a reconstrução,
/// o resíduo). dims[0] é o eixo mais lento (2D: [h,w]; 3D: [t,h,w]).
#[derive(Clone, Debug)]
pub struct Field {
    pub dims: Vec<usize>,
    pub data: Vec<f32>,
}

impl Field {
    pub fn zeros(dims: &[usize]) -> Self {
        let len = dims.iter().product();
        Field { dims: dims.to_vec(), data: vec![0.0; len] }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Conveniência 1D (áudio): f(i).
    pub fn from_fn_1d(n: usize, f: impl Fn(usize) -> f32) -> Self {
        Field { dims: vec![n], data: (0..n).map(f).collect() }
    }

    /// Conveniência 2D (imagem): f(x, y), dims [h, w].
    pub fn from_fn_2d(h: usize, w: usize, f: impl Fn(usize, usize) -> f32) -> Self {
        let mut data = Vec::with_capacity(w * h);
        for y in 0..h {
            for x in 0..w {
                data.push(f(x, y));
            }
        }
        Field { dims: vec![h, w], data }
    }

    /// Conveniência 3D (vídeo): f(x, y, t), dims [t, h, w].
    pub fn from_fn_3d(t: usize, h: usize, w: usize, f: impl Fn(usize, usize, usize) -> f32) -> Self {
        let mut data = Vec::with_capacity(t * h * w);
        for ti in 0..t {
            for y in 0..h {
                for x in 0..w {
                    data.push(f(x, y, ti));
                }
            }
        }
        Field { dims: vec![t, h, w], data }
    }

    /// Energia total (soma dos quadrados).
    pub fn energy(&self) -> f32 {
        self.data.iter().map(|v| v * v).sum()
    }

    pub fn mean(&self) -> f32 {
        self.data.iter().sum::<f32>() / self.data.len() as f32
    }

    /// Soma `k * patch` no campo (k = +amp pra aplicar, -amp pra remover).
    pub fn add_scaled_patch(&mut self, p: &Patch, k: f32) {
        let mut off = 0usize;
        for run in &p.runs {
            let dst = &mut self.data[run.field_start..run.field_start + run.len];
            let src = &p.data[off..off + run.len];
            for (d, s) in dst.iter_mut().zip(src) {
                *d += k * s;
            }
            off += run.len;
        }
    }
}

/// Converte índice flat → coordenadas por eixo (na ordem de dims).
pub fn unflatten(dims: &[usize], mut i: usize) -> Vec<usize> {
    let mut out = vec![0usize; dims.len()];
    for d in (0..dims.len()).rev() {
        out[d] = i % dims[d];
        i /= dims[d];
    }
    out
}

/// Um trecho contíguo do índice flat do campo.
#[derive(Clone, Debug)]
pub struct Run {
    pub field_start: usize,
    pub len: usize,
}

/// Patch esparso: contribuição local de um átomo como runs contíguos
/// (ordenados por field_start, sem sobreposição interna). `data` concatena
/// os valores na ordem dos runs.
#[derive(Clone, Debug)]
pub struct Patch {
    pub runs: Vec<Run>,
    pub data: Vec<f32>,
}

impl Patch {
    /// ⟨patch, campo restrito ao suporte⟩.
    pub fn dot_field(&self, f: &Field) -> f32 {
        let mut acc = 0.0f32;
        let mut off = 0usize;
        for run in &self.runs {
            let src = &f.data[run.field_start..run.field_start + run.len];
            let pd = &self.data[off..off + run.len];
            for (a, b) in pd.iter().zip(src) {
                acc += a * b;
            }
            off += run.len;
        }
        acc
    }

    /// ⟨self, other⟩ sobre a interseção dos suportes (0 se não sobrepõem).
    /// Two-pointer sobre runs ordenados por field_start.
    pub fn dot_patch(&self, other: &Patch) -> f32 {
        let mut acc = 0.0f32;
        let (mut i, mut j) = (0usize, 0usize);
        let (mut aoff, mut boff) = (0usize, 0usize);
        while i < self.runs.len() && j < other.runs.len() {
            let ra = &self.runs[i];
            let rb = &other.runs[j];
            let ae = ra.field_start + ra.len;
            let be = rb.field_start + rb.len;
            let s = ra.field_start.max(rb.field_start);
            let e = ae.min(be);
            if s < e {
                let ai = aoff + (s - ra.field_start);
                let bi = boff + (s - rb.field_start);
                for k in 0..(e - s) {
                    acc += self.data[ai + k] * other.data[bi + k];
                }
            }
            if ae <= be {
                aoff += ra.len;
                i += 1;
            }
            if be <= ae {
                boff += rb.len;
                j += 1;
            }
        }
        acc
    }

    /// ||patch||².
    pub fn norm2(&self) -> f32 {
        self.data.iter().map(|v| v * v).sum()
    }
}

/// Amostrador ponderado ∝ valor² do campo — a "comida sobrando": onde o
/// resíduo é alto, chove ser novo. Exato (prefix-sum + busca binária),
/// construir 1× por chuva.
pub struct HotspotSampler {
    prefix: Vec<f32>,
    total: f32,
}

impl HotspotSampler {
    pub fn new(f: &Field) -> Self {
        let mut prefix = Vec::with_capacity(f.len());
        let mut acc = 0.0f32;
        for v in &f.data {
            acc += v * v;
            prefix.push(acc);
        }
        HotspotSampler { prefix, total: acc }
    }

    /// Índice flat amostrado ∝ v² (uniforme se campo ~zerado).
    pub fn sample(&self, rng: &mut Rng) -> usize {
        let n = self.prefix.len();
        if self.total <= 1e-12 || n == 0 {
            return rng.usize(n.max(1));
        }
        let r = rng.f32() * self.total;
        self.prefix.partition_point(|p| *p <= r).min(n - 1)
    }
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

    fn rect_patch(x0: usize, y0: usize, w: usize, h: usize, fw: usize, val: f32) -> Patch {
        let mut runs = Vec::new();
        let mut data = Vec::new();
        for y in y0..y0 + h {
            runs.push(Run { field_start: y * fw + x0, len: w });
            for _ in 0..w {
                data.push(val);
            }
        }
        Patch { runs, data }
    }

    #[test]
    fn patch_roundtrip_zeroes() {
        let mut f = Field::zeros(&[16, 16]);
        let p = rect_patch(3, 5, 4, 4, 16, 0.7);
        f.add_scaled_patch(&p, 2.5);
        f.add_scaled_patch(&p, -2.5);
        assert!(f.data.iter().all(|v| v.abs() < 1e-6));
    }

    #[test]
    fn dot_patch_overlap_2d() {
        // interseção 1 pixel (canto) → 1*3
        let a = rect_patch(0, 0, 2, 2, 16, 1.0);
        let b = rect_patch(1, 1, 2, 2, 16, 3.0);
        assert!((a.dot_patch(&b) - 3.0).abs() < 1e-6);
        let c = rect_patch(10, 10, 2, 2, 16, 1.0);
        assert_eq!(a.dot_patch(&c), 0.0);
        // consistência com dot_field
        let mut f = Field::zeros(&[16, 16]);
        f.add_scaled_patch(&b, 1.0);
        assert!((a.dot_field(&f) - a.dot_patch(&b)).abs() < 1e-6);
    }

    #[test]
    fn dot_patch_1d_interval() {
        let a = Patch { runs: vec![Run { field_start: 5, len: 4 }], data: vec![2.0; 4] };
        let b = Patch { runs: vec![Run { field_start: 7, len: 4 }], data: vec![1.5; 4] };
        // interseção = índices 7,8 → 2*(2.0*1.5)
        assert!((a.dot_patch(&b) - 6.0).abs() < 1e-6);
    }

    #[test]
    fn unflatten_roundtrip() {
        let dims = [3usize, 4, 5];
        assert_eq!(unflatten(&dims, 0), vec![0, 0, 0]);
        assert_eq!(unflatten(&dims, 59), vec![2, 3, 4]);
        assert_eq!(unflatten(&dims, 27), vec![1, 1, 2]); // 1*20 + 1*5 + 2
    }

    #[test]
    fn hotspot_prefers_energy() {
        // energia concentrada no canto inferior direito
        let f = Field::from_fn_2d(32, 32, |x, y| if x >= 24 && y >= 24 { 10.0 } else { 0.0 });
        let mut rng = Rng::new(1);
        let sampler = HotspotSampler::new(&f);
        let mut hits = 0;
        for _ in 0..200 {
            let c = unflatten(&f.dims, sampler.sample(&mut rng));
            if c[1] >= 24 && c[0] >= 24 {
                hits += 1;
            }
        }
        assert_eq!(hits, 200, "sampler exato deve cair só onde há energia");
    }
}
