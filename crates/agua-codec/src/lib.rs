//! agua-codec — a população vira BYTES (e volta).
//!
//! Os parâmetros dos seres vivos SÃO a codificação do sinal. Em f32 cru cada
//! Gabor2D custa 28 B — generoso demais pra comparar com JPEG. Aqui eles
//! viram um bitstream compacto: cada campo recebe só os bits que a sua
//! sensibilidade justifica (posição precisa; fase/orientação não).
//!
//! Perfis:
//!   `F32`     — 32 bits/campo (referência, sem perda)
//!   `Fp16`    — IEEE half em tudo (2× menor, perda desprezível)
//!   `Q12`     — pos 12b, escalas log 8b, ângulos 8b, amp fp16 (~9 B/átomo)
//!   `Q8`      — pos 10b, freq/sigma 6b, ângulos 6b, amp 10b (~6 B/átomo)
//!
//! Formato: header ("AGUA" | v | modality | ndims | dims | dc | n | perfil)
//! seguido do bitstream dos átomos. Roundtrip testado bit-exato na estrutura.

use agua_atoms::{AtomBounds, Gabor1D, Gabor2D};
use agua_core::{Field, Patch};

pub const MAGIC: &[u8; 4] = b"AGUA";
pub const VERSION: u8 = 1;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Profile {
    F32,
    Fp16,
    Q12,
    Q8,
}

impl Profile {
    pub fn from_str(s: &str) -> Option<Profile> {
        match s {
            "f32" => Some(Profile::F32),
            "fp16" => Some(Profile::Fp16),
            "q12" => Some(Profile::Q12),
            "q8" => Some(Profile::Q8),
            _ => None,
        }
    }
    pub fn tag(self) -> u8 {
        match self {
            Profile::F32 => 0,
            Profile::Fp16 => 1,
            Profile::Q12 => 2,
            Profile::Q8 => 3,
        }
    }
    pub fn from_tag(t: u8) -> Option<Profile> {
        match t {
            0 => Some(Profile::F32),
            1 => Some(Profile::Fp16),
            2 => Some(Profile::Q12),
            3 => Some(Profile::Q8),
            _ => None,
        }
    }
    /// bits: (posição, freq, ângulo, sigma, amp)
    fn bits(self) -> (u32, u32, u32, u32, u32) {
        match self {
            Profile::F32 => (32, 32, 32, 32, 32),
            Profile::Fp16 => (16, 16, 16, 16, 16),
            Profile::Q12 => (12, 8, 8, 8, 16),
            Profile::Q8 => (10, 6, 6, 6, 10),
        }
    }
}

// ─────────────────────────── fp16 (IEEE half, sem deps) ───────────────────────────

pub fn f32_to_f16(v: f32) -> u16 {
    let b = v.to_bits();
    let sign = ((b >> 16) & 0x8000) as u16;
    let mut exp = ((b >> 23) & 0xff) as i32 - 127 + 15;
    let mant = b & 0x007f_ffff;
    if exp <= 0 {
        // subnormal ou zero
        if exp < -10 {
            return sign;
        }
        let mant = mant | 0x0080_0000;
        let shift = (14 - exp) as u32;
        let half = (mant >> shift) as u16;
        // arredondamento ao mais próximo
        let round = ((mant >> (shift - 1)) & 1) as u16;
        return sign | (half + round);
    }
    if exp >= 0x1f {
        return sign | 0x7c00; // inf
    }
    let mut half_mant = (mant >> 13) as u16;
    if (mant & 0x1000) != 0 {
        half_mant += 1;
        if half_mant == 0x400 {
            half_mant = 0;
            exp += 1;
            if exp >= 0x1f {
                return sign | 0x7c00;
            }
        }
    }
    sign | ((exp as u16) << 10) | half_mant
}

pub fn f16_to_f32(h: u16) -> f32 {
    let sign = ((h & 0x8000) as u32) << 16;
    let exp = ((h >> 10) & 0x1f) as u32;
    let mant = (h & 0x03ff) as u32;
    if exp == 0 {
        if mant == 0 {
            return f32::from_bits(sign);
        }
        // subnormal → normaliza. valor = mant·2^-24; após s deslocamentos
        // m = 1.xxx·2^10 ⇒ valor = 1.xxx·2^(-14-s) ⇒ expoente f32 = 113 - s.
        let mut s = 0i32;
        let mut m = mant;
        while (m & 0x400) == 0 {
            m <<= 1;
            s += 1;
        }
        m &= 0x3ff;
        let exp32 = (113 - s) as u32;
        return f32::from_bits(sign | (exp32 << 23) | (m << 13));
    }
    if exp == 0x1f {
        return f32::from_bits(sign | 0x7f80_0000 | (mant << 13));
    }
    let exp32 = exp + 127 - 15;
    f32::from_bits(sign | (exp32 << 23) | (mant << 13))
}

// ─────────────────────────── bit io ───────────────────────────

#[derive(Default)]
pub struct BitWriter {
    pub bytes: Vec<u8>,
    acc: u64,
    nbits: u32,
}

impl BitWriter {
    pub fn put(&mut self, value: u32, bits: u32) {
        debug_assert!(bits <= 32);
        debug_assert!(bits == 32 || value < (1u32 << bits), "valor {value} não cabe em {bits} bits");
        self.acc |= (value as u64) << self.nbits;
        self.nbits += bits;
        while self.nbits >= 8 {
            self.bytes.push((self.acc & 0xff) as u8);
            self.acc >>= 8;
            self.nbits -= 8;
        }
    }
    pub fn finish(mut self) -> Vec<u8> {
        if self.nbits > 0 {
            self.bytes.push((self.acc & 0xff) as u8);
        }
        self.bytes
    }
}

pub struct BitReader<'a> {
    data: &'a [u8],
    pos: usize,
    acc: u64,
    nbits: u32,
}

impl<'a> BitReader<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        BitReader { data, pos: 0, acc: 0, nbits: 0 }
    }
    pub fn get(&mut self, bits: u32) -> u32 {
        while self.nbits < bits {
            let b = *self.data.get(self.pos).unwrap_or(&0) as u64;
            self.pos += 1;
            self.acc |= b << self.nbits;
            self.nbits += 8;
        }
        let mask = if bits == 32 { u32::MAX as u64 } else { (1u64 << bits) - 1 };
        let v = (self.acc & mask) as u32;
        self.acc >>= bits;
        self.nbits -= bits;
        v
    }
}

// ─────────────────────────── quantizadores de campo ───────────────────────────

fn q_lin(v: f32, lo: f32, hi: f32, bits: u32) -> u32 {
    if bits == 32 {
        return v.to_bits();
    }
    if bits == 16 {
        return f32_to_f16(v) as u32;
    }
    let levels = ((1u64 << bits) - 1) as f32;
    (((v - lo) / (hi - lo).max(1e-9)).clamp(0.0, 1.0) * levels).round() as u32
}

fn dq_lin(q: u32, lo: f32, hi: f32, bits: u32) -> f32 {
    if bits == 32 {
        return f32::from_bits(q);
    }
    if bits == 16 {
        return f16_to_f32(q as u16);
    }
    let levels = ((1u64 << bits) - 1) as f32;
    lo + (q as f32 / levels) * (hi - lo)
}

fn q_log(v: f32, lo: f32, hi: f32, bits: u32) -> u32 {
    if bits >= 16 {
        return q_lin(v, lo, hi, bits);
    }
    q_lin(v.max(1e-9).ln(), lo.ln(), hi.ln(), bits)
}

fn dq_log(q: u32, lo: f32, hi: f32, bits: u32) -> f32 {
    if bits >= 16 {
        return dq_lin(q, lo, hi, bits);
    }
    dq_lin(q, lo.ln(), hi.ln(), bits).exp()
}

// ─────────────────────────── átomos codificáveis ───────────────────────────

/// Contrato: como um átomo vira bits (e volta). Cada modalidade implementa.
pub trait Codable: Sized {
    fn encode(&self, amp: f32, w: &mut BitWriter, p: Profile, dims: &[usize], b: &AtomBounds);
    fn decode(r: &mut BitReader, p: Profile, dims: &[usize], b: &AtomBounds) -> (Self, f32);
    fn render_unit_dyn(&self, dims: &[usize]) -> Patch;
}

const PI: f32 = std::f32::consts::PI;
const TAU: f32 = std::f32::consts::TAU;

impl Codable for Gabor2D {
    fn encode(&self, amp: f32, w: &mut BitWriter, p: Profile, dims: &[usize], b: &AtomBounds) {
        let (pb, fb, ab, sb, mb) = p.bits();
        w.put(q_lin(self.cy, 0.0, dims[0] as f32, pb), pb);
        w.put(q_lin(self.cx, 0.0, dims[1] as f32, pb), pb);
        w.put(q_log(self.freq, b.freq_min, b.freq_max, fb), fb);
        w.put(q_lin(self.theta.rem_euclid(PI), 0.0, PI, ab), ab);
        w.put(q_lin(self.phase.rem_euclid(TAU), 0.0, TAU, ab), ab);
        w.put(q_log(self.sigma, b.sigma_min, b.sigma_max, sb), sb);
        w.put(q_lin(amp, -AMP_RANGE, AMP_RANGE, mb), mb);
    }

    fn decode(r: &mut BitReader, p: Profile, dims: &[usize], b: &AtomBounds) -> (Self, f32) {
        let (pb, fb, ab, sb, mb) = p.bits();
        let cy = dq_lin(r.get(pb), 0.0, dims[0] as f32, pb);
        let cx = dq_lin(r.get(pb), 0.0, dims[1] as f32, pb);
        let freq = dq_log(r.get(fb), b.freq_min, b.freq_max, fb);
        let theta = dq_lin(r.get(ab), 0.0, PI, ab);
        let phase = dq_lin(r.get(ab), 0.0, TAU, ab);
        let sigma = dq_log(r.get(sb), b.sigma_min, b.sigma_max, sb);
        let amp = dq_lin(r.get(mb), -AMP_RANGE, AMP_RANGE, mb);
        (Gabor2D { cx, cy, freq, theta, phase, sigma }, amp)
    }

    fn render_unit_dyn(&self, dims: &[usize]) -> Patch {
        use agua_atoms::Atom;
        self.render_unit(dims)
    }
}

impl Codable for Gabor1D {
    fn encode(&self, amp: f32, w: &mut BitWriter, p: Profile, dims: &[usize], b: &AtomBounds) {
        let (pb, fb, ab, sb, mb) = p.bits();
        // áudio tem muitas amostras: posição precisa de mais bits
        let pb = if pb < 32 && pb < 20 { pb.max(16) } else { pb };
        w.put(q_lin(self.c, 0.0, dims[0] as f32, pb), pb);
        w.put(q_log(self.freq, b.freq_min, b.freq_max, fb), fb);
        w.put(q_lin(self.phase.rem_euclid(TAU), 0.0, TAU, ab), ab);
        w.put(q_log(self.sigma, b.sigma_min, b.sigma_max, sb), sb);
        w.put(q_lin(amp, -AMP_RANGE, AMP_RANGE, mb), mb);
    }

    fn decode(r: &mut BitReader, p: Profile, dims: &[usize], b: &AtomBounds) -> (Self, f32) {
        let (pb, fb, ab, sb, mb) = p.bits();
        let pb = if pb < 32 && pb < 20 { pb.max(16) } else { pb };
        let c = dq_lin(r.get(pb), 0.0, dims[0] as f32, pb);
        let freq = dq_log(r.get(fb), b.freq_min, b.freq_max, fb);
        let phase = dq_lin(r.get(ab), 0.0, TAU, ab);
        let sigma = dq_log(r.get(sb), b.sigma_min, b.sigma_max, sb);
        let amp = dq_lin(r.get(mb), -AMP_RANGE, AMP_RANGE, mb);
        (Gabor1D { c, freq, phase, sigma }, amp)
    }

    fn render_unit_dyn(&self, dims: &[usize]) -> Patch {
        use agua_atoms::Atom;
        self.render_unit(dims)
    }
}

/// Amplitudes vivem em [-4, 4] na prática (sinais em [-1,1] zero-médios).
const AMP_RANGE: f32 = 4.0;

// ─────────────────────────── arquivo .agua ───────────────────────────

pub struct Encoded {
    pub bytes: Vec<u8>,
    pub atoms: usize,
}

pub fn encode<A: Codable>(
    atoms: &[(A, f32)],
    dims: &[usize],
    dc: f32,
    bounds: &AtomBounds,
    profile: Profile,
) -> Encoded {
    let mut head = Vec::new();
    head.extend_from_slice(MAGIC);
    head.push(VERSION);
    head.push(profile.tag());
    head.push(dims.len() as u8);
    for d in dims {
        head.extend_from_slice(&(*d as u32).to_le_bytes());
    }
    head.extend_from_slice(&f32_to_f16(dc).to_le_bytes());
    head.extend_from_slice(&(atoms.len() as u32).to_le_bytes());
    // bounds fazem parte do arquivo (o decodificador precisa das escalas)
    for v in [bounds.freq_min, bounds.freq_max, bounds.sigma_min, bounds.sigma_max] {
        head.extend_from_slice(&v.to_le_bytes());
    }

    let mut w = BitWriter::default();
    for (a, amp) in atoms {
        a.encode(*amp, &mut w, profile, dims, bounds);
    }
    let mut bytes = head;
    bytes.extend_from_slice(&w.finish());
    Encoded { bytes, atoms: atoms.len() }
}

pub struct Decoded<A> {
    pub atoms: Vec<(A, f32)>,
    pub dims: Vec<usize>,
    pub dc: f32,
}

pub fn decode<A: Codable>(bytes: &[u8]) -> Result<Decoded<A>, String> {
    if bytes.len() < 8 || &bytes[0..4] != MAGIC {
        return Err("não é arquivo .agua".into());
    }
    if bytes[4] != VERSION {
        return Err(format!("versão {} desconhecida", bytes[4]));
    }
    let profile = Profile::from_tag(bytes[5]).ok_or("perfil inválido")?;
    let nd = bytes[6] as usize;
    let mut p = 7;
    let mut dims = Vec::with_capacity(nd);
    for _ in 0..nd {
        dims.push(u32::from_le_bytes(bytes[p..p + 4].try_into().unwrap()) as usize);
        p += 4;
    }
    let dc = f16_to_f32(u16::from_le_bytes(bytes[p..p + 2].try_into().unwrap()));
    p += 2;
    let n = u32::from_le_bytes(bytes[p..p + 4].try_into().unwrap()) as usize;
    p += 4;
    let mut bv = [0f32; 4];
    for b in bv.iter_mut() {
        *b = f32::from_le_bytes(bytes[p..p + 4].try_into().unwrap());
        p += 4;
    }
    let bounds =
        AtomBounds { freq_min: bv[0], freq_max: bv[1], sigma_min: bv[2], sigma_max: bv[3] };

    let mut r = BitReader::new(&bytes[p..]);
    let mut atoms = Vec::with_capacity(n);
    for _ in 0..n {
        atoms.push(A::decode(&mut r, profile, &dims, &bounds));
    }
    Ok(Decoded { atoms, dims, dc })
}

/// Reconstrói o campo (com DC) a partir dos átomos decodificados.
pub fn render<A: Codable>(d: &Decoded<A>) -> Field {
    let mut f = Field::zeros(&d.dims);
    for (a, amp) in &d.atoms {
        let unit = a.render_unit_dyn(&d.dims);
        f.add_scaled_patch(&unit, *amp);
    }
    for v in f.data.iter_mut() {
        *v += d.dc;
    }
    f
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fp16_roundtrip() {
        for v in [0.0f32, 1.0, -1.0, 0.5, 3.14159, -0.001, 65504.0, 1e-5] {
            let back = f16_to_f32(f32_to_f16(v));
            let err = (back - v).abs() / v.abs().max(1e-6);
            assert!(err < 1e-2, "fp16 {v} → {back}");
        }
        assert_eq!(f16_to_f32(f32_to_f16(0.0)), 0.0);
    }

    #[test]
    fn bit_io_roundtrip() {
        let mut w = BitWriter::default();
        let vals = [(5u32, 3u32), (1000, 12), (63, 6), (7, 16), (0, 8)];
        for (v, b) in vals {
            w.put(v, b);
        }
        let bytes = w.finish();
        let mut r = BitReader::new(&bytes);
        for (v, b) in vals {
            assert_eq!(r.get(b), v, "bits={b}");
        }
    }

    #[test]
    fn atom_roundtrip_fp16_precise() {
        let bounds =
            AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
        let g = Gabor2D { cx: 33.7, cy: 12.2, freq: 0.08, theta: 1.1, phase: 2.0, sigma: 5.5 };
        let enc = encode(&[(g.clone(), 0.42)], &[64, 64], 0.5, &bounds, Profile::Fp16);
        let dec: Decoded<Gabor2D> = decode(&enc.bytes).unwrap();
        assert_eq!(dec.atoms.len(), 1);
        let (g2, amp) = &dec.atoms[0];
        assert!((g2.cx - g.cx).abs() < 0.05, "cx {} vs {}", g2.cx, g.cx);
        assert!((g2.freq - g.freq).abs() < 1e-3);
        assert!((amp - 0.42).abs() < 1e-3);
        assert!((dec.dc - 0.5).abs() < 1e-3);
    }

    #[test]
    fn bytes_per_atom_shrinks() {
        let bounds =
            AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
        let g = Gabor2D { cx: 1.0, cy: 1.0, freq: 0.1, theta: 0.0, phase: 0.0, sigma: 3.0 };
        let atoms: Vec<_> = (0..100).map(|_| (g.clone(), 0.3)).collect();
        let sizes: Vec<usize> = [Profile::F32, Profile::Fp16, Profile::Q12, Profile::Q8]
            .iter()
            .map(|p| encode(&atoms, &[64, 64], 0.0, &bounds, *p).bytes.len())
            .collect();
        // header ~31 B; o que importa é a tendência estritamente decrescente
        assert!(sizes[0] > sizes[1] && sizes[1] > sizes[2] && sizes[2] > sizes[3], "{sizes:?}");
        // Gabor2D em Q8: 10+10+6+6+6+6+10 = 54 bits = 6.75 B/átomo
        let per_atom_q8 = (sizes[3] - 31) as f32 / 100.0;
        assert!(per_atom_q8 < 7.0, "q8 = {per_atom_q8} B/átomo");
        let per_atom_f32 = (sizes[0] - 31) as f32 / 100.0;
        assert!(per_atom_f32 > 27.0, "f32 = {per_atom_f32} B/átomo");
    }

    #[test]
    fn q8_decodes_close_enough() {
        let bounds =
            AtomBounds { freq_min: 0.015, freq_max: 0.45, sigma_min: 1.2, sigma_max: 24.0 };
        let g = Gabor2D { cx: 20.0, cy: 40.0, freq: 0.2, theta: 0.7, phase: 3.0, sigma: 4.0 };
        let enc = encode(&[(g.clone(), -0.9)], &[64, 64], 0.0, &bounds, Profile::Q8);
        let dec: Decoded<Gabor2D> = decode(&enc.bytes).unwrap();
        let (g2, amp) = &dec.atoms[0];
        assert!((g2.cx - g.cx).abs() < 0.2, "cx {}", g2.cx);
        assert!((g2.cy - g.cy).abs() < 0.2);
        assert!((amp + 0.9).abs() < 0.02);
    }

    #[test]
    fn rejects_garbage() {
        let e: Result<Decoded<Gabor2D>, _> = decode(b"NOPE1234");
        assert!(e.is_err());
    }
}
