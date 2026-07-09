//! agua-life — a física da sobrevivência: seres escoando como água sobre o resíduo.
//!
//! GENÉRICA sobre o átomo (`World<A: Atom>`): a MESMA física roda imagem
//! (Gabor2D), áudio (Gabor1D), vídeo (Gabor3D) e qualquer modalidade futura.
//! Este crate não sabe o que é pixel, amostra ou frame — só Field/Patch.
//!
//! Mundo = sinal-alvo S. Reconstrução = Σ (amp_i · ψ_i) dos seres vivos.
//! Resíduo r = S − reconstrução. Comida do ser = energia que o SEU átomo
//! explica do resíduo-sem-ele (utilidade marginal). Sem comida → morre e o
//! patch dele volta pro resíduo (vira comida pros outros). Onde sobra
//! resíduo, chove ser novo e os bem-alimentados reproduzem com mutação.
//!
//! O metabolismo é o dial taxa-distorção: um ser só sobrevive se explica
//! mais energia que o custo de existir → `metabolism/food_gain` é o piso
//! de utilidade por átomo. Metabolismo alto = codificação enxuta e grossa;
//! baixo = população grande e detalhista.

use agua_atoms::{fit_amp, Atom, AtomBounds};
use agua_core::{psnr, unflatten, Field, HotspotSampler, Patch, Rng};
use serde::{Deserialize, Serialize};

/// Parâmetros da dinâmica de vida (vem do TOML do experimento).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LifeConfig {
    /// população inicial (primeira chuva)
    pub init: usize,
    /// teto de população (competição por vaga)
    pub cap: usize,
    /// custo de existir por tick (o dial taxa-distorção)
    pub metabolism: f32,
    /// conversão energia-explicada → comida
    pub food_gain: f32,
    /// energia ao nascer
    pub initial_energy: f32,
    /// teto de energia (evita imortal rico)
    pub max_energy: f32,
    /// energia mínima pra reproduzir
    pub repro_threshold: f32,
    /// custo pago pelo pai ao reproduzir
    pub repro_cost: f32,
    /// intensidade da mutação (filhos e ajuste-de-vida)
    pub mutation_sigma: f32,
    /// ajuste-de-vida ligado? (busca local na geometria durante a vida)
    pub life_tuning: bool,
    /// tenta ajustar geometria a cada N ticks
    pub tune_every: u64,
    /// chuva de imigrantes a cada N ticks
    pub rain_every: u64,
    /// quantos imigrantes por chuva
    pub rain_count: usize,
}

impl Default for LifeConfig {
    fn default() -> Self {
        LifeConfig {
            init: 48,
            cap: 400,
            metabolism: 0.02,
            food_gain: 1.0,
            initial_energy: 2.0,
            max_energy: 10.0,
            repro_threshold: 8.0,
            repro_cost: 3.0,
            mutation_sigma: 0.25,
            life_tuning: true,
            tune_every: 3,
            rain_every: 20,
            rain_count: 6,
        }
    }
}

/// Um ser vivo: átomo + patch unitário cacheado + amplitude + energia.
#[derive(Debug)]
pub struct Being<A: Atom> {
    pub id: u64,
    pub atom: A,
    pub unit: Patch,
    /// ||ψ||² cacheado
    pub n2: f32,
    pub amp: f32,
    pub energy: f32,
    pub age: u64,
}

/// Snapshot serializável da população — a própria CODIFICAÇÃO do sinal.
#[derive(Serialize)]
pub struct BeingRecord<A: Atom> {
    pub id: u64,
    pub atom: A,
    pub amp: f32,
    pub energy: f32,
    pub age: u64,
}

pub struct World<A: Atom> {
    pub target: Field,
    pub recon: Field,
    pub residual: Field,
    pub beings: Vec<Being<A>>,
    pub rng: Rng,
    pub tick: u64,
    pub cfg: LifeConfig,
    pub bounds: AtomBounds,
    pub births: u64,
    pub deaths: u64,
    next_id: u64,
}

impl<A: Atom> World<A> {
    /// `target` deve estar zero-médio (o DC é side-info de 4 bytes; átomos
    /// oscilatórios têm média ~0). PSNR é idêntico com ou sem o DC.
    pub fn new(target: Field, cfg: LifeConfig, bounds: AtomBounds, seed: u64) -> Self {
        let mut world = World {
            recon: Field::zeros(&target.dims),
            residual: target.clone(),
            target,
            beings: Vec::new(),
            rng: Rng::new(seed),
            tick: 0,
            cfg,
            bounds,
            births: 0,
            deaths: 0,
            next_id: 0,
        };
        world.rain(world.cfg.init);
        world
    }

    /// Chuva: nascem seres em posições ∝ resíduo² (a água cai onde falta).
    fn rain(&mut self, count: usize) {
        if count == 0 || self.beings.len() >= self.cfg.cap {
            return;
        }
        let sampler = HotspotSampler::new(&self.residual);
        for _ in 0..count {
            if self.beings.len() >= self.cfg.cap {
                return;
            }
            let idx = sampler.sample(&mut self.rng);
            let pos: Vec<f32> = unflatten(&self.target.dims, idx)
                .into_iter()
                .map(|v| v as f32)
                .collect();
            let atom = A::rain_at(&pos, &self.bounds, &mut self.rng);
            self.spawn(atom);
        }
    }

    /// Nasce um ser: amp ajustada em forma fechada contra o resíduo atual.
    fn spawn(&mut self, atom: A) -> bool {
        let unit = atom.render_unit(&self.target.dims);
        let n2 = unit.norm2();
        if n2 < 1e-4 {
            return false;
        }
        let Some((amp, _)) = fit_amp(&unit, &self.residual) else {
            return false;
        };
        self.recon.add_scaled_patch(&unit, amp);
        self.residual.add_scaled_patch(&unit, -amp);
        self.next_id += 1;
        self.births += 1;
        self.beings.push(Being {
            id: self.next_id,
            atom,
            unit,
            n2,
            amp,
            energy: self.cfg.initial_energy,
            age: 0,
        });
        true
    }

    /// Remove o ser i: o patch dele volta pro resíduo (vira comida).
    fn kill(&mut self, i: usize) {
        let b = self.beings.swap_remove(i);
        self.recon.add_scaled_patch(&b.unit, -b.amp);
        self.residual.add_scaled_patch(&b.unit, b.amp);
        self.deaths += 1;
    }

    /// Busca local sem gradiente: propõe geometria mutada; troca se explicar
    /// mais o resíduo-sem-ele (com histerese de 2% pra não thrashar).
    fn try_tune(&mut self, i: usize) {
        let (cand, u_cur, amp_old) = {
            let b = &self.beings[i];
            let d_cur = b.unit.dot_field(&self.residual) + b.amp * b.n2;
            let u_cur = d_cur * d_cur / b.n2;
            let cand = b.atom.mutate(
                &mut self.rng,
                self.cfg.mutation_sigma,
                &self.bounds,
                &self.target.dims,
            );
            (cand, u_cur, b.amp)
        };
        let cu = cand.render_unit(&self.target.dims);
        let cn2 = cu.norm2();
        if cn2 < 1e-4 {
            return;
        }
        // ⟨r_sem_eu, ψc⟩ = ⟨residual, ψc⟩ + amp·⟨ψ_velho, ψc⟩ (overlap dos patches)
        let cd = cu.dot_field(&self.residual) + amp_old * cu.dot_patch(&self.beings[i].unit);
        let u_new = cd * cd / cn2;
        if u_new > u_cur * 1.02 {
            // remove o corpo velho, veste o novo, re-ajusta amp
            let old_unit = std::mem::replace(&mut self.beings[i].unit, cu);
            self.recon.add_scaled_patch(&old_unit, -amp_old);
            self.residual.add_scaled_patch(&old_unit, amp_old);
            let new_amp = cd / cn2;
            let b = &mut self.beings[i];
            b.atom = cand;
            b.n2 = cn2;
            b.amp = new_amp;
            let unit = b.unit.clone();
            self.recon.add_scaled_patch(&unit, new_amp);
            self.residual.add_scaled_patch(&unit, -new_amp);
        }
    }

    /// Um tick de vida: ajustar → comer → morrer → reproduzir → chover.
    pub fn step(&mut self) {
        self.tick += 1;
        let n = self.beings.len();

        // 1. ajuste-de-vida (geometria, busca local)
        if self.cfg.life_tuning && self.tick % self.cfg.tune_every == 0 {
            let order = self.rng.shuffled_indices(n);
            for i in order {
                if i < self.beings.len() {
                    self.try_tune(i);
                }
            }
        }

        // 2. comer: re-ajusta amp (relaxação Gauss-Seidel — a harmonização com
        //    os vizinhos) e ganha comida ∝ energia explicada
        let order = self.rng.shuffled_indices(self.beings.len());
        for i in order {
            let b = &mut self.beings[i];
            let d = b.unit.dot_field(&self.residual) + b.amp * b.n2;
            let new_amp = d / b.n2;
            let delta = new_amp - b.amp;
            if delta.abs() > 1e-7 {
                let unit = b.unit.clone();
                b.amp = new_amp;
                self.recon.add_scaled_patch(&unit, delta);
                self.residual.add_scaled_patch(&unit, -delta);
            }
            let b = &mut self.beings[i];
            let explained = d * d / b.n2;
            b.energy = (b.energy + self.cfg.food_gain * explained - self.cfg.metabolism)
                .min(self.cfg.max_energy);
            b.age += 1;
        }

        // 3. morte (fome): varre de trás pra frente por causa do swap_remove
        let mut i = self.beings.len();
        while i > 0 {
            i -= 1;
            if self.beings[i].energy <= 0.0 {
                self.kill(i);
            }
        }

        // 4. reprodução: pais ricos geram filhos mutados perto de si
        let parents: Vec<usize> = (0..self.beings.len())
            .filter(|&i| self.beings[i].energy >= self.cfg.repro_threshold)
            .collect();
        for i in parents {
            if self.beings.len() >= self.cfg.cap {
                break;
            }
            self.beings[i].energy -= self.cfg.repro_cost;
            let child = self.beings[i].atom.mutate(
                &mut self.rng,
                self.cfg.mutation_sigma,
                &self.bounds,
                &self.target.dims,
            );
            self.spawn(child);
        }

        // 5. chuva de imigrantes nos hotspots do resíduo
        if self.cfg.rain_every > 0 && self.tick % self.cfg.rain_every == 0 {
            self.rain(self.cfg.rain_count);
        }

        // 6. anti-deriva float: reconstrói recon/residual do zero periodicamente
        if self.tick % 500 == 0 {
            self.rebuild();
        }
    }

    /// Recalcula recon = Σ patches e residual = target − recon (mata deriva f32).
    pub fn rebuild(&mut self) {
        let mut recon = Field::zeros(&self.target.dims);
        for b in &self.beings {
            recon.add_scaled_patch(&b.unit, b.amp);
        }
        let mut residual = self.target.clone();
        for (r, c) in residual.data.iter_mut().zip(&recon.data) {
            *r -= c;
        }
        self.recon = recon;
        self.residual = residual;
    }

    /// Máx |target − recon − residual| — invariante de consistência.
    pub fn invariant_error(&self) -> f32 {
        self.target
            .data
            .iter()
            .zip(&self.recon.data)
            .zip(&self.residual.data)
            .map(|((t, c), r)| (t - c - r).abs())
            .fold(0.0f32, f32::max)
    }

    pub fn psnr(&self, peak: f32) -> f32 {
        let mse = self.residual.energy() / self.target.len() as f32;
        psnr(mse, peak)
    }

    /// Custo em bytes da codificação atual (params f32 por ser).
    pub fn bytes(&self) -> usize {
        self.beings.len() * A::PARAMS * 4
    }

    /// Histograma da população por banda de frequência característica.
    /// Mede se a pirâmide grosso→fino emerge sozinha.
    pub fn freq_bands(&self, edges: &[f32]) -> Vec<usize> {
        let mut bands = vec![0usize; edges.len() + 1];
        for b in &self.beings {
            let f = b.atom.freq();
            let mut k = edges.len();
            for (i, e) in edges.iter().enumerate() {
                if f < *e {
                    k = i;
                    break;
                }
            }
            bands[k] += 1;
        }
        bands
    }

    pub fn population_records(&self) -> Vec<BeingRecord<A>> {
        self.beings
            .iter()
            .map(|b| BeingRecord {
                id: b.id,
                atom: b.atom.clone(),
                amp: b.amp,
                energy: b.energy,
                age: b.age,
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use agua_atoms::{Gabor1D, Gabor2D};

    fn synthetic_2d(size: usize) -> Field {
        // alvo multi-escala zero-médio: blob suave + listras finas
        let s = size as f32;
        let mut f = Field::from_fn_2d(size, size, |x, y| {
            let (xf, yf) = (x as f32, y as f32);
            let blob = 0.6
                * (-((xf - s * 0.3).powi(2) + (yf - s * 0.3).powi(2)) / (2.0 * (s * 0.15).powi(2)))
                    .exp();
            let stripes = if xf > s * 0.5 && yf > s * 0.5 {
                0.3 * (std::f32::consts::TAU * 0.2 * xf).sin()
            } else {
                0.0
            };
            blob + stripes
        });
        let m = f.mean();
        for v in f.data.iter_mut() {
            *v -= m;
        }
        f
    }

    fn synthetic_1d(n: usize) -> Field {
        let mut f = Field::from_fn_1d(n, |i| {
            let t = i as f32;
            let low = 0.4 * (std::f32::consts::TAU * 0.004 * t).sin();
            let burst = if i > n / 3 && i < n / 2 {
                0.3 * (std::f32::consts::TAU * 0.08 * t).sin()
            } else {
                0.0
            };
            low + burst
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
    fn invariant_holds_over_ticks() {
        let mut w: World<Gabor2D> =
            World::new(synthetic_2d(32), LifeConfig::default(), bounds(), 11);
        for _ in 0..60 {
            w.step();
        }
        assert!(w.invariant_error() < 1e-3, "deriva: {}", w.invariant_error());
    }

    #[test]
    fn survival_improves_psnr_image() {
        let mut w: World<Gabor2D> =
            World::new(synthetic_2d(48), LifeConfig::default(), bounds(), 22);
        let p0 = w.psnr(1.0);
        for _ in 0..300 {
            w.step();
        }
        let p1 = w.psnr(1.0);
        assert!(w.beings.len() > 4, "população colapsou: {}", w.beings.len());
        assert!(p1 > p0 + 1.0, "psnr não subiu: {p0:.2} → {p1:.2}");
    }

    /// O teste do AGNOSTICISMO: a MESMA física resolve um sinal 1D (áudio).
    #[test]
    fn survival_improves_psnr_audio() {
        let mut w: World<Gabor1D> =
            World::new(synthetic_1d(2048), LifeConfig::default(), bounds(), 33);
        let p0 = w.psnr(1.0);
        for _ in 0..300 {
            w.step();
        }
        let p1 = w.psnr(1.0);
        assert!(w.beings.len() > 4, "população colapsou: {}", w.beings.len());
        assert!(p1 > p0 + 1.0, "psnr não subiu: {p0:.2} → {p1:.2}");
        assert!(w.invariant_error() < 1e-3);
    }

    #[test]
    fn deterministic_same_seed() {
        let mk = || {
            let mut w: World<Gabor2D> =
                World::new(synthetic_2d(32), LifeConfig::default(), bounds(), 77);
            for _ in 0..100 {
                w.step();
            }
            (w.psnr(1.0), w.beings.len(), w.births, w.deaths)
        };
        assert_eq!(mk(), mk());
    }
}
