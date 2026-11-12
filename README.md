# iara-agua

A population of Gabor "beings" lives on the residual of a signal: each one eats the part of the signal it
explains and pays a metabolic cost to exist. Beings that explain less than they cost die. The sum of the
survivors is the reconstruction, and the number of atoms is not chosen by anyone; it emerges from the cost.

Headless, deterministic Rust workspace (CPU only, no GPU). Every run writes a manifest with the git
revision, a `metrics.jsonl` and PNG/WAV outputs.

## Highlights

- **Hyperparameter-free denoising.** With metabolism set by the law `met = 40·σ̂²` (σ̂ estimated from the
  noisy signal by MAD), survival beats matching pursuit run at the *oracle* atom count in **737/737** cases
  of an 8-hour battery (**+2.73 dB** mean) and the classic VisuShrink Haar baseline in **81/81** (+3.84 dB).
  It does *not* beat a Gaussian filter tuned with the oracle σ: 340/713.
- **A measured metabolism law.** The critical metabolism follows `met_c ≈ σ²·(1 + 2·dof)`: one σ² for an
  atom that does not search, about 2σ² per degree of freedom, which is the AIC penalty. Fit on the 8-hour
  battery: `c = 2.253·dof + 0.854`, R² = 0.984, n = 990. The slope saturates at ~2.2 rather than exactly 2,
  and it depends on the search budget.
- **Discrepancy principle emerges on its own.** The population stops when `‖r‖²/n ÷ σ²` = 0.99 to 1.03
  (σ = 0.10 to 0.30, across image and audio) without being told σ.
- **Inpainting by routing between experts.** A smooth expert and a texture expert compete per region. In
  6 test scenes the result is never below Laplace diffusion and matches the per-scene oracle in 5 of 6.
- **4679 law validations, 0 errors** in the 8-hour battery (`scripts/battery8h.py`).

## How it works

```
signal (image / audio / video)
   |
   v
 loop:  spawn a being where the residual is large   ("rain" on hotspots)
        it eats the patch it covers                  (amplitude in closed form)
        it adjusts its geometry                      (analytic gradient or mutation)
        it dies if it explains < its metabolic cost  (met = k·σ²)
        stop when the residual reaches the noise level
   |
   v
 survivors = the code; their count N emerges
```

The same loop runs on any atom type behind `trait Atom` (1D/2D/3D Gabor, blob, edge, a non-linear
occlusion atom), so audio, images and video share one engine (`World<A>::step`).

## Results

Numbers from `docs/` and the project result log. "MP" = matching pursuit.

| test | result | baseline |
|---|---|---|
| zero-shot denoise (`denoise --auto`) | wins 737/737, +2.73 dB | MP at oracle N |
| same | wins 81/81, +3.84 dB | VisuShrink (Haar) |
| metabolism as a rate-distortion dial | met 0.005 to 1.0 gives 3000 to 129 atoms, 41.7 to 28.8 dB, monotonic | n/a |
| law 2, `met_c` vs dof | `c = 2.253·dof + 0.854`, R² 0.984, n = 990 | AIC predicts slope 2 |
| law 4, Morozov | `‖r‖²/n÷σ²` = 0.99 to 1.03 | 1.00 |
| law 6, predation phase transition | jump at θ* ≈ 0.50 (battery: 0.48 to 0.56) | n/a |
| polyglot routing (Gabor / blob / edge) | +5.53 dB on edges, +1.57 smooth, +0.40 texture | single Gabor dictionary |
| inpaint, routed experts + best-of-N | photo_texture 24.21 dB (= oracle), photo_smooth 33.38 dB | Laplace diffusion, fixed expert |
| non-linear occlusion atom | 6/9 scenes, +0.69 dB | best linear atom |
| synapse = Gram matrix `G_ij = ⟨ψi,ψj⟩` | +1.0 at d=0, −0.84 at d=4, +0.51 at d=8 (Mexican hat); susceptibility 74 → 120289 as coherence μ goes 0.15 → 0.60 | n/a |

## What didn't work (and why)

- **Survival loses to greedy pursuit at equal compute.** In 192 runs, matching pursuit with the same budget
  wins in every modality (synthetic −4.6 dB, photo −5.1, speech −12.9, trumpet −26.7). MP picks the best
  of C candidates against the updated residual; survival spends evaluations on rejected mutations. An
  early "beats MP on audio and video" result was an MP run with ~30× less budget.
- **It is not a compressor.** With quantized parameters (Q8) the codec needs 6633 B for 28.1 dB; JPEG q30
  needs 1259 B for 31.1 dB. JPEG is 5 to 18× better along the whole curve.
- **The pruning law does not transfer to LLM perplexity.** Metabolic pruning of a 40M byte-level LM's MLP
  neurons gives 1.604 bpb at 30% sparsity vs Wanda's 1.192 (full model 1.122), and loses at 50% and 70%
  too. It had won on per-layer reconstruction; that metric was a misleading proxy because errors compound
  across layers.
- **Pure selection cannot inpaint what is not observed.** Without a prior, filling a smooth hole gives
  −0.28 dB and a texture hole 6.08 dB, vs Laplace diffusion 18.72 / 15.41. Selection needs signal to eat.
- **Predation does not buy quality at fixed N.** +0.29 dB at N = 120 did not replicate (−2.72 dB, 0/3 on
  another target). It only makes the code use fewer atoms (53 vs 90).

Nine of my own hypotheses were falsified along the way (universal γ, Donoho-Johnstone scaling, continued
fractions for quantization, imposed Morozov replacing k, and others); each is written up in `docs/` and
`JOURNAL.md`.

## Crates

```
agua-core       N-dim Field, sparse Patch, deterministic RNG, hotspot sampler
agua-atoms      trait Atom + Gabor 1D/2D/3D, analytic gradient (gradcheck: cosine 100/100 image, 99/100 audio)
agua-life       the loop, generic: World<A>::step (spawn / adjust / die / predate)
agua-metrics    baselines: random, matching pursuit, JPEG
agua-codec      population -> bytes (.agua): fp16, Q12, Q8
agua-telemetry  run manifest + git rev, metrics.jsonl, PNG/WAV export
agua-lab        CLI harness (subcommands below)
```

## Build & run

Requires a stable Rust toolchain. No GPU, no special environment.

```bash
cargo test --release                       # 36 unit tests
cargo run --release --bin agua-lab -- demo                 # image (also: demo-audio, demo-video)
cargo run --release --bin agua-lab -- run --config configs/exp001_image_survival.toml
cargo run --release --bin agua-lab -- denoise --target synthetic --sigma 0.2 --auto
cargo run --release --bin agua-lab -- law2                 # metabolism law (also law3 .. law9)
cargo run --release --bin agua-lab -- polyglot             # emergent routing between atom types
cargo run --release --bin agua-lab -- inpaint --poliglota  # routed experts (--prior for the VQ prior)
cargo run --release --bin agua-lab -- synapse              # Gram-matrix synapse laws
```

Other subcommands: `compare`, `generate`, `occlusion`, `resonance`, `probe`, `law`. Run outputs go to
`runs/` (not committed). The Python files in `scripts/` only orchestrate sweeps and summarize
results (`battery8h.py`, `analise8h.py`, `lei_garimpo.py`).

## Docs

The write-ups are in Portuguese.

- `docs/LAWS.md`, `docs/LAW_METABOLISM.md`: the laws, each with its prediction recorded before the run
- `docs/LEI_DO_GARIMPO.md`: the "prospecting" penalty `2σ²·ln M` (RIC): finance, variable selection, pruning
- `docs/LEIS_SINAPSE.md`: Gram matrix as synapse, criticality
- `docs/IARA_INFERENCIA.md`, `docs/PONTE_LINGUAS.md`: inpainting and routing between representations
- `docs/EXPERIMENT_PROTOCOL.md`: the honesty rules (multi-seed, iso-budget baselines, controls)
- `docs/RESULTS_F1.md`, `docs/F2_RESSONANCIA.md`, `docs/F4_GERACAO.md`, `docs/F5_ESCALA.md`: phase reports
- `JOURNAL.md`: chronological lab notebook

## Status

Phases F0 to F6 of `docs/ROADMAP.md` are closed. This is a research codebase. Its contribution is a set of
measured laws (metabolism ≈ AIC, prospecting ≈ RIC, synapse = Gram matrix) and two working demos (denoising
and inpainting by routing). It is not a state-of-the-art compressor or LLM pruner.

Author: [leonardoosouza1](https://github.com/leonardoosouza1)

License: TBD (`Cargo.toml` declares MIT, but no LICENSE file is committed yet).
