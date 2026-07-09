# iara-água — solver universal por sobrevivência

> **A IARA madura:** joga-se uma população de seres-átomos sobre um sinal (imagem, áudio, texto…); tentando **sobreviver** — comida = resíduo local — eles se auto-organizam **como água achando o caminho** e a soma dos vivos **É** a reconstrução/geração. Visão completa: `../IARA_VISAO_SERES_COMO_AGUA.md`.

## Quickstart

```bash
cargo test --release                                        # 23 testes
cargo run --release --bin agua-lab -- demo                  # imagem ~4s
cargo run --release --bin agua-lab -- demo-audio            # áudio (gera .wav ouvível)
cargo run --release --bin agua-lab -- demo-video            # vídeo ~16s
cargo run --release --bin agua-lab -- run \
    --config configs/exp001_image_survival.toml             # 3 seeds
# exp002 = áudio (WAV real ou sintético) · exp003 = vídeo
```

Cada run gera `runs/<exp>/run-<ts>-s<seed>/` (manifest, metrics.jsonl, snapshots por modalidade — PNG/waveform/WAV/frames —, population.json = a codificação, summary.json com baselines) e appenda 1 linha no `JOURNAL.md`.

## Modalidades — a MESMA física

`agua-life` é genérico (`World<A: Atom>`) e não sabe o que é pixel, amostra ou frame. Nova modalidade = implementar o trait `Atom` (render/mutate/rain) — zero mudança na física.

| Modalidade | Átomo | Alvo | Status |
|---|---|---|---|
| Imagem | `Gabor2D` | PNG ou sintético | ✅ validado |
| Áudio | `Gabor1D` (grão) | WAV PCM16 ou sintético | ✅ validado |
| Vídeo | `Gabor3D` (espaço-temporal) | sintético (frames PNG na bateria) | ✅ validado |
| Texto/embedding | átomo em espaço latente | prior importado | F4 (roadmap) |

## Crates

| Crate | Papel | Deps |
|---|---|---|
| `agua-core` | Field N-dim, Patch por runs, RNG determinístico, hotspot sampler, PSNR | zero |
| `agua-atoms` | trait `Atom` + Gabor 1D/2D/3D, fit de amplitude em forma fechada | core |
| `agua-life` | A física (genérica): comer/morrer/reproduzir/chover/ajustar | core, atoms |
| `agua-metrics` | Baselines obrigatórios (genéricos): random, matching pursuit; JPEG | core, atoms |
| `agua-telemetry` | RunDir, manifest, JSONL, PNG/WAV/waveform, JOURNAL.md | core |
| `agua-lab` | Runner: config TOML → dispatch de modalidade → run rastreável | todos |

## Documentos

- `docs/ARCHITECTURE.md` — design, invariantes, decisões
- `docs/ROADMAP.md` — fases F0→F6 com critérios de aceite
- `docs/EXPERIMENT_PROTOCOL.md` — regras de honestidade (multi-seed, baselines, falso-positivo)
- `CLAUDE.md` — guia pro agente que continuar este projeto
