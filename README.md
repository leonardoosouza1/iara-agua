# iara-água — solver universal por sobrevivência

> **A IARA madura:** joga-se uma população de seres-átomos sobre um sinal (imagem, áudio, texto…); tentando **sobreviver** — comida = resíduo local — eles se auto-organizam **como água achando o caminho** e a soma dos vivos **É** a reconstrução/geração. Visão completa: `../IARA_VISAO_SERES_COMO_AGUA.md`.

## Quickstart

```bash
cargo test --release                                        # 15 testes
cargo run --release --bin agua-lab -- demo                  # smoke ~5s
cargo run --release --bin agua-lab -- run \
    --config configs/exp001_image_survival.toml             # 3 seeds
```

Cada run gera `runs/<exp>/run-<ts>-s<seed>/` (manifest, metrics.jsonl, snapshots PNG, population.json = a codificação, summary.json com baselines) e appenda 1 linha no `JOURNAL.md`.

## Crates

| Crate | Papel | Deps |
|---|---|---|
| `agua-core` | Field2D, Patch esparso, RNG determinístico, hotspot sampling, PSNR | zero |
| `agua-atoms` | Gabor2D (ser/alfinete), mutação, fit de amplitude em forma fechada | core |
| `agua-life` | A física: comer/morrer/reproduzir/chover/ajustar; World tick | core, atoms |
| `agua-metrics` | Baselines obrigatórios: random-atoms, matching pursuit, JPEG | core, atoms |
| `agua-telemetry` | RunDir, manifest, JSONL, snapshots PNG, JOURNAL.md | core |
| `agua-lab` | Runner: config TOML → run rastreável | todos |

## Documentos

- `docs/ARCHITECTURE.md` — design, invariantes, decisões
- `docs/ROADMAP.md` — fases F0→F6 com critérios de aceite
- `docs/EXPERIMENT_PROTOCOL.md` — regras de honestidade (multi-seed, baselines, falso-positivo)
- `CLAUDE.md` — guia pro agente que continuar este projeto
