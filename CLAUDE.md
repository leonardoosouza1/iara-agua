# CLAUDE.md — iara-água (solver universal por sobrevivência)

> Guia pro agente continuando este projeto. Fundação (F0) feita por Fable 5 em 2026-07-09.
> Leia nesta ordem: este arquivo → `docs/ROADMAP.md` (sua fase atual) → `docs/EXPERIMENT_PROTOCOL.md` (OBRIGATÓRIO antes de qualquer bateria) → `docs/ARCHITECTURE.md` (antes de mexer no código).

## O que é

Materialização da visão `../IARA_VISAO_SERES_COMO_AGUA.md`: população de seres-átomos (Gabors) sobrevive comendo o RESÍDUO de um sinal; a soma dos vivos é a reconstrução. Metabolismo = dial taxa-distorção. Tudo Rust headless, determinístico, rastreável.

## Regras do Leonardo (inegociáveis)

- **Honestidade brutal:** multi-seed, baselines sempre, matar o próprio falso-positivo, números exatos. Ver `docs/EXPERIMENT_PROTOCOL.md`.
- **Executar quando pedido** ("pode fazer, só não pode errar" = validar antes de dizer pronto). Só dizer "pronto" com tudo rodando.
- **Baterias contínuas:** terminou um lote → registra no JOURNAL → começa o próximo. Não esperar o Leonardo pra validar.
- **Rust puro; Python só orquestração/export offline.** Visualização = Bevy nativo (NUNCA browser).
- **Reusar, não recriar:** GA/niching de `../iara-swarm/crates/iara-evolution`; CfC do `../Universe/src/cfc.rs` (com teste de paridade); prior byte-model via export offline (`../bytebrain/`).

## Comandos

```bash
cargo test --release                                   # 15 testes, ~2s
cargo run --release --bin agua-lab -- demo             # smoke ~5s
cargo run --release --bin agua-lab -- run --config configs/exp001_image_survival.toml
# overrides: --seed N --ticks N --out dir
```

Sem GPU/env especial — CPU pura. (Se algum dia tocar torch/GPU: `HSA_OVERRIDE_GFX_VERSION=10.3.0`, pip sempre `--no-deps` — ver memória `reference_torch_rocm_setup`.)

## Estado (2026-07-09, fim do F0)

- 6 crates compilam, 15/15 testes. Demo: PSNR 15.8→26.9 dB em 800 ticks (400 seres, 4.8s), bate random 21.6, MP em 33.4 (gap 6.5 dB a fechar na F1), JPEG melhor por byte (params sem quantizar — F5).
- Ecologia real: nascimentos/mortes contínuos; bandas finas sub-povoadas ([284,90,16,6,4]) — primeira pergunta da F1.
- Cada run gera `runs/<exp>/run-*/` completo + linha no `JOURNAL.md`.

## Onde mexer pra quê

- Nova dinâmica de vida → `crates/agua-life/src/lib.rs` (`World::step`)
- Novo tipo de átomo → `crates/agua-atoms` (manter `fit_amp` em forma fechada)
- Novo baseline → `crates/agua-metrics`
- Novo experimento → copiar `configs/exp001_*.toml`, mudar `[experiment].name`
- NUNCA quebrar os invariantes: consistência target=recon+residual, determinismo por seed, amp em forma fechada (testes cobrem os três).
