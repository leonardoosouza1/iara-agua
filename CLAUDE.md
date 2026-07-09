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
cargo test --release                                   # 23 testes, ~3s
cargo run --release --bin agua-lab -- demo             # imagem ~4s
cargo run --release --bin agua-lab -- demo-audio       # áudio (exporta .wav ouvível)
cargo run --release --bin agua-lab -- demo-video       # vídeo ~16s
cargo run --release --bin agua-lab -- run --config configs/exp001_image_survival.toml
# exp002 = áudio · exp003 = vídeo · overrides: --seed N --ticks N --out dir
```

Sem GPU/env especial — CPU pura. (Se algum dia tocar torch/GPU: `HSA_OVERRIDE_GFX_VERSION=10.3.0`, pip sempre `--no-deps` — ver memória `reference_torch_rocm_setup`.)

## Estado (2026-07-09, fim do F0 — fundação MULTIMODAL)

- 6 crates, 23/23 testes. Física genérica `World<A: Atom>` (Field N-dim, Patch por runs) — `agua-life` não sabe o que é modalidade.
- Smokes (mesma física/config, 1 seed): imagem 25.9 dB (random 22.3, MP 35.6) · áudio 14.3 dB **bateu MP 12.8** · vídeo 29.0 dB **bateu MP 27.2**. Nota: MP é single-pass sem refit; a relaxação contínua da sobrevivência paga em áudio/vídeo — bateria F1 confirma multi-seed.
- Bandas finas sub-povoadas na imagem ([284,90,16,6,4]) — primeira pergunta da F1.
- Cada run gera `runs/<exp>/run-*/` completo (imagem: PNGs; áudio: waveform+WAV; vídeo: frames) + linha no `JOURNAL.md`.

## Onde mexer pra quê

- Nova dinâmica de vida → `crates/agua-life/src/lib.rs` (`World::step`) — GENÉRICA, não pode conhecer modalidade
- Nova MODALIDADE → implementar trait `Atom` em `crates/agua-atoms` (render_unit/mutate/rain_at/freq) + dispatch/IO em `agua-lab` — física intocada
- Novo baseline → `crates/agua-metrics` (genérico sobre `A: Atom`)
- Novo experimento → copiar `configs/exp00N_*.toml`, mudar `[experiment].name` e `[world].modality`
- NUNCA quebrar os invariantes: consistência target=recon+residual, determinismo por seed, amp em forma fechada, vida agnóstica (testes cobrem os quatro).
