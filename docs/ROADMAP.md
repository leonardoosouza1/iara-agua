# Roadmap — iara-água (fases com critério de aceite)

> Regra: cada fase só fecha com critério de aceite VALIDADO e registrado no JOURNAL.md.
> Modelo caro (Fable) fez F0; fases seguintes são executáveis por modelo mais barato
> seguindo `docs/EXPERIMENT_PROTOCOL.md`.

## ✅ F0 — Fundação (2026-07-09, Fable)
Workspace 6 crates, física de sobrevivência v0, baselines, telemetria, 15 testes, smoke validado.
**Aceite (batido):** `cargo test` verde; demo 800 ticks: PSNR 15.8→26.9 dB, bate random (+5.4 dB), ecologia com turnover (781 nasc/381 mortes), artefatos completos, JOURNAL appendado.

## F1 — Ciência: imagem-por-sobrevivência (bateria ~50 runs) ← PRÓXIMA
A pergunta da visão §5: seleção sem gradiente chega perto do matching pursuit? A pirâmide emerge?
- Sweep: `metabolism` (dial taxa-distorção) × `cap` × `mutation_sigma` × `life_tuning on/off` × `tune_every` × ticks longos (20k+).
- `life_tuning off` = evolução pura (só chuva+reprodução) — mede quanto o ajuste-de-vida paga (liga com lote 4 da ponte: aprender>evoluir).
- Alvos: sintético + ≥2 PNGs reais 64×64 e 128×128 (fotos com textura fina).
- Pirâmide: plotar `freq_bands` ao longo do tempo — grosso primeiro, fino depois? (no demo v0 as bandas finas ficaram sub-povoadas [284,90,16,6,4] — investigar: metabolismo alto demais pro átomo fino? chuva enviesada?)
- **Aceite:** multi-seed (≥3) com veredito honesto: gap médio pro MP em dB; curva taxa-distorção sobrevivência vs random vs MP vs JPEG; resposta sim/não da pirâmide. Registrar TUDO no JOURNAL.md.

## F2 — Órgãos temporais + harmonia explícita
- Portar `CfCCell` (teste de paridade numérica vs `Universe/src/cfc.rs`) → átomo com dinâmica temporal (sinal que muda no tempo; a ressonância freq↔decay já validada nas baterias 10/12 entra aqui).
- Termo de cooperação: a pergunta aberta da visão §7 — harmonia emerge sozinha (overlap negativo já pune via competição por resíduo) ou precisa de termo explícito? Medir tiling.
- Reuso: GA/niching de `iara-swarm/crates/iara-evolution`.
- **Aceite:** átomo temporal reconstrói sinal 2D+t melhor que átomo estático (multi-seed).

## F3 — Agnosticismo de modalidade (áudio 1D)
- `Field1D` + grão de Gabor 1D; `agua-life` INTACTO (é o teste: a física não pode mudar).
- Alvo: trechos de áudio reais (fala PT + música), PSNR/SNR vs MP 1D.
- **Aceite:** mesma física resolve 2 modalidades; diff de `agua-life` = zero linhas.

## F4 — Prior importado → GERAÇÃO
- O diferencial não-diferenciável: átomo/comida definidos por um MODELO (byte-model/latent), onde gradiente não serve.
- Ponte offline Python→Rust (mesmo padrão do `world_knowledge.json`): exportar prior (ex.: distribuição do latent-GPT do resonance_forest, ou ByteBrain 40M) → comida = casar com o prior → população AMOSTRA sinais novos.
- **Aceite:** geração de amostras novas (não reconstrução) qualitativamente plausíveis + métrica honesta (bpb/FID-like simples) vs amostrar o prior direto.

## F5 — Escala & performance
- rayon, profiling, 256×256+, populações 10k+; quantização int8/int4 dos params (bytes honestos na curva RD — hoje 28 B/átomo f32 é generoso com o JPEG).
- GPU (HIP/WGSL) SÓ se o profiling provar gargalo (regra: pesquisar antes de otimizar).
- **Aceite:** 256×256 em <60s/1k ticks OU justificativa medida de por que não.

## F6 — Visibilidade viva (com aval do Leonardo)
- Observatório NATIVO Bevy no Universe (nunca browser): população escoando sobre o sinal, seres coloridos por banda de freq, clique inspeciona átomo.
- Alternativa barata antes: `agua-lab report` gerando painel PNG/MD por run.
- **Aceite:** Leonardo vê a água escoar e aprova.

## Perguntas abertas da visão (§7) — onde respondem
1. Seleção ≈ matching pursuit? → F1
2. Quantos seres pra que dB (Pareto/enxugar)? → F1 (dial metabolismo)
3. Harmonia emerge ou precisa termo explícito? → F2
4. Vale vs gradiente? O caso onde gradiente NÃO serve → F4
5. Multimodal com resíduo conjunto → pós-F4 (F3 prova o caminho)
