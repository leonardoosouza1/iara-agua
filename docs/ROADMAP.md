# Roadmap — iara-água (fases com critério de aceite)

> Regra: cada fase só fecha com critério de aceite VALIDADO e registrado no JOURNAL.md.
> Modelo caro (Fable) fez F0; fases seguintes são executáveis por modelo mais barato
> seguindo `docs/EXPERIMENT_PROTOCOL.md`.

## ✅ F0 — Fundação multimodal (2026-07-09, Fable)
Workspace 6 crates, física de sobrevivência GENÉRICA (`World<A: Atom>`, Field N-dim, Patch por runs), átomos Gabor 1D/2D/3D, IO WAV/PNG/frames, baselines genéricos, telemetria, 23 testes, smoke das 3 modalidades.
**Aceite (batido):** `cargo test` 23/23; demos com MESMA física e MESMA config default:
- imagem 800 ticks: 17.5→25.9 dB (random 22.3, MP 35.6)
- áudio 800 ticks: 11.3→14.3 dB — **bateu MP (12.8)** e random (11.6); .wav ouvível exportado
- vídeo 400 ticks: 18.7→29.0 dB — **bateu MP (27.2)** e random (20.8)
(nota honesta: MP aqui é single-pass sem refit; sobrevivência tem relaxação contínua — em áudio/vídeo isso já paga. 1 seed; bateria confirma.)

## ✅ F1 — Bateria multimodal (192 runs, 2026-07-09) — VEREDITO NEGATIVO honesto
Resultados completos em `docs/RESULTS_F1.md`. Resumo: MP guloso iso-budget vence a sobrevivência em TODA modalidade/regime (−4.6 a −26.7 dB); o "áudio/vídeo batem MP" do F0 era falso-positivo (MP com 24 candidatos). Positivos reais: metabolismo = dial taxa-distorção (monotônico em 3 ordens de grandeza), pirâmide grosso→fino emerge (16/20 runs), adaptação online a sinal trocado (+1.8 dB sem reiniciar).
**Consequência:** F2 (harmonia/CfC) e F5 (perf) NÃO são o gargalo — otimizar busca não fecha 5–26 dB. **Pular pra F4** (átomo não-linear / prior importado), onde MP não se aplica.

### F1-bis (o que sobrou, opcional)
- Vídeo com frames PNG reais; quantização int8 dos params (pré-requisito pra comparar com JPEG: hoje 28 B/átomo f32, JPEG q30 faz 30.2 dB com 485 B).

## ~~F1 (plano original)~~
A pergunta da visão §5: seleção sem gradiente chega perto (ou bate) o matching pursuit? A pirâmide emerge?
- **Imagem (exp001):** sweep `metabolism` (dial taxa-distorção) × `cap` × `mutation_sigma` × `life_tuning on/off` × `tune_every` × ticks longos (20k+). Alvos: sintético + ≥2 PNGs reais 64×64 e 128×128 (fotos com textura fina).
- **Áudio (exp002):** WAVs reais (fala PT + música, PCM16) + sintético; bounds de grão longo (sigma_max 400). O smoke já bateu MP — validar multi-seed e contra MP com mais candidatos.
- **Vídeo (exp003):** carregar sequência de frames PNG (implementar no lab, ~20 linhas); medir se seres "surfam" o movimento (freq_t > 0 na população) vs viram estáticos.
- `life_tuning off` = evolução pura (só chuva+reprodução) — mede quanto o ajuste-de-vida paga (liga com lote 4 da ponte: aprender>evoluir).
- Pirâmide: plotar `freq_bands` ao longo do tempo — grosso primeiro, fino depois? (no demo imagem as bandas finas ficaram sub-povoadas [284,90,16,6,4] — investigar: metabolismo alto demais pro átomo fino? chuva enviesada?)
- **Aceite:** multi-seed (≥3) com veredito honesto POR MODALIDADE: gap médio pro MP em dB; curva taxa-distorção vs random/MP/JPEG; resposta sim/não da pirâmide. Registrar TUDO no JOURNAL.md.

## F2 — Órgãos temporais + harmonia explícita
- Portar `CfCCell` (teste de paridade numérica vs `Universe/src/cfc.rs`) → átomo com dinâmica temporal (sinal que muda no tempo; a ressonância freq↔decay já validada nas baterias 10/12 entra aqui).
- Termo de cooperação: a pergunta aberta da visão §7 — harmonia emerge sozinha (overlap negativo já pune via competição por resíduo) ou precisa de termo explícito? Medir tiling.
- Reuso: GA/niching de `iara-swarm/crates/iara-evolution`.
- **Aceite:** átomo temporal reconstrói sinal 2D+t melhor que átomo estático (multi-seed).

## ✅ F3 (antecipada pra F0) — Agnosticismo de modalidade
Feito na fundação: Field N-dim, Patch por runs, `World<A: Atom>` genérico, Gabor1D/2D/3D, WAV IO, alvos sintéticos das 3 modalidades, demos validados. `agua-life` tem ZERO conhecimento de modalidade (teste `survival_improves_psnr_audio` prova).
**Resta (vai pra F1):** alvos REAIS de áudio (WAV fala/música) e vídeo (frames PNG) nas baterias.

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
