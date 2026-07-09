# F1 — Bateria multimodal: resultados honestos (2026-07-09)

**192 runs**, 8 lotes, 3 modalidades, alvos sintéticos e reais, 3–5 seeds por variante.
Harness: `scripts/sweep.sh` (paralelo) + `scripts/collect.sh` (média ± spread por variante).
Todo run em `runs/<lote>-<variante>/run-<ts>-s<seed>/` + linha automática no `JOURNAL.md`.

---

## ⚠️ Falso-positivo morto (o mais importante)

No smoke do F0 reportei "**áudio e vídeo batem matching pursuit**". **Era artefato de um MP fraco** (24 candidatos por passo, orçamento ~30× menor que o da sobrevivência).

Ao dar ao MP o **mesmo orçamento de avaliações** (`mp_candidates = 0` → iso-budget), o resultado **inverte em todas as modalidades e todos os regimes**. A vitória não existia.

---

## Resultado central (negativo, robusto)

> **Para aproximação esparsa com dicionário paramétrico (Gabor), o matching pursuit guloso vence a sobrevivência em orçamento igual — sempre, em toda modalidade, todo regime de esparsidade e todo orçamento testado.**

| Lote | Modalidade / alvo | Melhor variante | PSNR (5 seeds) | vs MP iso-budget |
|---|---|---|---|---|
| L4 | imagem sintética, esparso | met=1.5 (69 átomos) | 27.02 ± 1.68 | **−4.58 dB** |
| L5 | imagem, orçamento baixo (t=100) | 5.9k evals | 22.91 ± 1.74 | **−6.46 dB** |
| L6 | áudio real (fala) | met=0.2 (1738 grãos) | 29.25 ± 3.08 | **−12.86 dB** |
| L6 | áudio real (trompete) | met=0.2 | 35.15 ± 0.90 | **−26.71 dB** |
| L8 | foto real 64×64 | rainHeavy | 23.33 ± 0.47 | **−5.09 dB** |

A sobrevivência **bate o random-atoms** com folga em todos os casos exceto um (`L6 speech-met08`, só 8 seres sobreviveram — colapso populacional, marcado ⚠ pelo collect). Ela é um otimizador legítimo; só não é competitivo com o guloso.

**Por que:** MP escolhe o **melhor de C candidatos** com resíduo atualizado a cada passo (descida coordenada gulosa). A sobrevivência gasta a maior parte dos evals propondo mutações locais que são **rejeitadas** (histerese de 2%), e só explora globalmente na chuva. Consistente: a variante que mais fecha o gap (`rainHeavy`, chuva contínua, −5.09) é justamente a que mais **se parece com MP** (mais amostragem global de candidatos novos).

## Lote 1 — metabolismo É o dial taxa-distorção (positivo limpo)

Varrendo `metabolism` (imagem sintética, cap folgado 3000, 3 seeds):

| metabolism | vivos | PSNR | bytes |
|---|---|---|---|
| 0.005 | 3000 (teto) | 41.73 | 84000 |
| 0.06 | 2759 | 39.78 | 77243 |
| 0.15 | 1325 | 35.42 | 37100 |
| 0.40 | 490 | 31.45 | 13729 |
| 1.00 | 129 | 28.82 | 3603 |

Monotônico e previsível ao longo de 3 ordens de grandeza. **A população de equilíbrio não é hiperparâmetro — emerge do custo de existir.** Isso o MP não faz: ele precisa que você diga N de antemão.

## Lote 2/8 — o que ajuda dentro da sobrevivência

- **Ajuste-de-vida (busca local) é essencial:** desligá-lo (`noTune`) piora o gap de −9.93 → −17.37 dB. Bate com a ponte survival↔inteligência (lote 4: aprender-durante-a-vida > evoluir do zero).
- **Mutação fraca é ruim** (`mut010`: −17.04). Mutação forte (0.4–0.6) é melhor por eval.
- **Chuva contínua > chuva rara** (`rainHeavy` −5.09 vs `base` −9.06): exploração global paga mais que refino local.
- **Mais ticks não fecha o gap** (`ticks12k`: gap piora pra −14.24) — a população cresce e MP escala melhor.

## Lote 7 — adaptação a sinal que muda (parcialmente positivo)

Alvo troca no tick 3000 (foto deslocada 4px + contraste), população **não reinicia**:

| orçamento pós-troca | PSNR sem adaptar | adaptado | MP a frio (iso-budget da fase B) |
|---|---|---|---|
| 37k evals | 19.63 | 21.14 | **28.97** |
| 205k evals | 19.63 | 21.43 | **31.80** |

A população **se re-harmoniza** (+1.5 a +1.8 dB) sem reiniciar — a água reescoa. Mas o MP resolvendo B **do zero** com o mesmo orçamento ainda ganha por ~8–10 dB. **A vantagem de adaptação existe mas não compensa a ineficiência de busca.**

## Pirâmide multi-escala: emerge (sinal real, não universal)

Fração da população em bandas finas (freq ≥ 0.12 ciclos/px), início vs fim do run:

| lote | delta médio | seeds positivas |
|---|---|---|
| L3-met080 | **+7.2 pp** | 4/5 |
| L5-t9000 | **+5.8 pp** | 4/5 |
| L8-rainHeavy | **+3.2 pp** | 5/5 |
| L8-base | +1.6 pp | 3/5 |

**16/20 runs** mostram grosso→fino. A seed 33 inverte consistentemente (nasce já com 19.8% fina). Veredito: **a pirâmide emerge da pressão de sobreviver, mas é tendência estatística, não lei.** (Ninguém programou "comece pelo grosso": átomos grossos explicam mais resíduo cedo, os finos só encontram comida depois que o grosso foi comido.)

## Custo em bytes: JPEG esmaga (esperado, não corrigido)

Foto 64×64: sobrevivência 23.3 dB com 5169 B (185 átomos × 7 f32). JPEG q30: **30.2 dB com 485 B**. Os parâmetros estão em f32 cru — 28 B/átomo. Quantização (F5) é pré-requisito pra qualquer alegação de codec.

---

## O que isto significa pra visão-água (honesto)

1. **A física funciona.** Determinística, converge, ecologia estável, bate random, descobre N sozinha, emerge pirâmide, adapta a mudança. Não há bug: `invariant_error < 1e-3` em todos os 192 runs.
2. **Ela não é competitiva onde o gradiente/greedy se aplica.** Em aproximação linear com átomo ajustável em forma fechada, MP é imbatível a igual orçamento. Insistir aqui é p-hacking.
3. **Portanto o diferencial da IARA-água tem que vir de onde MP NÃO se aplica** — exatamente o que a visão §4 previu e que ainda não testamos:
   - átomo **não-diferenciável / não-linear** (composição por oclusão/max, átomo = byte-model, átomo = pedaço de outro sinal);
   - **prior importado** → geração (F4), onde não há alvo pra fazer resíduo linear;
   - resíduo **conjunto multimodal** (sem produto interno único).
4. Os três achados positivos (metabolismo=dial RD, pirâmide emergente, adaptação online) são propriedades que **MP não tem** — mas custam caro. O caminho é preservá-las num regime onde o guloso não pode jogar.

## Recomendação de rota (revisa o ROADMAP)

- ~~F2 harmonia/CfC~~ e ~~F5 perf~~ **não são o gargalo**. Otimizar a busca não fecha 5–26 dB.
- **Pular direto pra F4 (prior importado / átomo não-linear)** — é o único lugar onde a tese pode ganhar de verdade. Primeiro experimento sugerido: composição por **max/oclusão** (MP não tem resíduo linear ⇒ não se aplica) no mesmo alvo, comparando com o melhor baseline aplicável.
- Antes disso, quantizar os params (F5-lite) pra que qualquer curva RD seja comparável a JPEG.
