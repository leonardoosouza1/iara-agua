# POC do modelo byte-hormonal — resultados honestos

> `scripts/hormone_poc.py`. Testa se a economia hormonal ganha ONDE deveria: num
> ambiente NÃO-ESTACIONÁRIO (stream de bytes que troca de idioma no meio, t=60k).

## Teste 1 — economia hormonal COMPLETA (naive): FALHOU

Hormonal = lr gateada por dopamina (surpresa) + cortisol (estresse) + decaimento
metabólico. Contra SGD online fixo.

| | antes da troca | após a troca |
|---|---|---|
| SGD fixo | 3.49 bpb | 3.51 |
| hormonal naive | 4.65 | 5.48 |

❌ **Pior em tudo, e INSTÁVEL** (bpb oscila 3.3→8.2). Amplificar a lr na surpresa
causa overshoot. A economia hormonal montada de forma ingênua é só um schedule de
lr ruim. Falhou.

## Teste 2 — o hormônio ISOLADO certo (metabolismo = esquecer): paga, modesto

lr estável nos dois; o hormonal só ESQUECE devagar os logits que não são
reforçados (o metabolismo aplicado à memória).

| | mundo A estável | choque de adaptação (pós−pré) |
|---|---|---|
| SGD (não esquece) | 3.49 bpb | **+0.022** |
| + metabolismo (esquece) | 3.64 | **−0.002** |

✅ **Esquecer AJUDA a adaptar** (choque menor). ⚠️ Mas **custa 0.15 bpb no mundo
estável** — trade-off honesto: esquecer paga quando o mundo muda, atrapalha quando
é estável. Ganho marginal, não dramático.

## Teste 3 — neuromodulação LOCAL (gateia ONDE, não a lr global): CONSERTA o naive

`scripts/hormone_poc2.py`. A ressalva do Teste 1 ("gateia ONDE aprender, não a lr
global") virou experimento. GATE por-contexto ∈ [0.1,1]: cada linha de logits
(contexto=byte anterior) tem plasticidade ALTA enquanto erra (`rl[ctx] > média
global`) e é CONGELADA quando já prevê bem — nunca passa da lr base (sem overshoot).
Multi-realização: 5 seeds sintéticos + 3 offsets reais (PT-literatura → instruct/código).

**Real (PT→código, o testbed que discrimina; sintético markov é ~incompressível 6.3 bpb):**

| arm | estável(pré) | choque(pós−pré) | instab(std do traço) |
|---|---|---|---|
| SGD fixo | 3.328 | **+0.079** | 0.369 |
| naive-global | 5.337 | +0.590 | **2.180** ❌ |
| **local-gate** | **3.323** | +0.087 | **0.373** |
| local+forget (blanket) | 3.483 | +0.052 | 0.349 |
| **local+selforget** | **3.384** | **+0.052** | 0.365 |

- ✅ **local-gate CONSERTA o naive:** instabilidade 2.180 → 0.373, estável 5.337 →
  3.323. A ressalva do doc estava CERTA — gatear ONDE (local) em vez de amplificar a
  lr global remove a instabilidade **por inteiro**. Mas **EMPATA com o SGD** (3.323 vs
  3.328; choque +0.087 vs +0.079) → neuromodulação bem-feita é INÓCUA aqui, não dá
  vantagem de adaptação sobre o gradiente puro (o gradiente já é maior onde erra mais).
- ✅ **local+selforget é o melhor trade-off achado:** esquece SÓ os contextos
  congelados (bem-aprendidos), escala `(1−gate)`. Mesmo ganho de adaptação do
  esquecer-tudo (choque **+0.052**, bate SGD +0.079) por **~1/3 do imposto no estável**
  (+0.056 vs +0.155 do blanket). Reduz o custo do metabolismo do Teste 2 em ~3×.
- ⚠️ Ainda **não DOMINA o SGD:** paga o estável pra comprar adaptação (num stream só
  metade não-estacionário, o SGD ganha na média). É um trade-off MELHOR, não uma vitória
  limpa. Sintético markov não discrimina (loss satura ~6.3 → surpresa fica baixa → nem o
  naive explode lá; explode no texto real, onde a loss varia 3↔8).

## Veredito da visão

- **Aterrado:** vários hormônios têm fórmula desta sessão (metabolismo=AIC,
  dopamina=RPE, estresse=susceptibilidade sináptica).
- **A economia COMPLETA naive falha** (amplificar lr por surpresa = instável), MAS o
  **design cuidadoso que o doc pediu EXISTE e funciona:** (a) o gate LOCAL remove a
  instabilidade e empata o SGD; (b) o esquecer SELETIVO corta o imposto do metabolismo
  em ~3× mantendo o ganho de adaptação. A neuromodulação real **gateia ONDE**, e isso
  foi validado.
- **Um benefício real e pequeno:** esquecer (metabolismo) ajuda a adaptar em ambiente
  não-estacionário; seletivo é ~3× mais barato que blanket. Custa no estável — trade-off,
  não almoço grátis. Não dominou o SGD num stream só-metade-não-estacionário.
- **Sem manifesto sem número.** O número diz: promissor em pedaços, modesto onde
  funciona, e o conserto do naive (gate local + esquecer seletivo) está MEDIDO. A visão
  vive como hipótese com evidência agora mais forte (o design certo é estável e o
  trade-off melhorou), não como vitória sobre o baseline.
