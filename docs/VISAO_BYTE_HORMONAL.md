# Visão fora da curva — o modelo-organismo byte-hormonal

> Ideia do Leonardo (2026-07-10): e se o modelo fosse TODO bytes, com órgãos,
> hormônios, e o "treino" INDUZISSE dopamina, estresse, ambiente? Não gradiente
> descendente numa loss escalar — mas um organismo que VIVE num ambiente de dados.
>
> Isto é fora da curva. Mas NÃO é mística: vários hormônios já têm fórmula EXATA
> desta sessão. Onde é hipótese, está marcado como hipótese.

## A troca de paradigma

Modelo de hoje: minimiza entropia cruzada — **uma escalar, um impulso, gradiente**.
O ambiente é um dataset estático; o modelo é uma função.

Modelo-organismo: **vive num ambiente de bytes**, movido por uma economia de
hormônios (cada um um sinal concreto). O "treino" não é uma loss — é sobreviver.

## Os hormônios, com a fórmula que a gente MEDIU

| hormônio | o que é | fórmula / status |
|---|---|---|
| **Metabolismo (saciedade)** | quantas peças o organismo sustenta | `met_c = σ²(1+2·dof)` — **PROVADO** = penalidade do AIC. Faminto → cresce capacidade; saciado → poda. |
| **Dopamina** | "minha predição melhorou mais que o esperado?" (RPE) | = a comida marginal = energia explicada a mais. **JÁ é o que calculamos** (a utilidade do ser). Gateia plasticidade. Universe `metabolism.rs` já tem dopamina=RPE. |
| **Cortisol (estresse)** | ameaça / surpresa sustentada (loss alta) | mapeia no ceticismo do garimpo `2σ²·ln M`: mais estresse = sarrafo mais alto, exploração estreita. **Aterrado** na Lei do Garimpo. |
| **Estabilidade (homeostase)** | fugir da beira do caos | = a susceptibilidade sináptica (Lei da Sinapse 2). O organismo se auto-organiza LONGE do crítico — **MEDIDO** (μ~0.03). |
| **Ambiente** | dados como pressão, não dataset | stream não-estacionário; o organismo prevê pra "comer". **A adaptação foi medida** (recupera 95% sob troca de sinal sem reiniciar). |

## Os órgãos = as línguas (a ponte poliglota, já construída)

Cada órgão é uma representação boa numa coisa; o roteador EMERGE da competição
(`docs/PONTE_LINGUAS.md`, medido):
- **byte** — robustez, estrutura fina, typo local (a "borda" da linguagem)
- **token-chunk** — padrões frequentes (a "textura" comum)
- **CfC** — ritmo lento no tempo (o "liso" temporal)

Cada um reivindica a região onde é insubstituível. O custo de ter todos é o
imposto de busca da Lei do Garimpo — só liga os que a diversidade paga.

## Por que isto pode valer (a aposta honesta)

A entropia cruzada é um bom impulso num mundo ESTÁTICO. A economia hormonal deve
ganhar num mundo NÃO-ESTACIONÁRIO (que muda) — porque:
- a **dopamina** foca a plasticidade onde há surpresa positiva (aprende o novo, não
  re-decora o velho);
- o **metabolismo** poda o que ficou obsoleto (esquece o que não paga mais);
- o **cortisol** sobe o sarrafo quando o mundo vira ameaça (não overfita o choque).

É onde o gradiente numa loss fixa sofre: mudança de distribuição, aprendizado
contínuo, ambiente vivo. **Essa é a hipótese testável.**

## O primeiro experimento (não poesia — número)

Um preditor de bytes minúsculo num STREAM que MUDA de distribuição no meio:
- **baseline:** SGD online numa entropia cruzada fixa.
- **hormonal:** taxa de aprendizado gateada por dopamina (RPE=surpresa) + poda
  metabólica de features obsoletas + sarrafo de cortisol sob choque.
- **métrica honesta:** quão RÁPIDO cada um se recupera da troca de distribuição
  (bpb ao longo do tempo). A vantagem prevista do hormonal é a adaptação.

Se o hormonal se adapta mais rápido à troca, a economia hormonal ganha ONDE deve
ganhar (mundo vivo). Se não, a ideia é bonita e morre — reportar assim.
**RODADO (docs/HORMONE_POC.md):** a economia completa naive FALHOU (instável); o hormônio isolado metabolismo=esquecer AJUDA a adaptar (choque −0.002 vs +0.022) mas custa 0.15 bpb no estável. Trade-off honesto, ganho marginal.
**TESTE 3 (scripts/hormone_poc2.py, multi-seed):** a ressalva "gateia ONDE, não a lr global" virou número. GATE LOCAL por-contexto ∈[0.1,1] CONSERTA o naive (instab 2.18→0.37, estável 5.34→3.32) e EMPATA o SGD (3.323 vs 3.328) — neuromodulação bem-feita é inócua, não vence o gradiente puro. Esquecer SELETIVO (só o já-aprendido) é o melhor trade-off: mesmo ganho de adaptação (choque +0.052 vs SGD +0.079) por ~1/3 do imposto (+0.056 vs +0.155). O design certo existe e é estável; ainda NÃO domina o SGD num stream só-metade-não-estacionário.

## A honestidade obrigatória (a lição do ecossistema)

O `resonance_forest` tinha ~10 manifestos grandiosos e código modesto — o abismo
manifesto↔código. Esta visão só vale se cada hormônio virar um NÚMERO. Metade já
virou (metabolismo, dopamina=RPE, susceptibilidade, adaptação). A outra metade
(a economia hormonal COMPLETA treina um byte-LM coerente?) é **hipótese não
testada** — e fica marcada como hipótese até um experimento honesto dizer sim ou
não. Sem manifesto sem número.
