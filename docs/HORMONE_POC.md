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

## Veredito da visão

- **Aterrado:** vários hormônios têm fórmula desta sessão (metabolismo=AIC,
  dopamina=RPE, estresse=susceptibilidade sináptica).
- **Um benefício real e pequeno:** o metabolismo (esquecer) ajuda em ambiente
  não-estacionário, com custo no estável.
- **A economia COMPLETA naive falha:** amplificar plasticidade por surpresa é
  instável. Precisa de design cuidadoso (neuromodulação real gateia ONDE aprender,
  não amplifica a lr global).
- **Sem manifesto sem número.** O número diz: promissor em pedaços, modesto onde
  funciona, instável onde foi ingênuo. A visão vive como hipótese com uma primeira
  evidência (fraca), não como resultado.
