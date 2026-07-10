# F4 — Geração por prior importado

> `agua-lab generate`. A promessa da visão §4: a sobrevivência deixa de
> RECONSTRUIR um alvo e passa a GERAR casando com um prior aprendido.

## O mecanismo

- **Prior importado:** codebook VQ (k-means) de patches 8×8 de uma imagem REAL.
  É a ponte offline do ecossistema no caso auto-contido — aprende de dado real.
- **Geração por sobrevivência:** tela nova, seres = ladrilhos (posição + índice
  do codebook). Comida = concordância com os vizinhos que se sobrepõem (costura
  coerente). Sem alvo. Sobrevive quem se harmoniza com o entorno.
- **O MP não se aplica:** não há sinal-alvo pra fazer resíduo. Só existe "isto é
  plausível localmente?".

## Resultado

| método | costura (↓) | novidade (↑) | copiou? |
|---|---|---|---|
| sobrevivência (k=32/64/128) | ~0.55 | ~0.04 | não |
| ladrilho aleatório | ~4.3 | ~0.04 | — |

- **Coerência: −87% de energia de costura** vs aleatório. A textura gerada é suave
  e estruturada; a aleatória é ruído.
- **Novidade > 0 em todos** (bem acima do limiar de cópia): gerou textura nova,
  não recortou o original.

## Veredito honesto

- ✅ **Provado:** a sobrevivência GERA — casa com um prior importado, sem alvo,
  onde o MP não pode jogar. Produz sinal novo e localmente plausível.
- ⚠️ **Limite conhecido:** é síntese de TEXTURA, não de CENA. Captura a estatística
  LOCAL (coerência de patch) mas não a estrutura GLOBAL — o gerado parece uma
  textura da montanha, não uma montanha nova. É o limite clássico de síntese
  baseada em patch sem um modelo global; para cena nova precisaria de um prior
  hierárquico (o próximo degrau, ligado ao latent-GPT do resonance_forest).
- **O que a peça acrescenta:** fecha o arco reconstrução→geração. A física é a
  mesma (nascer/harmonizar/morrer); muda só a comida (resíduo → plausibilidade).
