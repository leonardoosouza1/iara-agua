# A ponte entre línguas — o roteador que emerge (tese IARA-Router)

> `agua-lab polyglot`. A pergunta do Leonardo: e se cada representação é boa numa
> coisa? bytes em X, tokens em Y, CfC em Z... um tradutor que trabalha em várias.

## A ideia, tornada concreta

O seletor por sobrevivência JÁ é um roteador — átomos competem pelo resíduo. A
generalização: os candidatos vêm de FAMÍLIAS diferentes (línguas), cada região do
sinal é reivindicada por quem a explica mais barato. **O tradutor não é um modelo
que converte — é uma arena onde as línguas competem por região, e o roteamento
emerge.** Nenhum roteador é treinado.

3 línguas complementares: `gabor` (textura periódica), `blob` (regiões lisas),
`borda` (descontinuidade — o que o Gabor aproxima mal: uma soma de senos RINGA
numa borda dura).

## Resultado — o valor é LOCAL, não global

Alvo heterogêneo (3 faixas: lisa, textura, degrau). PSNR **por região**:

| região | poliglota | só gabor | Δ |
|---|---|---|---|
| degrau (borda dura) | **20.26** | 14.73 | **+5.53 dB** |
| liso (gradiente) | 22.49 | 20.92 | +1.57 dB |
| textura (listras) | 16.67 | 16.27 | +0.40 dB |

**O PSNR GLOBAL quase não mexe** (poliglota 19.60 vs gabor 19.43) — porque a média
dilui a vitória, que está concentrada na região onde uma língua especializada
domina. Olhar só o global esconde o ganho de +5.5 dB na borda.

## O roteamento emerge

A língua dominante por faixa, sem treinar roteador:
- degrau → **borda** ✓
- textura → **gabor** ✓
- liso → gabor domina a contagem, mas o PSNR mostra que o blob ajuda (+1.57)

## O regime de validade (a Lei do Garimpo cobra)

Alvo HOMOGÊNEO (só textura): poliglota **16.90 vs gabor 19.56** — o poliglota
PERDE. Misturar 3 línguas quando só uma serve desperdiça 2/3 do orçamento de busca
em famílias inúteis. É a Lei do Garimpo medida: você paga `2σ²·ln(M·C)` de imposto
de busca por ter M línguas. Só vale se a diversidade das regiões comprar mais que
o imposto.

## A lei que sai disso

> **Um dicionário poliglota ganha onde o sinal é HETEROGÊNEO e cada língua tem uma
> região onde é insubstituível. Perde onde o sinal é homogêneo (imposto de busca à
> toa). O roteamento por região emerge da competição por resíduo — não precisa de
> roteador treinado.**

Contraste honesto com o M74 do resonance_forest (generalista > swarm-de-especialistas
quando os domínios COMPARTILHAM estrutura): aqui é o oposto complementar — as línguas
são COMPLEMENTARES (não sobrepostas), então especializar POR REGIÃO paga. A regra
unificada: especialização paga quando as unidades cobrem coisas DISJUNTAS; transferência
(generalista) paga quando cobrem coisas COMPARTILHADAS.

## Ligação com o ecossistema (a resposta à pergunta)

bytes/tokens/CfC/Gabor como línguas: uma arena de sobrevivência onde cada uma
reivindica o que faz melhor, POR REGIÃO do sinal/sequência, e o roteador emerge.
- bytes: robustez, universalidade, estrutura fina (typo local) — a "borda" da linguagem
- tokens: compressão de padrões frequentes — a "textura" comum
- CfC: dinâmica temporal lenta — o "liso" no tempo
O custo de ter todas é o imposto de busca; a regra acima diz quando ligar quais.
