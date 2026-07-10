# A Lei do Garimpo — a penalidade é 2σ²·ln(M)

> 2026-07-10. `scripts/lei_garimpo.py`. Predições registradas antes de cada medida.
> Um erro de desenho meu está documentado abaixo, junto do conserto.

## A tese

A lei do metabolismo (`docs/LAWS.md`) mediu, dentro do solver de sinais:

$$\text{met}_c = c\,\sigma^2, \qquad c \approx 1 + 2\cdot\text{dof}, \qquad c \approx a\ln M_{\text{busca}}$$

Isso **não é uma propriedade de átomos Gabor**. É uma lei sobre qualquer sistema onde unidades competem por resíduo e você escolhe as vencedoras. E ela tem nome na literatura estatística: é o **RIC** (Risk Inflation Criterion, Foster & George 1994).

| Critério | Penalidade por parâmetro | O que ele assume |
|---|---|---|
| **AIC** | `2σ²` | você **não buscou** — o modelo estava dado |
| **BIC** | `σ²·ln n` | consistência assintótica em n amostras |
| **RIC** | `2σ²·ln p` | você **garimpou** entre p candidatos |

O AIC é o caso degenerado `M ≈ e`. A nossa lei mede `2σ²·ln(M_busca)` porque os seres **buscam** — cada grau de liberdade otimizado é uma dimensão a mais de garimpo.

Testei em três domínios fora do sistema dos seres.

---

## A) Garimpo puro — o Sharpe que o acaso produz

M estratégias com **zero alfa** (retornos puro ruído). Predição: `E[max Sharpe] ≈ σ_S·√(2 ln M)`.

| M estratégias testadas | E[max Sharpe] | √(2 ln M) | razão |
|---|---|---|---|
| 8 | 1.41 | 2.04 | 0.69 |
| 128 | 2.64 | 3.12 | 0.85 |
| 2048 | 3.48 | 3.91 | 0.89 |
| 8192 | **3.84** | 4.25 | 0.90 |

✅ Confirmada na forma; a razão converge devagar para 1 (mesma correção `O(ln ln M)` que vimos no EVT dentro do solver).

**Consequência prática, dura:** um backtest que testou **8192 variantes** produz Sharpe de **3.8 sem nenhum alfa**. Isso é o Deflated Sharpe Ratio de Bailey & López de Prado, e é a mesma equação da nossa lei do metabolismo.

→ **Implicação direta para o `saci-trader`**: o histórico daquele projeto tem um falso-positivo célebre (89% de acurácia que era vazamento de rótulo). O número de configurações varridas **precisa** entrar no limiar de aceitação. Um Sharpe de 2.0 depois de 500 tentativas é ruído.

---

## B) Seleção de variáveis (p ≫ n)

n = 200, k = 10 variáveis verdadeiras, σ = 1. F1 do suporte recuperado, 6 sementes.

⚠️ **Primeiro desenho meu: quebrado.** Fixei β ≈ 1 com colunas de norma 1 e σ = 1 — ou seja, coloquei o sinal **abaixo do limiar de detecção**. RIC e Lasso não selecionavam nada, e até o OMP com k-oráculo dava F1 = 0.30. Isso não é "o RIC perde", é um problema onde *nenhum* método pode funcionar. O parâmetro certo é β **relativo** ao limiar universal `σ√(2 ln p)`.

Refeito, com β = `mult · σ·√(2 ln p)`:

| p | mult | AIC (2σ²) | BIC | **RIC (2σ²·ln p)** | OMP (k-oráculo) |
|---|---|---|---|---|---|
| 1000 | 0.7 | 0.122 | 0.206 | **0.255** | 0.233 |
| 1000 | 1.0 | 0.349 | **0.473** | 0.463 | 0.583 |
| 1000 | 1.5 | 0.426 | 0.718 | **0.975** | 0.950 |
| 1000 | 2.5 | 0.434 | 0.726 | **1.000** | 1.000 |

✅ **Confirmada.** Acima do limiar de detecção, o RIC recupera o suporte perfeitamente e **iguala o OMP que recebeu o k verdadeiro**. O AIC satura em F1 ≈ 0.43 e **sobre-seleciona** (29 variáveis quando há 10) — exatamente o sintoma de quem não desconta o garimpo.

Abaixo do limiar (mult ≤ 1.0) o BIC ganha: quando o sinal é fraco demais para ser achado, a penalidade agressiva do RIC descarta os verdadeiros junto com o ruído. **Honesto: a lei tem um regime de validade.**

---

## C) Poda metabólica de rede ⭐

Um neurônio "sobrevive" se a perda que ele explica supera `2σ²·ln(W)`, onde W é a largura da camada — porque você garimpou W neurônios. Rede de random features, W = 200, σ = 0.5, 6 sementes.

| método | neurônios | MSE held-out |
|---|---|---|
| rede cheia | 200 | 1.1799 |
| **poda metabólica (2σ²·ln W)** | **12.7** | **0.6884** |
| poda por magnitude (mesmo k) | 12.7 | 1.1707 |
| AIC ingênuo (2σ², sem ln W) | 29.7 | 0.7146 |

✅ **Confirmada, 6/6 sementes.** A poda derivada da lei bate a poda por magnitude com folga no **mesmo orçamento de neurônios**, e bate o AIC ingênuo usando **menos da metade** dos neurônios. O `ln W` é o que faz a diferença.

Note também: a rede podada generaliza **melhor que a rede cheia** (0.69 vs 1.18). A poda não é só compressão — é regularização derivada de um princípio, sem threshold mágico.

---

## Síntese

$$\boxed{\;\text{penalidade por parâmetro} \;=\; 2\sigma^2 \ln M_{\text{alternativas}}\;}$$

Uma equação, quatro manifestações:

1. **Sinais** — o metabolismo crítico dos seres (`docs/LAWS.md`).
2. **Finanças** — o Sharpe que o garimpo produz (Deflated Sharpe).
3. **Estatística** — RIC vs AIC em seleção de variáveis com p ≫ n.
4. **Redes neurais** — poda por sobrevivência, sem threshold arbitrário.

O que muda entre elas é só **quem é M**: o orçamento de busca, o número de backtests, o número de candidatas, a largura da camada.

## Regimes de validade (o que não vale)

- **Sinal abaixo do limiar de detecção** (β < σ√(2 ln p)): a penalidade do garimpo mata os verdadeiros junto com o ruído. Ali o BIC é melhor.
- **Ruído esparso** (poucas amostras por unidade): o CLT não vale, e `c` se descola (medido: −55% com sal-e-pimenta e suporte pequeno).
- **A constante 2** é consistente com o AIC mas nossas medidas dão 2.0–2.3; não distinguimos.

## Próximo teste que eu faria

Aplicar (C) em pesos **reais** — as sub-redes do Qwen já mapeadas em `project_surviving_subnetworks` — e comparar a poda metabólica com magnitude/Wanda no mesmo orçamento de neurônios. σ ali é estimável pelo MAD dos resíduos de ativação. É o mesmo experimento, com um alvo que importa.

---

## Adendo — poda em pesos REAIS (ByteBrain 40M)

Apliquei (C) nos pesos treinados de verdade (`bytebrain/research/poda_metabolica.py`): MLP 640→2560→640, reconstruir a saída da camada em entradas com a estatística real dos embeddings. Poda = escolher os neurônios ocultos.

**Duas ressalvas que eu mesmo levantei antes de reportar:**
1. A metabólica empata o OMP **exatamente** — porque ela *é* OMP com regra de parada. "Metabólica ≈ OMP no mesmo k" é tautológico, não resultado.
2. Uma camada é **determinística** — não há ruído, logo `2σ²·ln W` não tem σ bem-definido. A lei escolhe orçamento por piso-de-ruído; sem ruído, o σ̂ vira um ponto de corte arbitrário (aqui, agressivo demais: 6–12 de 2560). **A auto-seleção de N pela lei só é significativa quando o alvo tem ruído genuíno** (denoise, regressão ruidosa). Para compressão exata de camada, o orçamento vem do erro aceitável.

**O resultado que sobra é a seleção a orçamento fixo** — qual método escolhe os melhores k neurônios (todos refazem o least-squares da saída; a diferença é só *quais* neurônios). MSE relativo da reconstrução, 4 blocos:

| orçamento | greedy-resíduo bate Wanda | bate magnitude |
|---|---|---|
| k=64 (agressiva) | **4/4** | **4/4** |
| k=256 | 3/4 | 4/4 |
| k=512 (leve) | 3/4 | 3/4 |

**Padrão limpo: quanto mais agressiva a poda, mais o critério de resíduo bate os heurísticos de importância** (Wanda, magnitude). No bloco 6, k=64: 0.556 vs 0.631 (Wanda). No regime leve todos convergem e o Wanda ganha de raspão por arredondamento.

**Honesto:** isto é reconstrução da SAÍDA DA CAMADA, não perplexidade do modelo inteiro — medir a perplexidade exigiria re-rodar o forward completo com as camadas podadas, o próximo passo real. E "greedy-resíduo > importância-heurística" é um resultado conhecido de sparse approximation; o que a lei acrescenta é o *critério de parada* (que precisa de ruído pra valer). A contribuição honesta: (a) confirma que o critério de resíduo vale em pesos reais, não só random features; (b) delimita quando a auto-seleção de N funciona.

---

## Adendo 2 — poda em PERPLEXIDADE REAL do LM: NEGATIVO (2026-07-10)

`bytebrain/research/poda_perplexidade.py`. Apliquei a poda metabólica nos
neurônios ocultos das MLPs do ByteBrain 40M e medi a bpb REAL em bytes held-out
(não reconstrução de camada), contra Wanda (SOTA one-shot) e magnitude.

| esparsidade | metabólica | Wanda | modelo cheio |
|---|---|---|---|
| 30% | 1.604 | **1.192** | 1.122 |
| 50% | 2.594 | **1.413** | |
| 70% | 3.808 | **2.046** | |

❌ **A poda metabólica PERDE do Wanda em todas as esparsidades**, decisivamente.
O refit dos pesos de saída não ajudou (e minha implementação de refit era no-op).

**Por que o resultado da reconstrução de camada (adendo 1, onde eu ganhava) NÃO
transferiu:**
1. Na perplexidade, TODAS as camadas são podadas juntas e os erros COMPÕEM pela
   profundidade.
2. Reconstruir o output L2 de uma camada isolada ≠ preservar o que as camadas
   SEGUINTES precisam. O greedy-resíduo escolhe um subconjunto pior pra função do
   modelo inteiro.
3. Meu desconto-de-redundância heurístico ativamente atrapalhou.

**Lição honesta e útil:** *a MSE de reconstrução de camada é um proxy ENGANOSO
para a qualidade de poda de LLM.* Um método pode ganhar reconstruindo uma camada e
perder na perplexidade real. A importância por-neurônio ciente-de-ativação do
Wanda é melhor pro problema real. A Lei do Garimpo continua válida como TEORIA
(RIC, finanças, seleção de variáveis) — mas sua aplicação à poda de LLM não bate o
SOTA. A aposta prática falhou; registro sem maquiar.
