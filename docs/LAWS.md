# As Leis da População — resultados consolidados

> 2026-07-09. Código: `agua-lab law2|law3|law4|law5|law6`. Cada lei tem predição feita **antes** da medida, e as falsificações estão registradas com o mesmo destaque dos acertos.

---

## LEI 1 — Metabolismo crítico

$$\text{met}_c \;=\; c\,\sigma^2$$

Um ser vive sse arranca do resíduo mais energia do que custa existir. Em ruído puro, o valor esperado do que ele arranca define o limiar.

| Predição | Medido | Veredito |
|---|---|---|
| c independe de σ | c = 12.94 ± 0.55, **CV 4.2%** ao longo de 20× em σ | ✅ |
| c independe da resolução | c plano (13–17) ao longo de **28× em N_pixels** | ✅ |
| c ∝ ln(orçamento de busca) | c = 1.62·ln(M) − 0.72, R² 0.87 | ✅ forma; constante ≠ 2 |
| c ∝ ln(N_pixels), Donoho–Johnstone | plano ⇒ **não é DJ** | ❌ falsificada |

**Por que não é Donoho–Johnstone:** o limiar universal (λ² = 2σ²·ln N) tem a mesma forma, mas o N dele é o dicionário inteiro, porque thresholding examina *todos* os coeficientes. Nossos seres fazem **busca local**. Consequência prática: **met_c não depende da resolução da imagem** — calibra uma vez, serve para qualquer tamanho.

---

## LEI 2 — c = 1 + a·dof, com a ≈ 2 (⟺ penalidade do **AIC**) ⭐

Medida limpa dentro da fase harmônica (ver Lei 6), congelando parâmetros um a um. A **primeira medida** (200 ticks) deu:

$$c \;=\; 2{,}007\,\cdot\,\text{dof} \;+\; 1{,}067 \qquad (R^2 = 0{,}9912)$$

⚠️ **Correção obrigatória (não celebrar cedo demais).** `c` também cresce com o orçamento de busca (Lei 1). Então testei se a *inclinação* é robusta ou se aquele 2.007 foi acidente do orçamento:

| ticks | inclinação a | intercepto b | R² |
|---|---|---|---|
| 50 | 1.370 | 1.598 | 0.969 |
| 200 | **2.023** | 1.142 | 0.993 |
| 800 | **2.243** | 0.889 | 0.992 |
| 3000 | **2.187** | 0.980 | 0.993 |

**O que sobrevive:**
- **Intercepto → 1.0**, em *qualquer* orçamento (1.14, 0.89, 0.98). É a âncora teórica: um átomo que não busca captura 1·σ², e mais tempo não muda isso. Medido isoladamente: 1.13. ✅
- **Inclinação satura em ~2.2** (2.02–2.24 na busca convergida). Cada grau de liberdade otimizado custa ≈2σ².

$$\boxed{\;\text{met}_c \;\approx\; \sigma^2\,\bigl(1 + 2\,\text{dof}\bigr)\;}$$

Dois por parâmetro é a penalidade do **AIC** (`−2 log L + 2k`) — ninguém programou isso; emerge de fome e morte. **Ressalva honesta: a medida (2.02–2.24) é consistente com 2, mas não o distingue de 2.2.** O valor exato 2.007 foi de um orçamento específico e não deve ser citado como "o" número.

---

## LEI 3 — Equipartição (preço sombra / KKT)

No equilíbrio, todo ser vivo explica ≈ met. Razões `E[explicado]/met` medidas: **1.04, 1.04, 1.00, 1.07** na faixa útil.

Se um ser explicasse muito mais que met, vizinhos nasceriam e dividiriam o resíduo; se menos, morreria. **Utilidade marginal = custo marginal** — a condição KKT emergindo sem otimizador.

---

## LEI 4 — Caça ao ruído = princípio de Morozov ⭐

Com met dado pela lei, a população para **exatamente no nível do ruído**, sem que ninguém lhe informe σ:

| σ | ‖r‖²/n final ÷ σ² |
|---|---|
| 0.10 | **0.99** |
| 0.15 | **1.01** |
| 0.20 | **1.03** |
| 0.30 | **1.02** |

É o princípio da discrepância de Morozov, implementado por metabolismo. É *por isso* que o denoise zero-shot funciona (`docs/LAW_METABOLISM.md`).

---

## LEI 5 — Universalidade (CLT), e onde ela quebra

O ser mede `⟨r,ψ̂⟩` — soma ponderada de centenas de pixels. Pelo CLT, isso é gaussiano **seja qual for a distribuição do ruído**.

| ruído | c (suporte grande) | vs Gauss |
|---|---|---|
| Gauss | 11.78 | — |
| Laplace (cauda pesada) | 11.63 | −1.2% |
| Uniforme (cauda leve) | 11.93 | +1.3% |
| Sal-e-pimenta (esparso) | 16.42 | **+39.5%** |

✅ **Universalidade confirmada** para ruído denso, qualquer cauda.
❌ **Minha predição de quebra estava errada na direção**: previ que sal-e-pimenta dispararia no suporte *pequeno*. Deu o oposto — dispara no grande (+39.5%) e **despenca** no pequeno (−54.7%), porque um átomo pequeno geralmente não pega **nenhum** impulso.

**O parâmetro de controle é `k = p·|suporte|`** (impulsos por átomo), não a cauda. Aumentando a densidade p, o desvio converge a ±5%:

| p | 0.01 | 0.03 | 0.10 | 0.30 | 0.60 | 1.00 |
|---|---|---|---|---|---|---|
| desvio (suporte pequeno) | −51.6% | −35.1% | +5.1% | +14.5% | −3.3% | −5.7% |

---

## LEI 6 — Transição de fase da harmonia ⭐

**H2 (falsificada):** eu previ que a sobrevivência ortogonalizaria a população. Ela faz o **oposto** — a população viva é **+96% a +167% mais coerente** que uma aleatória, com pares de μ ≈ 0.999 (quase-duplicatas). Eles nascem nos hotspots do resíduo, logo se aglomeram.

Essa colinearidade **infla os coeficientes** (como preditores colineares em regressão) e quebra a âncora: `c(dof=0)` dá 11.87 em vez de 1.

**Predação por nicho** (exclusão competitiva de Gause: quem tem μ > θ é devorado pelo mais forte, e o cadáver vira comida) restaura a lei — e a transição é **abrupta**:

| θ | 0.30 | 0.38 | 0.46 | **0.50** | 0.58 | 0.70 |
|---|---|---|---|---|---|---|
| c(dof=0) | 1.57 | 2.29 | 2.70 | 4.20 | 11.06 | 9.96 |
| fase | harmônica | harmônica | harmônica | ⟵ **θ\*** | colinear | colinear |

**θ\* ≈ 0.50** (60° entre átomos). Abaixo: população harmônica, âncora recuperada, lei válida. Acima: colinearidade, coeficientes inflados. Análogo às transições de fase de Donoho–Tanner em compressed sensing.

**A harmonia não emerge — precisa ser imposta por predação.** Foi essa a intuição do Leonardo ("seres que eliminam outros ou viram energia") e ela estava certa.

---

## LEI 8 — dof efetivo: hipótese **refutada**

A Lei 2 bate no Gabor2D (dof 6) mas o Gabor1D (dof 4) fica ~17% abaixo do previsto, sistematicamente. Hipótese: os parâmetros do 1D são degenerados (mover o centro ≈ mudar a fase), então o que conta é o posto do Jacobiano.

**Medida independente** (razão de participação dos autovalores da Gram de ∂ψ/∂p):

| átomo | dof nominal | dof_eff medido | (c−1)/2 esperado | erro |
|---|---|---|---|---|
| Gabor2D | 6 | 3.39 | 6.10 | **−44%** |
| Gabor1D | 4 | 2.71 | 3.25 | −17% |

❌ **Refutada.** E não troquei de métrica até alguma bater (isso seria p-hacking). O espectro da Gram mostra que o **1D tem posto 4 cheio** ([1.00, 0.51, 0.50, 0.015]) e é até *menos* desigual que o 2D ([1.00, 0.38, 0.37, 0.35, 0.099, 0.008]) — a degenerescência não existe.

**Explicação alternativa, essa sim sustentada:** o 1D só converge mais devagar. Com mais orçamento, `c_1D` sobe 7.50 → 7.87 → 8.22 → **8.34**, rumo a 9. O desvio era falta de convergência, não falta de graus de liberdade.

---

## LEI 7 — A predação paga? (teste iso-N)

Menos átomos com menos PSNR pode ser só andar na curva taxa-distorção. O teste honesto força o **mesmo N** (metabolismo calibrado por bisseção) e compara.

| N alvo | sem predação | com predação (θ=0.85) | Δ | vitórias |
|---|---|---|---|---|
| 120 | 30.75 dB | **31.05 dB** | **+0.29 dB** | 3/3 |

Ganho real de representação, não só de taxa — mas **modesto** (+0.29 dB). Em eficiência bruta o efeito é maior: PSNR/átomo sobe de 0.1117 para 0.2402, e o denoise mantém a qualidade com 53 átomos em vez de 90 (**−41%**).

*(N = 250 e 500 ainda rodando quando isto foi escrito; atualizar.)*

---

## Síntese

A população implementa, sem que ninguém programe:

1. **AIC** — cada grau de liberdade custa 2σ² (Lei 2, R² 0.99).
2. **Morozov** — para exatamente no nível de ruído (Lei 4).
3. **KKT** — utilidade marginal iguala custo marginal (Lei 3).
4. **Gause** — mas só se a predação for imposta; a harmonia não é espontânea (Lei 6).

O que a sobrevivência *não* faz é o ajuste contínuo — isso é do gradiente. A divisão de trabalho medida: **gradiente resolve o contínuo, seleção resolve o discreto** (quantos, quais, onde). É a mesma conclusão da ponte survival↔inteligência, agora com equação.

## Honestidade

- Cada lei foi medida em **1 seed** para varreduras exploratórias; as afirmações fortes (denoise, iso-N) usam 3–5 seeds. Marcado onde não.
- `c` na fase colinear é um número **inflado** e sem significado teórico — só a fase harmônica dá a lei.
- A constante γ = met_ótimo/met_c ≈ 2.7 continua **não explicada**.
- Nada aqui é comparado com BM3D/wavelet. **Zero alegação de estado-da-arte.**
