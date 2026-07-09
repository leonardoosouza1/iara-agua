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

## LEI 2 — c = 1 + 2·dof ⟺ **Critério de Akaike** ⭐

Medida limpa dentro da fase harmônica (ver Lei 6), congelando parâmetros um a um:

$$c \;=\; 2{,}007\,\cdot\,\text{dof} \;+\; 1{,}067 \qquad (R^2 = 0{,}9912)$$

- **Intercepto 1.07**: um átomo que não busca nada captura exatamente 1·σ² (a âncora teórica; medido isoladamente: 1.13).
- **Inclinação 2.007**: cada grau de liberdade otimizado custa **2σ²**.

Dois por parâmetro é literalmente a penalidade do **AIC** (`−2 log L + 2k`). Ninguém programou isso. Portanto:

$$\boxed{\;\text{met}_c \;=\; \sigma^2\,\bigl(1 + 2\,\text{dof}\bigr)\;}$$

**O metabolismo crítico é o AIC, emergindo de fome e morte.** Seleção natural fazendo seleção de modelo.

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
