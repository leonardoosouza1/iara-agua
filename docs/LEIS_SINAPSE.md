# As leis da sinapse — acoplamento, alcance, criticalidade

> `agua-lab synapse`. Pergunta do Leonardo: existe lei de sinapse? disparo? raio
> de distância? como um neurônio que dispara atinge outro? por que o outro dispara?

## A sinapse já estava lá, implícita

A reconstrução ótima resolve as equações normais **G·a = b**, onde
**G_ij = ⟨ψ_i, ψ_j⟩ é a matriz de Gram — e G É a matriz de pesos sinápticos.**
Quando o átomo j "dispara" (muda amplitude em Δ), o resíduo que i sente muda em
−Δ·G_ij, e a amplitude ótima de i se desloca em −Δ·G_ij/G_ii. **Essa é a
transmissão sináptica.** Como um neurônio atinge outro: pelo overlap dos suportes.
Por que o outro dispara: porque o resíduo mudou e sua amplitude ótima mudou junto.

## Lei 1 — o acoplamento É o campo receptivo (chapéu-mexicano)

G_ij entre dois átomos idênticos separados por distância d, ao longo da orientação:

| d (px) | G_ij | sinal |
|---|---|---|
| 0 | +1.00 | excita |
| 2 | +0.06 | excita |
| 3–6 | −0.58 a −0.13 | **inibe** |
| 8–10 | +0.51 a +0.11 | excita |
| 12 | −0.22 | inibe |

**O acoplamento ALTERNA excitatório/inibitório com a distância**, sob um envelope
gaussiano de comprimento λ ≈ 9.9 px ≈ σ√2. Isso é o **chapéu-mexicano / inibição
lateral** do córtex visual (center-surround). A sinapse entre dois átomos = a
autocorrelação deles = o próprio Gabor, que é um campo receptivo de V1. **A forma
do átomo determina a forma da sinapse.**

## Lei 2 — a susceptibilidade diverge no colinear (criticalidade)

Perturbar um átomo (disparo) e medir ‖Δamp‖/‖disparo‖ resolvendo G·a = b, em
populações de coerência controlada (raio R menor = átomos empilhados = μ maior):

| μ (coerência) | cond(G) | susceptibilidade |
|---|---|---|
| 0.147 | 205 | 74 |
| 0.253 | 273 | 782 |
| 0.381 | 295 | 7 932 |
| 0.503 | 241 | 4 908 |
| **0.600** | **2101** | **120 289** |

**A resposta a um disparo DISPARA ~1600× conforme μ cresce.** Átomos colineares ⇒
G quase-singular ⇒ um disparo minúsculo causa oscilações de amplitude gigantes que
se cancelam. É o análogo da criticalidade / avalanche neural: perto do ponto
crítico, a rede responde a qualquer cutucão de forma explosiva.

⚠️ **Honesto:** métricas ingênuas de avalanche (contar quantos átomos se moveram,
alcance espacial) NÃO mostram isso — são dominadas pela geometria. A métrica certa
é a **susceptibilidade** (condição da matriz sináptica), e ela é limpíssima.

## Uma história só (a conexão com a Lei 6)

Essa é **a mesma colinearidade que inflou o c na Lei 6**: alta coerência → G
mal-condicionada → amplitudes infladas (c=12 em vez de 1) → resposta a disparo
explosiva. Metabolismo, harmonia e sinapse são a mesma física vista de ângulos
diferentes.

**E a sobrevivência foge da criticalidade:** a dinâmica mantém μ ~ 0.03 (regime
harmônico, subcrítico), longe da borda instável — o **oposto** de self-organized
criticality. O seletor por sobrevivência é um ESTABILIZADOR: espalha os átomos pra
que um disparo não derrube a rede. Um dispara, poucos respondem, tudo estável.
