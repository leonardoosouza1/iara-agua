# A Lei do Metabolismo Crítico

> Derivada, testada e parcialmente falsificada em 2026-07-09. Código: `crates/agua-lab/src/law2.rs` (`agua-lab law2`), aplicação: `agua-lab denoise --auto`.

## Enunciado

Um ser que ajusta sua geometria contra ruído branco N(0,σ²) captura energia espúria. Sobreviver exige capturar mais que o custo de existir. Logo existe um **metabolismo crítico** abaixo do qual até o ruído paga a conta:

$$\boxed{\;\text{met}_c \;=\; c\,\sigma^2, \qquad c \;\approx\; a\,\ln M \;+\; b\;}$$

onde **M é o orçamento de busca do ser** (nº de geometrias que ele testa) — não o tamanho da imagem.

## Derivação

Para um átomo unitário fixo ψ̂ sobre ruído puro, ⟨r,ψ̂⟩ ~ N(0,σ²), e a energia explicada é E = ⟨r,ψ̂⟩². Portanto **E[E] = σ²** — um átomo qualquer captura exatamente uma unidade de variância.

Mas o ser **escolhe** a melhor entre M geometrias contra *aquela realização* do ruído. Então E/σ² é o máximo de M variáveis χ²₁, e a teoria de valores extremos dá E[max] ≈ 2 ln M. Daí `met_c = 2σ²·ln M`.

## O que foi medido

| Predição | Resultado | Veredito |
|---|---|---|
| Âncora: átomo **sem busca** captura 1·σ² | 1.13 (2D), 0.91 (1D) | ✅ confirmada |
| **Q2**: c independe de σ | c = 12.94 ± 0.55 (**CV 4.2%**) ao longo de 20× em σ | ✅ confirmada |
| **Q1**: c linear em ln(M), inclinação 2 | c = **1.62**·ln(M) − 0.72, R² = 0.87 | ⚠️ forma certa, inclinação 1.62 |
| **P3**: c cresce com graus de liberdade | 13.69 (6 geom) vs 7.25 (4 geom) | ✅ direção certa, razão 1.89 ≠ 1.50 |
| **Q3**: extinção abrupta em met_c | 600→334, 238→66, 131→31 | ⚠️ transição **suave**, não degrau |
| **Q4**: c ∝ ln(N_pixels) (Donoho–Johnstone) | c plano (13–17) ao longo de **28× em N** | ❌ **falsificada** |

**Por que a inclinação é 1.62 e não 2:** os passos de Adam são correlacionados, então o número *efetivo* de geometrias independentes é menor que `ticks × competidores`. A forma logarítmica sobrevive; a constante não.

**Por que Q3 é suave:** `met_c` é a **média** do explicado. A cauda da distribuição sobrevive acima dela. A extinção é probabilística: sobrevive a fração P(E > met).

**Por que Q4 falhou — e por que isso é bom:** o limiar universal de Donoho–Johnstone (λ = σ√(2 ln N)) tem exatamente a mesma forma em energia (λ² = 2σ² ln N), mas o `N` dele é o tamanho do dicionário, porque o thresholding examina **todos** os coeficientes. Nossos seres fazem **busca local** (gradiente + chuva em hotspots do resíduo). Consequência medida: **met_c independe da resolução** (c constante de 576 a 16.384 pixels). Calibra-se uma vez, vale pra qualquer tamanho — o que o limiar universal não permite.

⚠️ **Confound que quase me enganou:** na primeira medição fixei a população em 600 seres para todos os campos. Num campo de 24×24 isso significa 600 seres devorando 576 pixels de ruído; o "explicado" despenca por competição, não por escala. Deu inclinação 5.38 e R² 0.98 — uma lei falsa, bonita e errada. Com **densidade de seres constante**, a tendência some. Sempre fixar a densidade ao variar o tamanho do domínio.

⚠️ **Erro da v1 (`law.rs`, mantido no repo como registro):** medi "população viva após T ticks" como proxy do limiar. Mas o tempo até a morte é `energia_inicial/(met − explicado)` — com met pequeno o ser demora a morrer e *parece* sobreviver. Aquilo media tempo de morte, não limiar. Deu expoente 0.75 em vez de 2.

## Aplicação: denoise sem hiperparâmetro

σ é estimável do próprio sinal ruidoso (MAD do detalhe diagonal, Donoho–Johnstone 1994): σ̂ = MAD/0.6745. Erro medido: **~3%**.

Então `met = k·σ̂²` fecha o laço. A população decide sozinha **quantos** átomos usar — que é exatamente o que o matching pursuit **não** sabe fazer (ele exige que você diga N).

Calibrado em σ ∈ {0.05, 0.10, 0.20} → k ≈ 40. Testado em σ ∈ {0.15, 0.30} (**held-out**), 3 alvos × 3 seeds:

| alvo | σ | ruidoso | gauss@oráculo | **sobrevivência** | MP@N-oráculo | vida − MP |
|---|---|---|---|---|---|---|
| photo_smooth | 0.15 | 16.44 | 28.17 | 27.46 | 26.08 | **+1.38** |
| photo_smooth | 0.30 | 10.42 | 25.75 | 25.23 | 22.79 | **+2.45** |
| photo_texture | 0.15 | 16.44 | 24.44 | 23.88 | 22.62 | **+1.26** |
| photo_texture | 0.30 | 10.42 | 21.63 | 20.67 | 19.67 | **+1.00** |
| synthetic | 0.15 | 16.44 | 24.17 | **28.41** | 24.08 | **+4.33** |
| synthetic | 0.30 | 10.42 | 21.90 | **26.28** | 21.24 | **+5.04** |

- **18/18 runs** batem o MP mesmo quando o MP escolhe seu N **vendo a imagem limpa**. Ganho médio **+2.58 dB**.
- Contra o filtro gaussiano com σ_blur escolhido por oráculo: **6/18**. Ganha +4.3 dB onde o átomo casa com o sinal (oscilação); perde ~0.7 dB em foto natural, dominada por suavidade.
- O N escolhido pela lei fica a **36% (mediana)** do N-oráculo do MP, sem nunca ver o limpo. Num caso: 291 átomos contra 293 do oráculo.

**Honestidade obrigatória:** `gauss@oráculo` é um teto clássico *fraco* e ainda assim ganha nas fotos. **Não comparamos com BM3D, NL-means ou wavelet-thresholding** — não há nenhuma alegação de estado-da-arte aqui. O resultado defensável é estreito e específico: *a seleção automática de N por metabolismo bate o matching pursuit com N-oráculo, na mesma família de átomos.*

## Constante em aberto

O metabolismo **ótimo para PSNR** (k ≈ 40) fica ~2.7× acima do **crítico** medido (c ≈ 14.6 no orçamento do denoise). Faz sentido — acima do crítico você também poda átomos que pegam sinal fraco, e ganha por não ajustar ruído; é o mesmo motivo pelo qual o limiar universal de D–J é conservador e o SURE difere dele. Se `γ = met_ótimo/met_c` for constante entre sinais, vira a próxima lei. **Não testado.**

## Ligações com o resto do ecossistema

- **MDL/BIC**: "sobreviver = explicar mais energia que o custo dos próprios graus de liberdade" é seleção de modelo por comprimento de descrição. O metabolismo é o preço em bits.
- **Donoho–Johnstone**: mesma *forma* (2σ² ln M), `M` diferente (busca local vs dicionário inteiro). É a razão da independência de resolução.
- **[[project_bridge_survival_intelligence]]**: a sobrevivência não faz o ajuste contínuo (gradiente faz), mas decide o discreto — aqui, **quantos átomos**. É a mesma divisão de trabalho da ponte: gradiente/importação faz o contínuo, seleção faz a estrutura.
