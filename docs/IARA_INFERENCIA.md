# IARA — Inferência como SELEÇÃO (não geração)

> A tese-mãe: a IARA não pinta o sinal pixel a pixel. Ela **seleciona quais seres-átomo sobrevivem**,
> e a soma dos vivos **é** a resposta. A "inferência" dela é paralela, determinística (amplitude em
> forma fechada, sem alucinar — um ser só existe se explica algo real) e **agnóstica ao domínio**.
> Tudo aqui é medido com honestidade brutal (multi-seed, baselines, oráculo como teto).

Este doc mapeia os **modos de inferência** já implementados e validados, e a **fronteira** que os separa.

## Os modos (todos = seleção sobre uma restrição diferente)

| modo | restrição | comida do ser | resultado honesto | onde |
|---|---|---|---|---|
| **Reconstruir** | um alvo | come o resíduo | 24–31 dB, **população emerge** da lei do metabolismo | `run` / `demo*` |
| **Gerar** | um prior (codebook) | concordar com o vizinho (costura) | **87% mais coerente que aleatório, novidade>0** (gerou, não copiou) | `generate` |
| **Compor não-linear** | oclusão (não-diferenciável) | melhora da composição | vence o linear **6/9**, e o MP **nem se define** ali | `occlusion` |
| **Completar (inpaint)** | pixels visíveis | casar o visível + costurar | seleção pura **não gera** o não-visto; **com prior, gera** | `inpaint --prior` |
| **Rotear (poliglota)** | qual expert por região | verificador escolhe | **≥ difusão nos 6 cenários, fecha o oráculo** | `inpaint --poliglota` |

## A fronteira (o que a seleção pode e não pode)

**Reconstruir o observado, a seleção vence.** Come o resíduo, harmoniza (Gauss-Seidel), a população de
equilíbrio emerge da lei do metabolismo (AIC): `met_c = σ²(1+2·dof)`. A amplitude é sempre forma-fechada
`amp = ⟨r,ψ⟩/‖ψ‖²` — sem isso EXPLODE (foi um bug real). Sobreviver só esculpe a geometria.

**Gerar o não-observado do NADA, a seleção NÃO faz.** No inpaint por seleção pura o buraco vira lixo
(PSNR negativo): sem sinal não há comida, o corpo Gabor na borda despeja oscilação não-restringida. Isso
é *feature* (não alucina), e é a fronteira **seletor-não-gerador**.

**Gerar COM uma restrição, a seleção volta a vencer:** um prior (codebook), um vizinho (costura), uma
composição não-linear (oclusão). A geração dela **é** seleção sobre a restrição.

## Rotear = swarm + verificador (o inpaint poliglota)

Um SWARM de ~529 ladrilhos-seres seleciona de um codebook. Dois experts: **SUAVE** (difusão de Laplace)
e **TEXTURA** (quilt do codebook). Por componente do buraco, o **VERIFICADOR** reposiciona uma réplica
do buraco (mesma forma, textura da borda casada) numa região visível, cada expert a preenche, e roteia
pro que acerta o observável. Regra de segurança: **só troca do padrão (difusão) quando ele FALHA** (erro
alto na réplica) e a textura o bate com margem proporcional ao tamanho do buraco.

Resultado (média 3 seeds, PSNR no buraco, POLIGLOTA vs sempre-difusão vs oráculo):

| cenário | sint | photo_texture | photo_smooth |
|---|---|---|---|
| 2 buracos | 17.07 / 17.07 / 17.33 | **24.21 / 22.60 / 24.21** | 33.38 / 33.38 / 33.38 |
| ~40 buracos (--scatter) | 24.52 / 24.53 / 24.53 | 20.32 / 20.32 / 20.32 | 27.76 / 27.76 / 27.76 |

A poliglota **nunca cai abaixo da difusão** e **bate o oráculo em 5/6** (fecha +1.6 dB no photo_texture).

## Conhecimento de MODELO TREINADO nos seres

`inpaint --prior-json <cb.json>`: o codebook pode vir de um modelo de imagem treinado (SD-VAE via
`../bytebrain/research/extract_image_prior.py`). Achado honesto: o conhecimento importado **ajuda o
expert +1.0–1.7 dB exatamente quando o CONTEXTO é ESCASSO** (buraco grande), em todas as fotos; com muito
visível, a **auto-similaridade** (o codebook da própria imagem) vence. **Lei precisa:** quanto mais o ser
OBSERVA, menos precisa de prior importado; quanto menos observa, mais o prior treinado paga.

## Comandos

```bash
cargo run --release --bin agua-lab -- inpaint                          # seleção pura (a fronteira)
cargo run --release --bin agua-lab -- inpaint --prior                  # + prior (self-similaridade)
cargo run --release --bin agua-lab -- inpaint --poliglota [--scatter] [--img foto.png]
cargo run --release --bin agua-lab -- inpaint --poliglota --img foto.png --prior-json assets/prior_vae.json
cargo run --release --bin agua-lab -- generate                         # gerar por prior
cargo run --release --bin agua-lab -- occlusion                        # átomo não-linear
```

> Ver `JOURNAL.md` para a trilha completa (o que ganhou, o que morreu, e por quê).
