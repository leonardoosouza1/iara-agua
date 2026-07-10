# iara-água

> Uma população de átomos famintos come o resíduo de um sinal. Quem explica mais
> energia do que custa existir, sobrevive. A soma dos vivos **é** a resposta.
> Ninguém programa a solução — só a pressão de sobreviver.

Workspace Rust headless, determinístico, rastreável. Materializa a visão
`../IARA_VISAO_SERES_COMO_AGUA.md`. Este README é o mapa honesto: **o que ganha,
o que perde, e as leis que sobraram depois de eu matar as minhas próprias
hipóteses.**

## Em uma imagem

```
   sinal (imagem/áudio/vídeo)
        │
        ▼
   ┌─────────────────────────────────────────┐
   │  🌧️ nasce onde falta explicar (resíduo)  │
   │  🍽️ come o pedaço que cobre               │  ← o laço
   │  🔧 ajusta a geometria (gradiente/mutação)│
   │  💀 morre se explica < 2σ²·ln(alternativas)│
   │  ⏹️ para quando o resíduo = nível do ruído │  (Morozov)
   └─────────────────────────────────────────┘
        │
        ▼
   os sobreviventes SÃO a codificação
   (e quantos são EMERGE — ninguém escolhe N)
```

## O que GANHA ✅

| | resultado |
|---|---|
| **Denoise sem hiperparâmetro** | σ estimado do sinal, N emerge da lei. Bate matching pursuit com N-oráculo **737/737**, e VisuShrink clássico **81/81** (+3.8 dB). `agua-lab denoise --auto` |
| **Geração por prior** | codebook de imagem real → textura nova (não cópia), costura −87% vs aleatório. `agua-lab generate` |
| **Átomo não-linear** | oclusão, onde o matching pursuit *nem se define*. 6/9 vs melhor linear. `agua-lab occlusion` |
| **Roteador poliglota** | línguas (textura/liso/borda) competem, cada região vai pra língua certa, roteador **emerge**. +5.5 dB na borda. `agua-lab polyglot` |

## O que PERDE ❌ (medido, sem maquiar)

| | verdade |
|---|---|
| **Compressão** | JPEG ganha 5–18× em bytes, mesmo com o codec quantizado. Provado duas vezes. A sobrevivência **não é um compressor**. |
| **Poda de LLM** | Wanda ganha na perplexidade real do ByteBrain (1.19 vs 1.60 bpb). Reconstrução de camada era proxy enganoso. |

## As LEIS (cada uma com predição feita ANTES da medida)

- **Metabolismo crítico** — `met_c = σ²·(1 + 2·dof)`. Um átomo que não busca captura 1·σ²; cada grau de liberdade custa ≈2σ². **Dois por parâmetro é a penalidade do AIC.** A multidão redescobre seleção de modelo por fome. (`docs/LAWS.md`, R²=0.98, n=990)
- **Lei do Garimpo** — `penalidade = 2σ²·ln M`. É o RIC (Foster & George 1994). Vale pra finanças (8192 backtests → Sharpe 3.8 do acaso), seleção de variáveis e — como teoria — poda. (`docs/LEI_DO_GARIMPO.md`)
- **Morozov emergente** — a população para exatamente no nível do ruído sem saber σ. Cross-modal. (`docs/LAWS.md`)
- **Sinapse = matriz de Gram** — `G_ij = ⟨ψ_i,ψ_j⟩` é o peso sináptico; o acoplamento alterna excita/inibe com a distância = **chapéu-mexicano** (córtex visual). A susceptibilidade a um disparo dispara 1600× no regime colinear (criticalidade) — a mesma colinearidade que infla o metabolismo. (`docs/LEIS_SINAPSE.md`)
- **Transição de fase da harmonia** — a população fica coerente, não ortogonal; a predação impõe harmonia numa transição abrupta em θ*≈0.5. (`docs/LAWS.md`)

## O que morreu (o que dá confiança no resto)

γ universal · Donoho-Johnstone · predação que paga · Morozov-substitui-k · frações contínuas · alfabeto de frações · espelho espectral (era trivial) · ressonância linear · poda em LLM. **Nove hipóteses minhas, falsificadas com teste honesto.**

## Crates

```
agua-core       Field N-dim · Patch esparso · RNG determinístico · hotspot sampler
agua-atoms      trait Atom + Gabor 1D/2D/3D · gradiente analítico validado
agua-life       O LAÇO (genérico): World<A>::step — nascer/ajustar/morrer/predar
agua-metrics    baselines: random · matching pursuit · JPEG
agua-codec      população → bytes (.agua): fp16, Q12, Q8
agua-telemetry  cada run rastreável: manifest+git · metrics.jsonl · PNG/WAV
agua-lab        harness: run · denoise · generate · occlusion · polyglot · synapse
                · resonance · probe · law2..law9
```

## Rodar

```bash
cargo test --release                                    # 33 testes
cargo run --release --bin agua-lab -- denoise --target synthetic --sigma 0.2 --auto
cargo run --release --bin agua-lab -- law2              # a lei do metabolismo
cargo run --release --bin agua-lab -- synapse           # as leis da sinapse
cargo run --release --bin agua-lab -- polyglot          # o roteador que emerge
```

## Docs

- `docs/O_FLUXO.md` — o seletor por sobrevivência unificado (denoise=seleção=poda)
- `docs/LAWS.md` · `docs/LEI_DO_GARIMPO.md` · `docs/LEIS_SINAPSE.md` — as leis
- `docs/PONTE_LINGUAS.md` — a ponte entre representações (tese IARA-Router)
- `docs/ROADMAP.md` — F0..F6 todos fechados
- `docs/EXPERIMENT_PROTOCOL.md` — as regras de honestidade (cada uma tem um cadáver)
- `docs/VISAO_BYTE_HORMONAL.md` — a próxima fronteira (fora da curva, aterrada)

## O que isto é, honestamente

Uma **contribuição conceitual**: leis novas medidas (fome = AIC, garimpo = RIC,
sinapse = Gram) + duas demos que funcionam (denoise, geração). **Não** é um
compressor nem um podador que bate o estado-da-arte — e onde não é, está escrito.
A disciplina de matar o próprio falso-positivo é o método, não o acidente.
