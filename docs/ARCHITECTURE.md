# Arquitetura — iara-água

## O modelo mental (mapa visão → código)

| Visão (`IARA_VISAO_SERES_COMO_AGUA.md`) | Código |
|---|---|
| Mundo / relevo (AGNÓSTICO) | `Field { dims, data }` N-dim flat (zero-médio; DC = 4 bytes de side-info) |
| Comida / energia | utilidade marginal: energia que o átomo explica do resíduo-sem-ele |
| Ser / átomo / alfinete | `Being<A: Atom> { atom, unit: Patch, amp, energy }` — Gabor1D/2D/3D |
| Reconstrução | `recon = Σ amp_i · ψ_i` (mantida incremental, patches esparsos) |
| Sobreviver = reduzir erro | `energy += food_gain·explicado − metabolism`; ≤0 → morre, patch volta pro resíduo |
| Água escoando | chuva de imigrantes em hotspots ∝ resíduo² + reprodução com mutação |
| Harmonizar com vizinhos | re-fit de amp a cada tick = relaxação Gauss-Seidel (átomos se acomodam) |
| Resultado | `population.json` — os parâmetros dos vivos SÃO a codificação |

## Invariantes (testados)

1. **Consistência:** `target − recon − residual ≈ 0` sempre (`invariant_error()`, rebuild anti-deriva a cada 500 ticks e antes de números finais).
2. **Determinismo:** mesma seed → mesmo mundo, bit a bit (RNG próprio xorshift64*, zero `rand`/`Instant` na física; teste `deterministic_same_seed`).
3. **Amplitude nunca evolui:** `amp* = ⟨r,ψ⟩/‖ψ‖²` em forma fechada (o ser "bebe"). Sobrevivência só esculpe GEOMETRIA — comparação limpa com matching pursuit.
4. **Metabolismo = dial taxa-distorção:** ser sobrevive sse explica > `metabolism/food_gain` de energia. Sobe o metabolismo → codificação enxuta; desce → detalhista.

## Agnosticismo de modalidade (implementado)

- `Patch` = lista de **runs contíguos no índice flat** — uma linha de imagem, um intervalo de áudio, uma linha de um frame de vídeo. `dot_field`/`dot_patch`/`add_scaled_patch` funcionam pra qualquer dimensionalidade (two-pointer sobre runs ordenados).
- `agua-life::World<A: Atom>` é genérico: a física NÃO conhece modalidade. O trait `Atom` (em `agua-atoms`) exige `rain_at(pos)`, `mutate`, `render_unit(dims)`, `freq()`. Nova modalidade = novo átomo, física intocada (o teste `survival_improves_psnr_audio` prova: mesma config resolve 1D).
- Baselines `random_atoms`/`matching_pursuit` também genéricos — comparação sempre na MESMA família de átomos do experimento.
- Hotspot = amostrador exato ∝ v² (prefix-sum + busca binária), N-dimensional por construção.
- IO por modalidade vive em `agua-lab` (dispatch) + `agua-telemetry` (PNG, WAV PCM16 reader/writer próprio, waveform PNG). Vídeo v0: alvo sintético; frames PNG na bateria.
- Texto/latente (F4): átomo em espaço de embedding com comida = casar com prior importado — o trait já comporta (Field = matriz de embeddings flat).

## Decisões de design

- **Rust headless puro, zero GPU na v0.** 64×64 roda em segundos na CPU; ciência primeiro, aceleração quando o perf virar gargalo real (F5). GPU/HIP só com profiling provando necessidade.
- **Patch unitário cacheado por ser** (`unit`, `n2`): feeding/refit custa O(suporte), não O(imagem). Morte/mutação = subtrair/somar patch.
- **Busca local com histerese** (`try_tune`, aceita se ≥2% melhor): mata thrashing (lição do rollback D2 do ByteBrain).
- **Baselines no mesmo crate-família:** random e MP usam os MESMOS Gabors/fit — diferença medida é só o algoritmo de busca, não a família de funções.
- **Config = TOML versionado + cópia no run dir**: qualquer resultado reproduz com `agua-lab run --config runs/<...>/config.toml --seed <s>`.
- **Reuso do ecossistema:** GA/niching viRão de `iara-swarm/iara-evolution` (F2); CfC portado com teste de paridade vs `Universe` (F2); prior byte-model via export offline Python→arquivo (F4, mesma ponte do `world_knowledge.json`). NÃO recriar engines.

## Fluxo de um tick (`World::step`)

```
1. tune   (a cada tune_every): mutação de geometria, aceita se explica ≥1.02×
2. feed   : refit amp (relaxação) + energia += food − metabolismo
3. death  : energia ≤ 0 → patch devolvido ao resíduo (vira comida)
4. repro  : energia ≥ threshold → filho mutado perto do pai (se pop < cap)
5. rain   (a cada rain_every): imigrantes em hotspots do resíduo
6. rebuild (a cada 500): recon do zero, mata deriva f32
```

## Rastreabilidade (agua-telemetry)

```
runs/<exp>/run-<unix>-s<seed>/
├── manifest.json      # seed, git hash, alvo, quando
├── config.toml        # cópia exata
├── metrics.jsonl      # série temporal: psnr, vivos, births/deaths, bytes, bandas de freq
├── population.json    # A CODIFICAÇÃO (todos os átomos vivos)
├── summary.json       # final + baselines random/MP/JPEG + wall time + invariant_error
└── snapshots/*.png    # target, recon_NNNN, residual_NNNN (resíduo assinado, 0=cinza)
JOURNAL.md             # 1 linha honesta por run (veredito bate/não-bate)
```

## Escalabilidade planejada (não implementada — YAGNI até F5)

- `rayon` no loop de feed (seres independentes dado o resíduo do tick).
- Grade espacial de vizinhança quando houver interação ser↔ser explícita.
- Quantização dos 7 params (f32→int8/int4) pra bytes honestos — liga com o codec DNA/semente (piso 768B é do codec, mesmo problema).
- Modalidade nova = implementar `render_unit` 1D (grão de áudio) — a física de `agua-life` não muda (é o teste de agnosticismo, F3).
