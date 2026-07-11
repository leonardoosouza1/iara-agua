# Verificação independente do guia §4 (cérebro por órgãos)

> `scripts/brain_organs_test.py`. Reimplementei o cérebro-por-órgãos mínimo
> (sensor→spikes→valor, plasticidade de 3 fatores via REINFORCE) e MEDI os 4
> achados centrais do `GUIA_CEREBRO_VIVO_PARA_RUST.md §4`. 5 seeds, held-out.

Acaso = 0.25.

| claim | medido | veredito |
|---|---|---|
| 1. aprende e generaliza | held-out 0.907 ± 0.020 | ✅ (guia: 0.72; mais alto, mesma história) |
| 2. dopamina essencial (ablação) | 0.284 ± 0.125 | ✅ vira acaso sem RPE (guia: 0.49) |
| 3. ternário-no-repouso ≈ float | 0.885 vs 0.907 (−0.02) | ✅ precisão por propósito |
| 4. contexto anti-esquecimento | sem gate: A 0.54 / B 0.20; com gate: A 0.91 / B 0.91 | ✅ o gate preserva os dois |

**Bug honesto que peguei no caminho:** minha 1ª versão usava seleção GULOSA (sem
exploração) — não aprendia (0.27 ≈ acaso). Aprendizado por recompensa PRECISA de
exploração (amostrar a ação, REINFORCE com baseline). Não era o claim do guia que
falhava — era meu código. Consertado, tudo confirma.

**Conclusão:** os 4 pilares do porte pra Rust estão validados de forma
independente. A plasticidade de 3 fatores aprende, a dopamina é o que esculpe, o
ternário-no-repouso não custa nada, e o gate por contexto resolve o esquecimento
catastrófico. O guia é sólido como fundação.
