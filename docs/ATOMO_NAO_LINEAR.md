# O átomo não-linear — o primeiro terreno próprio da tese-água

> 2026-07-10. `agua-lab occlusion` (`crates/agua-lab/src/occlusion.rs`).
> A visão §4 sempre disse: o diferencial da sobrevivência é o átomo NÃO-linear,
> onde gradiente e guloso-linear não alcançam. Este é o primeiro teste disso.

## Por que aqui o matching pursuit não joga

- Soma linear: `recon(x) = Σ amp_i·ψ_i(x)`. O MP funciona porque o resíduo é
  linear e a amplitude tem forma fechada `⟨r,ψ⟩/‖ψ‖²`.
- **Oclusão** (algoritmo do pintor): `recon(x)` = valor do disco MAIS À FRENTE
  em x. Não há amplitude linear no valor (é um `max` por profundidade), tirar um
  disco tampado não muda nada, e não existe `⟨resíduo, átomo⟩` que preveja a
  melhoria. **O MP nem se define.**

A sobrevivência, sim: a comida de um disco = quanto a MSE piora se ele for
removido (utilidade marginal na composição não-linear). Isso é mensurável sem
gradiente e sem linearidade. Disco tampado dá comida ~0 → morre de fome.

## Resultado (alvos de oclusão: K discos opacos empilhados)

| K discos | oclusão (sobrevivência) | melhor linear (MP-soma) | Δ |
|---|---|---|---|
| 6 | 22–25 dB | 20–27 dB | **+1.4 / +1.7 / −1.3** |
| 12 | 21–24 dB | 21–24 dB | +0.34 / +2.03 / −0.07 |
| 20 | 20–22 dB | 20–20 dB | +0.01 / −0.11 / **+2.23** |

**6/9 vitórias, Δ médio +0.69 dB.**

## Veredito honesto

- **O que está provado:** a sobrevivência OTIMIZA uma composição não-diferenciável
  (max por oclusão) só com pressão de comer — onde o MP não tem como tentar. Esse
  é o ponto da tese, e ele se sustenta.
- **O que NÃO está provado:** que a oclusão é dramaticamente melhor. A margem é
  modesta e some com cenas densas (K=20), onde a soma linear aproxima bem uma
  imagem cheia. Em cenas esparsas (K=6, bordas duras isoladas) a oclusão ganha
  claro, porque a soma vaza/ringa nas bordas.
- **Ressalva de desenho:** um alvo feito por oclusão é naturalmente melhor servido
  por oclusão — isso é esperado e não é o mérito. O mérito é a OTIMIZABILIDADE sem
  gradiente, não a margem de PSNR.

## Onde isto aponta (a vitória real da tese, ainda por vir)

A oclusão é o caso mais simples de átomo não-linear. O caso que importa de verdade
(visão §4): o átomo = um MODELO importado (byte-model, pedaço de outro sinal, regra).
Aí a comida = "quanto casar com o prior melhora o todo", a composição pode ser
qualquer coisa, e a GERAÇÃO (não reconstrução) vira possível — amostrar sinais
novos casando com o prior. Isso é a F4 e continua intocado. Mas a peça-chave —
que a sobrevivência otimiza o não-diferenciável — está agora demonstrada.
