# F5 — Escala e performance (veredito honesto)

## Curva taxa-distorção: JPEG esmaga (a compressão é jogo perdido, agora medido)

Com o codec CORRIGIDO (bug da amplitude consertado) e bytes quantizados REAIS,
alvo foto real 96²:

| método | bytes | PSNR |
|---|---|---|
| sobrevivência Q8 (met=0.8) | 2462 | 25.64 |
| sobrevivência Q8 (met=0.3) | 6633 | 28.12 |
| sobrevivência fp16 (met=0.1) | 36095 | 34.75 |
| **JPEG q30** | **1259** | **31.12** |
| **JPEG q70** | **1982** | **35.13** |

**JPEG ganha por 5–18× em bytes no mesmo PSNR.** A F1 já dizia que compressão era
jogo perdido; agora está medido COM quantização — não é mais "os params não estão
quantizados". Mesmo quantizado, perde feio. Motivo: 7 params/átomo × muitos átomos
vs DCT+RLE+Huffman do JPEG. Sem surpresa, mas honestamente cravado.

**Conclusão definitiva: a sobrevivência NÃO é um compressor.** O valor dela está
nas LEIS, no denoise (seleção de N), na geração e na poda — não em bytes/PSNR.

## Performance e escala (medida, não otimizada às cegas)

| tamanho | wall (ticks) | evals/s | RSS pico |
|---|---|---|---|
| 96² | 51.6s (400) | 11.7k | ~42 MB |
| 128² | 100.6s (400) | 9.4k | ~50 MB |
| 192² | 176.7s (400) | 6.8k | — |
| 256² | 245.4s (300) | 4.8k | 145 MB |

Escala ~O(N·suporte) por tick. RSS modesto (dezenas de MB até 192²).

## Paralelização: rayon na chuva, honesto

Paralelizei a chuva competitiva (K candidatos, render+fit puros) preservando o
determinismo (sorteio sequencial do RNG, avaliação paralela). Teste de
determinismo por seed segue passando.

**MAS não move o hot path do denoise.** O gargalo é o feed sequencial
(Gauss-Seidel: cada ser ajusta amp contra o resíduo JÁ atualizado pelos anteriores
— dependência serial) e o tune por-gradiente. A chuva é fração pequena. Para
paralelizar o hot path seria preciso feed JACOBI (todos contra o mesmo snapshot,
aplicam juntos) — mudança de algoritmo com questão de convergência aberta. Regra
do ecossistema: medir antes de otimizar; medi, e o ganho do rayon aqui é pequeno.

## Conclusão de F5

Compressão fechada (negativa, definitiva). Performance medida e honesta. O rayon
é correto mas o hot path é sequencial por design. GPU não se justifica: o valor da
tese não está em throughput de reconstrução, está nas leis e nos regimes
não-lineares onde o guloso não joga.
