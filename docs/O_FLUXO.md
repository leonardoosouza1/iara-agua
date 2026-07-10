# O Fluxo — um seletor por sobrevivência, quatro problemas

> A síntese prática de tudo. Denoise, seleção de variáveis, poda de rede e
> compressão de sinal são o MESMO procedimento com nomes diferentes. Este é o
> fluxo, e o harness que o executa.

## O procedimento único

```
ENTRADA:  um alvo Y (sinal, resposta, ativação-de-saída)
          um banco de candidatos {ψ_j} (átomos, variáveis, neurônios)
          uma escala de ruído σ (medida, não chutada)

LAÇO (sobrevivência):
  1. resíduo r ← Y − Σ (vivos)
  2. cada candidato propõe: quanto de r ele explica?  g_j = ⟨r, ψ̂_j⟩²
  3. NASCE quem explica mais (chuva nos hotspots do resíduo)
  4. AJUSTA a geometria por gradiente (Adam) OU mutação
  5. MORRE quem explica menos que o custo:  g_j < 2σ²·ln M
  6. repete até o resíduo bater σ²  (Morozov — parada automática)

SAÍDA:    o conjunto dos vivos É a resposta
          o número deles EMERGE — ninguém escolhe N
```

O passo 5 é a Lei do Garimpo. O passo 6 é Morozov. Juntos eliminam **todos**
os hiperparâmetros: σ vem do MAD do sinal, M vem do orçamento de busca, e N sai
do equilíbrio. Não há threshold para ajustar.

## O mesmo fluxo, quatro problemas

| Problema | Alvo Y | Candidatos ψ | σ | O que emerge |
|---|---|---|---|---|
| **Denoise** | imagem/áudio ruidoso | átomos Gabor | MAD do detalhe | reconstrução limpa |
| **Seleção de variáveis** | vetor resposta | colunas de X | ruído da regressão | o suporte verdadeiro |
| **Poda de rede** | saída da camada | neurônios ocultos | MAD do resíduo | sub-rede esparsa |
| **Codec** | sinal-alvo | átomos Gabor | (tende a 0) | os parâmetros = os bytes |

Muda só o que são Y e ψ. A física — nascer, ajustar, morrer por custo, parar no
ruído — é literalmente o mesmo código (`agua-life::World<A: Atom>`).

## Onde cada peça mora no repo

```
agua-core         Field N-dim · Patch esparso · RNG determinístico · hotspot sampler
agua-atoms        trait Atom + Gabor 1D/2D/3D · gradiente analítico validado
agua-life         O LAÇO. World<A>::step() — nascer/ajustar/morrer/chover/predar
agua-metrics      baselines obrigatórios: random, matching pursuit, JPEG
agua-codec        a população → bytes (.agua) — fp16, Q12, Q8
agua-telemetry    cada run rastreável: manifest+git, metrics.jsonl, snapshots
agua-lab          o harness: run · denoise --auto · probe (10 kinds) · law2..law9
```

## O harness de baterias (reusável)

Toda medição primitiva é uma invocação `agua-lab probe --kind <X>` que emite **uma
linha JSON**. As baterias são só orquestração dessas linhas:

```
scripts/battery8h.py       fila de ~400 validações + refill de seeds até deadline
scripts/battery8h_extra.py fila pesada (escala, deriva longa, áudio real)
scripts/analise8h.py       agrega o results.jsonl → veredito por frente
```

Kinds do probe (cada um 1 JSON, componível):

```
c            fator de busca c = E[explicado]/σ²  (leis 1,2,4,5,6)
morozov      denoise met=k·σ̂² · PSNR + discrepância
equilibrium  N ≈ E_capturada/met  (Lei da População)
spectrum     assinatura de frequência da população
gradcheck    gradiente analítico vs diferenças finitas (100%)
drift        run longo · deriva do invariante (< 1e-3 até 50k ticks)
adapt        troca de sinal · velocidade de re-harmonização
perf         wall · evals/s · RSS de pico
codec        população real → .agua nos 4 perfis → PSNR
gamma_fdr    γ = met_ótimo/met_c como quantil de falsa-sobrevivência
```

## O protocolo de honestidade (o que impede o auto-engano)

Extraído das falsificações reais desta sessão — cada regra tem um cadáver por trás:

1. **Multi-seed sempre.** A "predação paga +0.29 dB" era ruído de config; morreu em outro alvo.
2. **Baseline com oráculo.** O denoise só vale porque o MP escolhe seu N *vendo a imagem limpa* e mesmo assim perde.
3. **Controle contra o trivial.** O "Espelho Espectral 0.992" batia igual com o alvo *errado* — a correlação era trivial. Só a assinatura discriminativa (18/21) sobrou.
4. **Regime de validade.** A Lei do Garimpo falha abaixo do limiar de detecção; documentado, não escondido.
5. **Matar a tautologia.** "Poda metabólica ≈ OMP" é trivial (é OMP com parada). O resultado real é a seleção a orçamento fixo.
6. **σ oráculo separa estimador de método.** O "Morozov substitui k" era viés do estimador, não a regra.

> Se a evidência não distingue duas hipóteses, o fluxo diz que não distingue.
> Sete hipóteses minhas morreram assim nesta sessão. É o que dá confiança nas quatro que sobraram.
