# Números como frações — veredito (IARA_NUMEROS_COMO_FRACOES.md)

> Testado 2026-07-09 em `scripts/fracoes.py`, sobre os **3000 átomos vivos** de um run real (PSNR f32 = 40.92 dB). Métrica end-to-end (PSNR após decodificar e re-renderizar), como o doc exige na §4. Predições registradas **antes** de rodar.

Seguindo a §1 do doc: **racional livre p/q não foi testado** (estoura o denominador — já sabido).

## Resultados (mesmo orçamento de bits)

| esquema | B/átomo | PSNR | Δ vs f32 |
|---|---|---|---|
| float32 (referência) | 28.00 | 40.92 | — |
| **float16** | 14.00 | **40.82** | **−0.10** |
| uniforme 8 bits (amp range fixo ±4) | 7.00 | 7.06 | −33.85 |
| **uniforme 8 bits (amp range adaptativo)** | 7.00 | **25.26** | −15.65 |
| CONTROLE: K-means **livre** K=16 (4 bits) | 3.50 | **17.48** | −23.44 |
| (A) alfabeto 16 **frações** q≤16 (4 bits) | 3.50 | 14.70 | −26.21 |
| (A) alfabeto 16 frações q≤8 (4 bits) | 3.50 | −12.06 | −52.98 |
| (B) fração contínua, 3 termos | 10.78 | 9.26 | −31.66 |
| (B) fração contínua, 5 termos | 19.12 | 27.36 | −13.56 |

## Vereditos

**(B) fração contínua / Stern-Brocot — ❌ REFUTADA.** Predição feita antes: perderia, porque os parâmetros geométricos são contínuos e não têm estrutura racional. Confirmado com folga: a 10.78 B/átomo entrega 9.26 dB, contra 26.90 dB da uniforme a 10.50 B. O float16 a 14 B entrega 40.82 dB. A hierarquia grosso→fino do Stern-Brocot é conceitualmente linda, mas os convergentes gastam bits em *exatidão racional* que ninguém pediu.

**(A) alfabeto de frações — ❌ REFUTADA pelo controle.** O alfabeto ganha muito da uniforme no mesmo orçamento (14.70 vs −3.03 dB a 3.5 B/átomo). **Mas o controle mostra que o mérito é do codebook, não das frações:** K-means *livre* com o mesmo K entrega **17.48 dB**, melhor que os 14.70 do K-means restrito a racionais. Restringir centróides a racionais só pode piorar — e piora. O que funciona é **quantização adaptada à distribuição**, que é product quantization clássico.

## O que a investigação rendeu (achado colateral, esse sim real)

**Bug de projeto no `agua-codec`.** O intervalo de amplitude estava fixo em ±4, mas as amplitudes reais são ~1e-3 — a quantização gastava todos os níveis em valores que nunca ocorrem. Com **escala adaptativa** (range real da população, 8 bytes no header):

- 8 bits/param: **7.06 dB → 25.26 dB** (+18.2 dB)
- 10 bits/param: 21.31 → 26.89 dB (+5.6 dB)

Corrigido em `agua-codec` (`AmpScale`), com teste de regressão (`amp_minuscula_sobrevive_q8`).

**Recomendação prática para o codec:** amplitude em fp16 (ou range adaptativo), geometria em 8–10 bits. A geometria é robusta; a amplitude não.

## Ressalvas

- 1 população, 1 seed, alvo sintético. As conclusões (B refutada, controle vence A) são grandes o bastante para não dependerem disso, mas o número exato do bug (+18.2 dB) é desse run.
- Não testei (C) — exatidão em recorrência longa do CfC. O doc já a classifica como nicho marginal, e (A)/(B) não animaram.
- Não testei nos **pesos de modelo** (bytebrain/Qwen), que é a §2(A) original. Aqui o alvo foram os parâmetros dos átomos, que é o gancho §2(B)-4 (o mais alinhado com a água). Para pesos, a distribuição é diferente (concentrada em zero) e o alfabeto pode se sair melhor — mas o controle K-means livre continuaria valendo como baseline obrigatório.
