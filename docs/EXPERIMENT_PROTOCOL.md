# Protocolo de Experimento — regras de honestidade

> Extraído das lições reais do ecossistema (falso-positivo H=8 morto em 5 seeds;
> word-ok 0.79 que era ruído; span invertido por memorização). NÃO são opcionais.

1. **Multi-seed obrigatório:** nenhuma conclusão com 1 seed. Mínimo 3; vitória fina (<1 dB) exige 5. Reportar média ± spread.
2. **Baseline sempre:** todo run compara com random-atoms (piso) e matching pursuit (teto guloso) com o MESMO nº de átomos e MESMA família, + JPEG em bytes reais. Não bater o random = ruído, reportar como ruído.
3. **Matar o próprio falso-positivo:** achou vitória? Rode mais seeds/janelas ANTES de anunciar. O journal registra a morte do falso-positivo também.
4. **Held-out:** conclusões sobre generalização usam alvos NÃO usados na calibração de hiperparâmetros.
5. **Wall time declarado:** todo resultado diz quanto tempo levou (`summary.json:wall_secs`).
6. **Números exatos, nunca redondos inventados.** Copiar do summary.json, não de memória.
7. **Invariante checado:** `invariant_error` no summary deve ser < 1e-3; maior = bug, resultado inválido.
8. **Um run = um registro:** nada de rodar 20 vezes e reportar o melhor. Toda tentativa vai pro JOURNAL.md (a linha automática cuida disso — não apagar linhas).
9. **Comparações justas de bytes:** átomo f32 = 28 B. Ao comparar com JPEG, dizer explicitamente que os params não estão quantizados (até F5).
10. **Negativo é resultado:** fase que falha o aceite fecha com veredito negativo documentado, não com silêncio.
