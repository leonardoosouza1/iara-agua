#!/usr/bin/env bash
# collect.sh — agrega os summary.json de um lote em tabela honesta (média±spread por variante).
# uso: scripts/collect.sh <prefixo-do-exp>     ex: scripts/collect.sh L1
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PREFIX="$1"

printf "%-28s %6s  %8s %8s %8s  %7s %7s  %7s  %s\n" \
  VARIANTE seeds "psnr_med" "spread" "vs_MP" "vivos" "bytes" "wall_s" "bandas(grosso→fino)"

for d in "$ROOT"/runs/${PREFIX}-*/; do
  [[ -d "$d" ]] || continue
  name=$(basename "$d")
  files=("$d"run-*/summary.json)
  [[ -e "${files[0]}" ]] || continue
  jq -s --arg n "$name" -r '
    map(select(.final_psnr_db != null)) as $r
    | ($r | map(.final_psnr_db)) as $p
    | ($r | map(.baselines.matching_pursuit[-1].psnr_db // 0)) as $mp
    | ($r | map(.baselines.random[-1].psnr_db // 0)) as $rd
    | ($p | add / length) as $pm
    | ($mp | add / length) as $mpm
    | ($rd | add / length) as $rdm
    | (($p | max) - ($p | min)) as $spread
    | ($r | map(.alive) | add / length) as $alive
    | ($r | map(.bytes) | add / length) as $bytes
    | ($r | map(.wall_secs) | add / length) as $wall
    | ($r | map(.invariant_error) | max) as $inv
    | [$n, ($r|length), $pm, $spread, ($pm-$mpm), $alive, $bytes, $wall, $rdm, $inv]
    | @tsv' "${files[@]}" 2>/dev/null | while IFS=$'\t' read -r n c pm sp vsmp alive bytes wall rdm inv; do
      # bandas do último metrics.jsonl da variante
      last=$(ls -d "$d"run-* | tail -1)
      bands=$(tail -1 "$last/metrics.jsonl" 2>/dev/null | jq -c '.freq_bands' 2>/dev/null || echo "-")
      warn=""
      awk -v i="$inv" 'BEGIN{exit !(i>1e-3)}' && warn=" ⚠INV"
      awk -v p="$pm" -v r="$rdm" 'BEGIN{exit !(p<r+0.5)}' && warn="$warn ⚠NÃO-BATE-RANDOM"
      printf "%-28s %6s  %8.2f %8.2f %+8.2f  %7.0f %7.0f  %7.1f  %s%s\n" \
        "${n#$PREFIX-}" "$c" "$pm" "$sp" "$vsmp" "$alive" "$bytes" "$wall" "$bands" "$warn"
    done
done
