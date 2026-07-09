#!/usr/bin/env bash
# sweep.sh — gera configs a partir de um template + overrides e roda em paralelo.
#
# uso:  scripts/sweep.sh <lote> <template.toml> <PARALLEL> <<'EOF'
#       name=<tag> KEY=VALUE KEY=VALUE ...      (uma variante por linha)
#       EOF
#
# Chaves aceitas: qualquer `key = ...` presente no template (life/atoms/world/experiment).
# Cada variante × cada seed do template vira um run com [experiment].name = <lote>-<tag>-s<seed>.
set -euo pipefail

LOTE="$1"; TEMPLATE="$2"; PAR="${3:-6}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GEN="$ROOT/runs/_configs/$LOTE"
mkdir -p "$GEN"
BIN="$ROOT/target/release/agua-lab"

# lê variantes do stdin
mapfile -t VARIANTS

CMDS="$GEN/cmds.txt"; : > "$CMDS"
for line in "${VARIANTS[@]}"; do
  [[ -z "${line// }" || "$line" == \#* ]] && continue
  tag=""; cfg="$GEN/tmp.toml"; cp "$TEMPLATE" "$cfg"
  for kv in $line; do
    k="${kv%%=*}"; v="${kv#*=}"
    if [[ "$k" == "name" ]]; then tag="$v"; continue; fi
    # substitui a PRIMEIRA ocorrência de `k = ...` (chaves são únicas nos templates)
    sed -i "0,\|^${k} = .*|s||${k} = ${v}|" "$cfg"
  done
  [[ -z "$tag" ]] && { echo "variante sem name=: $line" >&2; exit 1; }
  sed -i "0,\|^name = .*|s||name = \"${LOTE}-${tag}\"|" "$cfg"
  final="$GEN/${tag}.toml"; mv "$cfg" "$final"
  # 1 processo por seed (paralelismo real; cada run é determinístico)
  for seed in $(grep -oP '^seeds = \[\K[^]]*' "$final" | tr -d ' ' | tr ',' ' '); do
    echo "$BIN run --config $final --seed $seed --out $ROOT" >> "$CMDS"
  done
done

echo "▶ lote $LOTE: $(wc -l < "$CMDS") runs, paralelismo $PAR"
xargs -a "$CMDS" -P "$PAR" -I{} sh -c '{} > /dev/null 2>&1 || echo "FALHOU: {}"'
echo "✔ lote $LOTE concluído"
