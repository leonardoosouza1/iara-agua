#!/usr/bin/env bash
# denoise_sweep.sh — bateria de restauração por sobrevivência (zero-shot).
# uso: scripts/denoise_sweep.sh <tag> <PAR> <<'EOF'
#      --target X --sigma S --metabolism M --seed N ...   (uma linha por run)
#      EOF
set -euo pipefail
TAG="$1"; PAR="${2:-6}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/target/release/agua-lab"
OUT="$ROOT/runs/_denoise/$TAG"; mkdir -p "$OUT"
mapfile -t LINES
CMDS="$OUT/cmds.txt"; : > "$CMDS"
i=0
for line in "${LINES[@]}"; do
  [[ -z "${line// }" || "$line" == \#* ]] && continue
  i=$((i+1))
  echo "$BIN denoise $line > $OUT/run_$i.txt 2>&1" >> "$CMDS"
done
echo "▶ $TAG: $i runs, paralelismo $PAR"
xargs -a "$CMDS" -P "$PAR" -I{} sh -c '{}' || true
echo "✔ $TAG"
# agrega: última linha de cada run é o JSON
python3 - "$OUT" <<'PY'
import json,sys,glob,os,re,statistics as st
out=sys.argv[1]
rows=[]
for f in sorted(glob.glob(out+"/run_*.txt")):
    txt=open(f).read()
    hdr=re.search(r"== DENOISE (\S+) sigma=(\S+) seed=(\S+) met=(\S+) grad=(\S+) ==",txt)
    js=[l for l in txt.strip().split("\n") if l.startswith("[{")]
    if not hdr or not js: continue
    mod,sig,seed,met,grad=hdr.groups()
    tgt=re.search(r"--target (\S+)",open(out+"/cmds.txt").read().split("\n")[int(re.search(r"run_(\d+)",f).group(1))-1])
    data=json.loads(js[-1])
    d={r["metodo"]:r for r in data}
    rows.append(dict(mod=mod,sigma=float(sig),seed=int(seed),met=float(met),
        alvo=os.path.basename(tgt.group(1)) if tgt else "?",
        ruidoso=d["ruidoso (nada)"]["psnr_limpo"],
        gauss=d.get("gauss@oráculo",{}).get("psnr_limpo",float('nan')),
        vida=d["SOBREVIVÊNCIA"]["psnr_limpo"], vida_n=d["SOBREVIVÊNCIA"]["atomos"],
        mp_or=d["MP@N-oráculo"]["psnr_limpo"], mp_or_n=d["MP@N-oráculo"]["atomos"],
        mp_mo=d["MP@Morozov"]["psnr_limpo"], mp_mo_n=d["MP@Morozov"]["atomos"]))
json.dump(rows,open(out+"/rows.json","w"),indent=1)

key=lambda r:(r["alvo"],r["sigma"],r["met"])
groups={}
for r in rows: groups.setdefault(key(r),[]).append(r)
print(f"\n{'alvo':<20}{'sig':>5}{'met':>6}{'n':>3} | {'ruid':>6}{'gauss':>7}{'VIDA':>7}{'MP_or':>7}{'MP_mz':>7} | {'vida-MPor':>10}{'N_vida':>7}{'N_MPor':>7}")
for k in sorted(groups):
    g=groups[k]; m=lambda f:st.mean(r[f] for r in g)
    d=m("vida")-m("mp_or")
    star=" <<<" if d>0 else ""
    print(f"{k[0]:<20}{k[1]:>5}{k[2]:>6}{len(g):>3} | {m('ruidoso'):>6.2f}{m('gauss'):>7.2f}{m('vida'):>7.2f}{m('mp_or'):>7.2f}{m('mp_mo'):>7.2f} | {d:>+10.2f}{m('vida_n'):>7.0f}{m('mp_or_n'):>7.0f}{star}")
PY
