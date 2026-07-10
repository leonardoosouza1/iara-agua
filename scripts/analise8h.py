#!/usr/bin/env python3
"""Análise da bateria 8h. Honestidade: multi-seed, spread, e o teste
anti-circularidade do Espelho Espectral (contra a FFT do ALVO)."""
import json, sys, math, statistics as st, glob, os
from collections import defaultdict
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = [json.loads(l) for l in open(f"{ROOT}/runs/_battery8h/results.jsonl")]
ok = [r for r in rows if "error" not in r]
print(f"# BATERIA 8H — {len(rows)} validações, {len(rows)-len(ok)} erros\n")

def group(g): return [r for r in ok if r.get("_group") == g]
def mean(xs): return st.mean(xs) if xs else float("nan")
def sd(xs): return st.pstdev(xs) if len(xs) > 1 else 0.0

OUT = {}

# ── LEI 1: c independe de σ ────────────────────────────────────────────
g = [r for r in ok if r.get("kind") == "c" and r.get("_group", "").startswith("lei1")]
by = defaultdict(list)
for r in g: by[round(r["sigma"], 3)].append(r["c"])
print("## LEI 1 — c = met_c/σ² independe de σ?")
cs_all = []
for s in sorted(by):
    v = by[s]; cs_all += v
    print(f"  σ={s:<6} c = {mean(v):6.2f} ± {sd(v):4.2f}   (n={len(v)})")
cv = 100*sd([mean(by[s]) for s in by])/mean([mean(by[s]) for s in by])
print(f"  → CV entre σ = {cv:.1f}%  (n total {len(cs_all)})")
OUT["lei1"] = {"cv_pct": cv, "n": len(cs_all), "c_mean": mean(cs_all)}

# ── LEI 2: c = a·dof + b ───────────────────────────────────────────────
g = group("lei2_dof")
def dof_of(frozen): return 6 - bin(int(frozen)).count("1")
by = defaultdict(list)
for r in g: by[dof_of(r["frozen"])].append(r["c"])
xs, ys = [], []
print("\n## LEI 2 — c = a·dof + b  (a≈2 = penalidade do AIC; b=1 = âncora)")
for d in sorted(by):
    print(f"  dof={d}  c = {mean(by[d]):6.2f} ± {sd(by[d]):4.2f}  (n={len(by[d])})")
    for v in by[d]: xs.append(d); ys.append(v)
if len(set(xs)) > 2:
    A = np.polyfit(xs, ys, 1); pred = np.polyval(A, xs)
    r2 = 1 - np.sum((np.array(ys)-pred)**2)/np.sum((np.array(ys)-mean(ys))**2)
    print(f"  → c = {A[0]:.3f}·dof + {A[1]:.3f}   (R²={r2:.4f}, n={len(xs)})")
    OUT["lei2"] = {"a": float(A[0]), "b": float(A[1]), "r2": float(r2), "n": len(xs)}

# ── LEI 5: universalidade do ruído ─────────────────────────────────────
g = group("lei5_universal")
by = defaultdict(list)
for r in g: by[r["noise"]].append(r["c"])
print("\n## LEI 5 — c depende da FORMA do ruído? (CLT diz: não)")
base = mean(by.get("gauss", [1]))
uni = {}
for t in sorted(by):
    m = mean(by[t]); uni[t] = 100*(m-base)/base
    print(f"  {t:<12} c = {m:6.2f} ± {sd(by[t]):4.2f}   ({uni[t]:+6.1f}% vs gauss)  n={len(by[t])}")
OUT["lei5"] = uni

# ── LEI 6: transição de fase θ* ────────────────────────────────────────
g = group("lei6_fase")
by = defaultdict(list)
for r in g: by[round(r["pred"], 2)].append(r["c"])
print("\n## LEI 6 — transição de fase da harmonia (θ*)")
prev, theta_star = None, None
for t in sorted(by):
    m = mean(by[t])
    fase = "harmônica" if m < 4 else "colinear"
    salto = ""
    if prev and m/prev > 1.8: salto = " ⟵ SALTO"; theta_star = t
    print(f"  θ={t:<5} c = {m:6.2f} ± {sd(by[t]):4.2f}  {fase}{salto}")
    prev = m
OUT["lei6"] = {"theta_star": theta_star}

# ── DENOISE: vs MP@oráculo, gauss@oráculo, VisuShrink ──────────────────
den = [r for r in ok if r.get("kind") == "denoise" and r.get("vida") and r.get("mp_or")]
wav = {(r["target"], round(r["sigma"],3), r["seed"]): r["psnr_visushrink"] for r in ok if r.get("kind")=="math_wavelet"}
print(f"\n## DENOISE — {len(den)} runs")
w_mp = sum(1 for r in den if r["vida"] > r["mp_or"])
w_g  = sum(1 for r in den if r.get("gauss_or") and r["vida"] > r["gauss_or"])
print(f"  vs MP@N-oráculo    : {w_mp}/{len(den)} vitórias · Δ médio {mean([r['vida']-r['mp_or'] for r in den]):+.2f} dB")
gg = [r for r in den if r.get("gauss_or")]
print(f"  vs gauss@oráculo   : {w_g}/{len(gg)} vitórias · Δ médio {mean([r['vida']-r['gauss_or'] for r in gg]):+.2f} dB")
pares = []
for r in den:
    a = r["args"]
    try:
        key = (a[a.index("--target")+1], round(float(a[a.index("--sigma")+1]),3), int(a[a.index("--seed")+1]))
    except Exception: continue
    if key in wav: pares.append((r["vida"], wav[key]))
if pares:
    w_w = sum(1 for v,x in pares if v > x)
    print(f"  vs VisuShrink Haar : {w_w}/{len(pares)} vitórias · Δ médio {mean([v-x for v,x in pares]):+.2f} dB  (clássico!)")
    OUT["denoise"] = {"vs_mp": [w_mp, len(den)], "vs_gauss": [w_g, len(gg)], "vs_wavelet": [w_w, len(pares)],
                      "d_mp": mean([r['vida']-r['mp_or'] for r in den]),
                      "d_wav": mean([v-x for v,x in pares])}

# denoise em ESCALA (a lei diz: met_c independe da resolução)
esc = group("denoise_escala")
if esc:
    print("\n  denoise em ESCALA (mesma lei, sem recalibrar k):")
    by = defaultdict(list)
    for r in esc:
        a = r["args"]; size = int(a[a.index("--size")+1])
        by[size].append((r["vida"], r["mp_or"]))
    for s in sorted(by):
        v = by[s]
        print(f"    {s}×{s}: vida {mean([x for x,_ in v]):.2f} dB  vs MP@or {mean([y for _,y in v]):.2f}  → {sum(1 for x,y in v if x>y)}/{len(v)}")

# ── LEI DA POPULAÇÃO: N ≈ E_capturada/met ─────────────────────────────
eq = [r for r in ok if r.get("kind") == "equilibrium" and r["alive"] > 10]
print(f"\n## LEI DA POPULAÇÃO — N ≈ E_capturada / met   ({len(eq)} runs)")
by = defaultdict(list)
for r in eq: by[round(r["met"],3)].append(r["alive"]/r["n_pred"])
for m in sorted(by):
    print(f"  met={m:<6} N/(E/met) = {mean(by[m]):.3f} ± {sd(by[m]):.3f}  (n={len(by[m])})")
allr = [r["alive"]/r["n_pred"] for r in eq]
# capados no teto não valem (cap=3000)
livres = [r for r in eq if r["alive"] < 2900]
print(f"  todos: {mean(allr):.3f} ± {sd(allr):.3f}   |  sem saturar o cap: {mean([r['alive']/r['n_pred'] for r in livres]):.3f} ± {sd([r['alive']/r['n_pred'] for r in livres]):.3f} (n={len(livres)})")
OUT["lei_pop"] = {"ratio": mean([r['alive']/r['n_pred'] for r in livres]), "sd": sd([r['alive']/r['n_pred'] for r in livres])}

# ── LEI DO ESPELHO: teste ANTI-CIRCULARIDADE contra FFT do alvo ────────
def synth_image(size):
    s = float(size)
    x = np.arange(size)[None,:].astype(float); y = np.arange(size)[:,None].astype(float)
    TAU = 2*np.pi
    grad = 0.25*np.sin(TAU*x/s*0.5)*np.cos(TAU*y/s*0.5)
    blob = 0.45*np.exp(-(((x-s*0.3)**2)+((y-s*0.32)**2))/(2*(s*0.12)**2))
    stripes = np.where((x > s*0.55)&(y > s*0.55), 0.22*np.sin(TAU*0.22*x), 0.0)
    square = np.where((np.abs(x-s*0.7) < s*0.09)&(np.abs(y-s*0.25) < s*0.09), 0.3, 0.0)
    return np.clip(0.45+grad+blob+stripes+square, 0, 1)

def target_spectrum_bands(target, size, edges=(0.03,0.06,0.12,0.24)):
    if target.startswith("synthetic"): img = synth_image(size)
    else:
        from PIL import Image
        img = np.asarray(Image.open(f"{ROOT}/{target}").convert("L").resize((size,size)), dtype=float)/255.0
    img = img - img.mean()
    F = np.fft.fftshift(np.abs(np.fft.fft2(img))**2)
    fy = np.fft.fftshift(np.fft.fftfreq(size))[:,None]
    fx = np.fft.fftshift(np.fft.fftfreq(size))[None,:]
    rad = np.sqrt(fx**2 + fy**2)
    bands, lo = [], 0.0
    for e in list(edges)+[10.0]:
        bands.append(float(F[(rad>=lo)&(rad<e)].sum())); lo = e
    return bands

sp = [r for r in ok if r.get("kind") == "spectrum"]
print(f"\n## LEI DO ESPELHO — população espelha o espectro do sinal? ({len(sp)} runs)")
cor_recon, cor_fft = [], []
def pearson(a,b):
    a=np.array(a,dtype=float); b=np.array(b,dtype=float)
    if a.std()==0 or b.std()==0: return float("nan")
    return float(np.corrcoef(a,b)[0,1])
cache = {}
for r in sp:
    p = r["band_pop"]
    cor_recon.append(pearson(p, r["band_energy"]))
    key = (r["target"], r["size"])
    if key not in cache:
        try: cache[key] = target_spectrum_bands(r["target"], r["size"])
        except Exception: cache[key] = None
    if cache[key]: cor_fft.append(pearson(p, cache[key]))
cr = [c for c in cor_recon if not math.isnan(c)]
cf = [c for c in cor_fft if not math.isnan(c)]
print(f"  corr(pop, energia da PRÓPRIA recon) = {mean(cr):.3f} ± {sd(cr):.3f}   ← parcialmente circular")
print(f"  corr(pop, espectro FFT do ALVO)     = {mean(cf):.3f} ± {sd(cf):.3f}   ← teste honesto (n={len(cf)})")
OUT["espelho"] = {"circular": mean(cr), "honesto": mean(cf), "n": len(cf)}

# ── MATEMÁTICA PURA (fora dos seres) ──────────────────────────────────
evt = [r for r in ok if r.get("kind")=="math_evt"]
print(f"\n## MATEMÁTICA PURA — as leis fora do sistema ({len(evt)} EVT)")
by = defaultdict(list)
for r in evt: by[r["M"]].append(r["ratio"])
print("  EVT: E[max de M χ²₁] / (2·ln M)   → deve → 1 quando M cresce")
for M in sorted(by): print(f"    M={M:<6} razão = {mean(by[M]):.3f}")
aic = [r for r in ok if r.get("kind")=="math_aic"]
if aic:
    hits = sum(r["hit"] for r in aic)
    over = sum(1 for r in aic if r["chosen_k"] > r["true_k"])
    under = sum(1 for r in aic if r["chosen_k"] < r["true_k"])
    print(f"  AIC em regressão polinomial: acerta a ordem {hits}/{len(aic)} · superestima {over} · subestima {under}")
    OUT["math"] = {"aic_hit": hits, "aic_n": len(aic), "aic_over": over, "aic_under": under}

# ── INTELIGÊNCIA: adaptação ───────────────────────────────────────────
ad = [r for r in ok if r.get("kind")=="adapt"]
if ad:
    print(f"\n## ADAPTAÇÃO — sinal troca; a população reescoa ({len(ad)} runs)")
    print(f"  queda no switch: {mean([r['psnr_a']-r['psnr_on_switch'] for r in ad]):.2f} dB")
    print(f"  recuperação    : {mean([r['recovery_db'] for r in ad]):+.2f} dB  (±{sd([r['recovery_db'] for r in ad]):.2f})")
    frac = [ (r['psnr_final']-r['psnr_on_switch'])/max(r['psnr_a']-r['psnr_on_switch'],1e-6) for r in ad]
    print(f"  fração do dano recuperada: {100*mean(frac):.0f}%")
    OUT["adapt"] = {"recovery": mean([r['recovery_db'] for r in ad]), "frac": mean(frac)}

# ── PERFORMANCE / MEMÓRIA ─────────────────────────────────────────────
pf = [r for r in ok if r.get("kind")=="perf"]
if pf:
    print(f"\n## PERFORMANCE ({len(pf)} runs)")
    by = defaultdict(list)
    for r in pf: by[r["size"]].append(r)
    print(f"  {'lado':>6} {'evals/s':>10} {'ticks/s':>9} {'RSS MB':>8} {'PSNR':>7} {'vivos':>7}")
    for s in sorted(by):
        v = by[s]
        print(f"  {s:>6} {mean([x['evals_per_s'] for x in v]):>10.0f} {mean([x['ticks_per_s'] for x in v]):>9.1f} "
              f"{mean([x['rss_peak_mb'] for x in v]):>8.1f} {mean([x['psnr'] for x in v]):>7.2f} {mean([x['alive'] for x in v]):>7.0f}")
    OUT["perf"] = {str(s): {"evals_s": mean([x['evals_per_s'] for x in by[s]]),
                            "rss_mb": mean([x['rss_peak_mb'] for x in by[s]])} for s in sorted(by)}

dr = [r for r in ok if r.get("kind")=="drift"]
if dr:
    print(f"\n## MEMÓRIA / ESTABILIDADE NUMÉRICA ({len(dr)} runs longos)")
    for r in sorted(dr, key=lambda r: r["ticks"]):
        print(f"  {r['ticks']:>7} ticks · invariante = {r['invariant_error']:.2e} · PSNR {r['psnr']:.2f} · {r['wall_s']:.0f}s")
    OUT["drift"] = {"max_inv": max(r["invariant_error"] for r in dr), "max_ticks": max(r["ticks"] for r in dr)}

gc = [r for r in ok if r.get("kind")=="gradcheck"]
if gc:
    p = sum(r["pass"] for r in gc); f = sum(r["fail"] for r in gc)
    print(f"\n## CÁLCULO — gradiente analítico vs diferenças finitas: {p}/{p+f} átomos ({100*p/(p+f):.1f}%)")
    OUT["gradcheck"] = {"pass": p, "total": p+f}

# ── ARMAZENAMENTO ─────────────────────────────────────────────────────
cd = [r for r in ok if r.get("kind")=="codec"]
if cd:
    print(f"\n## ARMAZENAMENTO — codec .agua ({len(cd)} populações)")
    prof = defaultdict(list)
    for r in cd:
        for p in r["profiles"]:
            prof[p["profile"]].append((p["bytes"]/r["alive"], p["psnr"] - r["psnr_f32_world"]))
    print(f"  {'perfil':>8} {'B/átomo':>9} {'Δ PSNR vs f32':>15}")
    for name in ["f32","fp16","q12","q8"]:
        v = prof[name]
        print(f"  {name:>8} {mean([b for b,_ in v]):>9.2f} {mean([d for _,d in v]):>+15.2f} dB")
    OUT["codec"] = {n: {"bpa": mean([b for b,_ in prof[n]]), "dpsnr": mean([d for _,d in prof[n]])} for n in ["f32","fp16","q12","q8"]}
ent = [r for r in ok if r.get("kind")=="math_entropy" and "campos" in r]
if ent:
    e = ent[-1]
    print(f"  entropia ordem-0 dos params: {e['entropia_total']:.1f} bits vs {e['bits_fixos_total']:.0f} fixos → rANS ganharia {e['ganho_rans_pct']}%")
    OUT["entropia"] = {"ganho_rans_pct": e["ganho_rans_pct"]}

json.dump(OUT, open(f"{ROOT}/runs/_battery8h/summary.json","w"), indent=1)
print(f"\n→ resumo em runs/_battery8h/summary.json")
