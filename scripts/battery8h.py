#!/usr/bin/env python3
"""BATERIA 8H — 200+ validações em várias frentes, com deadline e refill.

Frentes (não preso à metáfora dos seres — inclui matemática pura):
  leis      — c(σ), c(dof), universalidade(ruído), transição de fase θ*, cross-modal
  matematica— as MESMAS leis testadas FORA do sistema (numpy): EVT, AIC em regressão,
              Morozov/VisuShrink em wavelet Haar — se a teoria não valer lá, não é teoria
  denoise   — vs MP@oráculo, gauss@oráculo E wavelet clássico (baseline que faltava)
  inteligencia — adaptação a troca de sinal (velocidade de re-harmonização)
  performance — wall/evals/s/RSS vs tamanho e população
  memoria   — deriva numérica em runs longos (20k-100k ticks) + gradcheck em massa
  armazenamento — codec 4 perfis em populações reais + entropia ordem-0 (ganho de rANS)
  novas-leis — Lei da População (N=E/met), Lei do Espelho Espectral (pop ∝ espectro)

Cada resultado = 1 linha JSON em runs/_battery8h/results.jsonl (flush imediato).
Deadline: para de LANÇAR após DEADLINE_H; refill de seeds preenche o tempo.
"""
import json, os, subprocess, sys, time, math, itertools, random
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = f"{ROOT}/target/release/agua-lab"
OUT = f"{ROOT}/runs/_battery8h"
os.makedirs(OUT, exist_ok=True)
RESULTS = open(f"{OUT}/results.jsonl", "a", buffering=1)
LOG = open(f"{OUT}/progress.log", "a", buffering=1)

DEADLINE_H = float(os.environ.get("BATTERY_HOURS", "7.5"))
T0 = time.time()
def left(): return DEADLINE_H*3600 - (time.time()-T0)
def log(m):
    stamp = f"[{(time.time()-T0)/3600:5.2f}h]"
    print(stamp, m); LOG.write(f"{stamp} {m}\n")

def emit(group, row, wall):
    row = dict(row); row["_group"] = group; row["_wall_s"] = round(wall, 2)
    row["_t_h"] = round((time.time()-T0)/3600, 3)
    RESULTS.write(json.dumps(row) + "\n")

def probe(group, args):
    t = time.time()
    try:
        r = subprocess.run([BIN, "probe"] + [str(a) for a in args],
                           capture_output=True, text=True, timeout=3600)
        line = [l for l in r.stdout.strip().split("\n") if l.startswith("{")]
        if not line:
            emit(group, {"error": r.stderr[-400:] or "sem saida", "args": args}, time.time()-t); return
        emit(group, json.loads(line[-1]), time.time()-t)
    except Exception as e:
        emit(group, {"error": str(e)[:400], "args": [str(a) for a in args]}, time.time()-t)

def denoise_cmd(group, args):
    t = time.time()
    try:
        r = subprocess.run([BIN, "denoise"] + [str(a) for a in args],
                           capture_output=True, text=True, timeout=3600)
        line = [l for l in r.stdout.strip().split("\n") if l.startswith("[{")]
        if not line:
            emit(group, {"error": r.stderr[-400:], "args": args}, time.time()-t); return
        rows = json.loads(line[-1])
        d = {x["metodo"]: x for x in rows}
        emit(group, {"kind": "denoise", "args": [str(a) for a in args],
                     "ruidoso": d.get("ruidoso (nada)", {}).get("psnr_limpo"),
                     "gauss_or": d.get("gauss@oráculo", {}).get("psnr_limpo"),
                     "vida": d.get("SOBREVIVÊNCIA", {}).get("psnr_limpo"),
                     "vida_n": d.get("SOBREVIVÊNCIA", {}).get("atomos"),
                     "mp_or": d.get("MP@N-oráculo", {}).get("psnr_limpo"),
                     "mp_or_n": d.get("MP@N-oráculo", {}).get("atomos"),
                     "mp_moroz": d.get("MP@Morozov", {}).get("psnr_limpo")}, time.time()-t)
    except Exception as e:
        emit(group, {"error": str(e)[:400]}, time.time()-t)

def pymath(group, fn):
    t = time.time()
    try:
        emit(group, fn(), time.time()-t)
    except Exception as e:
        emit(group, {"error": str(e)[:400]}, time.time()-t)

# ───────────────────────── matemática pura (numpy, sem seres) ─────────────────────────
import numpy as np

def math_evt(M, seed):
    """teoria de valores extremos: E[max de M χ²₁] ≈ 2 ln M + O(lnln). É a base da Lei 1."""
    rng = np.random.default_rng(seed)
    xs = rng.normal(size=(4000, M))**2
    emax = xs.max(axis=1).mean()
    return {"kind": "math_evt", "M": M, "seed": seed, "e_max": float(emax),
            "pred_2lnM": float(2*np.log(M)), "ratio": float(emax/(2*np.log(M)))}

def math_aic(n, sigma, seed):
    """AIC em regressão polinomial: recupera a ordem verdadeira? (a Lei 2 fora dos seres)"""
    rng = np.random.default_rng(seed)
    x = np.linspace(-1, 1, n); true_k = 3
    y = 1.2*x**3 - 0.8*x + 0.3 + rng.normal(0, sigma, n)
    best_k, best_aic = None, 1e18
    for k in range(0, 9):
        c = np.polyfit(x, y, k); res = y - np.polyval(c, x)
        aic = n*np.log(np.mean(res**2)) + 2*(k+1)
        if aic < best_aic: best_aic, best_k = aic, k
    return {"kind": "math_aic", "n": n, "sigma": sigma, "seed": seed,
            "true_k": true_k, "chosen_k": int(best_k), "hit": bool(best_k == true_k)}

def haar2(img):
    a = img[0::2, 0::2]; b = img[0::2, 1::2]; c = img[1::2, 0::2]; d = img[1::2, 1::2]
    return (a+b+c+d)/2, (a-b+c-d)/2, (a+b-c-d)/2, (a-b-c+d)/2

def ihaar2(ll, lh, hl, hh):
    h, w = ll.shape; out = np.zeros((h*2, w*2))
    a = (ll+lh+hl+hh)/2; b = (ll-lh+hl-hh)/2; c = (ll+lh-hl-hh)/2; d = (ll-lh-hl+hh)/2
    out[0::2, 0::2] = a; out[0::2, 1::2] = b; out[1::2, 0::2] = c; out[1::2, 1::2] = d
    return out

def wavelet_denoise(noisy, sigma, levels=3, mode="universal"):
    """VisuShrink Haar (soft-threshold λ=σ√(2lnN)) — o baseline clássico que faltava."""
    coeffs = []; ll = noisy.copy()
    for _ in range(levels):
        ll, lh, hl, hh = haar2(ll); coeffs.append((lh, hl, hh))
    lam = sigma*np.sqrt(2*np.log(noisy.size))
    soft = lambda v: np.sign(v)*np.maximum(np.abs(v)-lam, 0)
    for i in range(len(coeffs)):
        lh, hl, hh = coeffs[i]; coeffs[i] = (soft(lh), soft(hl), soft(hh))
    for lh, hl, hh in reversed(coeffs):
        ll = ihaar2(ll, lh, hl, hh)
    return ll

def synth_image(size):
    s = float(size)
    x = np.arange(size)[None, :].astype(float); y = np.arange(size)[:, None].astype(float)
    TAU = 2*np.pi
    grad = 0.25*np.sin(TAU*x/s*0.5)*np.cos(TAU*y/s*0.5)
    blob = 0.45*np.exp(-(((x-s*0.3)**2)+((y-s*0.32)**2))/(2*(s*0.12)**2))
    stripes = np.where((x > s*0.55) & (y > s*0.55), 0.22*np.sin(TAU*0.22*x), 0.0)
    square = np.where((np.abs(x-s*0.7) < s*0.09) & (np.abs(y-s*0.25) < s*0.09), 0.3, 0.0)
    return np.clip(0.45+grad+blob+stripes+square, 0, 1)

def load_gray(path, size):
    from PIL import Image
    im = Image.open(path).convert("L").resize((size, size))
    return np.asarray(im, dtype=float)/255.0

def math_wavelet(target, size, sigma, seed):
    """baseline clássico no MESMO protocolo do nosso denoise (honestidade)."""
    img = synth_image(size) if target == "synthetic" else load_gray(f"{ROOT}/{target}", size)
    rng = np.random.default_rng(seed)
    noisy = img + rng.normal(0, sigma, img.shape)
    den = wavelet_denoise(noisy, sigma)
    p = lambda r: 10*np.log10(1.0/np.mean((img-r)**2))
    return {"kind": "math_wavelet", "target": target, "size": size, "sigma": sigma, "seed": seed,
            "psnr_noisy": float(p(noisy)), "psnr_visushrink": float(p(den))}

def math_entropy(run_glob):
    """entropia ordem-0 dos params quantizados de uma população real — quanto um
    entropy coder (rANS) ainda ganharia sobre os bits fixos do codec."""
    import glob
    runs = sorted(glob.glob(f"{ROOT}/{run_glob}/run-*/population.json"))
    if not runs: return {"kind": "math_entropy", "error": "sem populacao"}
    pop = json.load(open(runs[-1]))
    A = np.array([[a["atom"]["cx"], a["atom"]["cy"], a["atom"]["freq"],
                   a["atom"]["theta"] % np.pi, a["atom"]["phase"] % (2*np.pi),
                   a["atom"]["sigma"]] for a in pop])
    amp = np.array([a["amp"] for a in pop])
    out = {}
    tot_fixed, tot_h = 0.0, 0.0
    names = ["cx", "cy", "freq", "theta", "phase", "sigma", "amp"]
    for j, name in enumerate(names):
        v = amp if name == "amp" else A[:, j]
        bits = 10 if name in ("cx", "cy", "amp") else 6
        lo, hi = v.min(), v.max()
        q = np.round((v-lo)/max(hi-lo, 1e-9)*(2**bits-1)).astype(int)
        _, cnt = np.unique(q, return_counts=True)
        pr = cnt/cnt.sum()
        h = float(-(pr*np.log2(pr)).sum())
        out[name] = {"bits_fixos": bits, "entropia": round(h, 3)}
        tot_fixed += bits; tot_h += h
    return {"kind": "math_entropy", "n": len(pop), "campos": out,
            "bits_fixos_total": tot_fixed, "entropia_total": round(tot_h, 2),
            "ganho_rans_pct": round(100*(1-tot_h/tot_fixed), 1)}

# ───────────────────────── fila de tarefas ─────────────────────────
def build_queue():
    q = []  # (prioridade, grupo, callable)
    S5 = [11, 22, 33, 44, 55]
    S3 = [11, 22, 33]

    # LEIS — réplicas multi-seed dos achados fortes (agora com 5 seeds)
    for sg, sd in itertools.product([0.02, 0.05, 0.1, 0.2, 0.4], S5):
        q.append((1, "lei1_c_sigma", lambda sg=sg, sd=sd: probe("lei1_c_sigma",
            ["--kind","c","--sigma",sg,"--ticks",800,"--comp",8,"--seed",sd])))
    masks = [(0b111111,0),(0b111110,1),(0b111100,2),(0b111000,3),(0b110000,4),(0b100000,5),(0,6)]
    for (m, dof), sd in itertools.product(masks, S3):
        q.append((1, "lei2_dof", lambda m=m, sd=sd, dof=dof: probe("lei2_dof",
            ["--kind","c","--frozen",m,"--comp",1,"--pred",0.3,"--rain",0,"--ticks",2000,"--seed",sd])))
    for tipo, sg, sd in itertools.product(["gauss","laplace","uniforme","sal_0.01","sal_0.3"],[0.05,0.1,0.2],S3):
        q.append((2, "lei5_universal", lambda t=tipo, sg=sg, sd=sd: probe("lei5_universal",
            ["--kind","c","--noise",t,"--sigma",sg,"--ticks",800,"--comp",8,"--seed",sd])))
    for th, sd in itertools.product([0.3,0.38,0.44,0.48,0.52,0.56,0.62,0.7,0.8], S3):
        q.append((2, "lei6_fase", lambda th=th, sd=sd: probe("lei6_fase",
            ["--kind","c","--frozen",0b111111,"--comp",1,"--pred",th,"--rain",0,"--ticks",800,"--seed",sd])))

    # MATEMÁTICA PURA — as leis fora dos seres
    for M, sd in itertools.product([4,16,64,256,1024,4096], S3):
        q.append((1, "math_evt", lambda M=M, sd=sd: pymath("math_evt", lambda: math_evt(M, sd))))
    for n, sg, sd in itertools.product([50,200],[0.05,0.2,0.5],S5):
        q.append((1, "math_aic", lambda n=n, sg=sg, sd=sd: pymath("math_aic", lambda: math_aic(n, sg, sd))))

    # DENOISE — agora com o baseline wavelet clássico no mesmo protocolo
    alvos = ["synthetic","assets/photo_texture.png","assets/photo_smooth.png"]
    for alvo, sg, sd in itertools.product(alvos, [0.1,0.15,0.2,0.3], S3):
        q.append((2, "denoise", lambda a=alvo, sg=sg, sd=sd: denoise_cmd("denoise",
            ["--target",a,"--size",96,"--sigma",sg,"--seed",sd,"--ticks",800,"--auto"])))
        q.append((1, "math_wavelet", lambda a=alvo, sg=sg, sd=sd: pymath("math_wavelet",
            lambda: math_wavelet(a, 96, sg, sd))))
    for sg, sd in itertools.product([0.1,0.2],S3):  # áudio real
        q.append((2, "denoise_audio", lambda sg=sg, sd=sd: denoise_cmd("denoise_audio",
            ["--modality","audio","--target","assets/speech.wav","--size",8192,"--sigma",sg,"--seed",sd,"--ticks",800,"--auto","--k",25])))

    # γ(SNR) — mapa do k ótimo (2 modalidades)
    for k, sg, sd in itertools.product([10,20,30,40,60,90],[0.05,0.1,0.2,0.3],[11,22]):
        q.append((3, "gamma_map_img", lambda k=k, sg=sg, sd=sd: probe("gamma_map_img",
            ["--kind","morozov","--k",k,"--sigma",sg,"--seed",sd,"--ticks",800,"--target","assets/photo_texture.png"])))
    for k, sg, sd in itertools.product([8,15,25,40,60],[0.1,0.2],[11,22]):
        q.append((3, "gamma_map_aud", lambda k=k, sg=sg, sd=sd: probe("gamma_map_aud",
            ["--kind","morozov","--modality","audio","--target","assets/speech.wav","--size",8192,"--k",k,"--sigma",sg,"--seed",sd,"--ticks",800])))

    # NOVAS LEIS
    for met, alvo, sd in itertools.product([0.05,0.1,0.2,0.4,0.8], ["synthetic","assets/photo_texture.png"], S3):
        q.append((2, "lei_populacao", lambda m=met, a=alvo, sd=sd: probe("lei_populacao",
            ["--kind","equilibrium","--met",m,"--target",a,"--seed",sd,"--ticks",1500])))
    for alvo, sd in itertools.product(["synthetic","assets/photo_texture.png","assets/photo_smooth.png"], S3):
        q.append((2, "lei_espelho", lambda a=alvo, sd=sd: probe("lei_espelho",
            ["--kind","spectrum","--target",a,"--met",0.2,"--seed",sd,"--ticks",1500])))

    # INTELIGÊNCIA — adaptação
    for met, sd in itertools.product([0.3,0.5,1.0], S3):
        q.append((2, "adapt", lambda m=met, sd=sd: probe("adapt",
            ["--kind","adapt","--met",m,"--seed",sd,"--ticks-a",1500,"--ticks-b",600])))

    # PERFORMANCE + MEMÓRIA
    for size in [64,96,128,192,256]:
        for sd in [11,22]:
            q.append((2, "perf", lambda s=size, sd=sd: probe("perf",
                ["--kind","perf","--size",s,"--ticks",600,"--seed",sd])))
    for cap in [500,1500,3000,6000]:
        q.append((3, "perf_cap", lambda c=cap: probe("perf_cap",
            ["--kind","perf","--size",96,"--cap",c,"--met",0.02,"--ticks",600,"--seed",11])))
    for t, sd in itertools.product([20000,50000],[11,22]):
        q.append((3, "drift", lambda t=t, sd=sd: probe("drift",
            ["--kind","drift","--ticks",t,"--seed",sd])))
    for mod, sd in itertools.product(["image","audio"], S3):
        q.append((1, "gradcheck", lambda m=mod, sd=sd: probe("gradcheck",
            ["--kind","gradcheck","--modality",m,"--n",100,"--seed",sd])))

    # ARMAZENAMENTO — codec em populações reais + entropia
    for alvo, met, sd in itertools.product(["synthetic","assets/photo_texture.png"],[0.1,0.5],S3):
        q.append((2, "codec", lambda a=alvo, m=met, sd=sd: probe("codec",
            ["--kind","codec","--target",a,"--met",m,"--seed",sd,"--ticks",1200])))
    q.append((3, "entropia", lambda: pymath("entropia", lambda: math_entropy("runs/exp001-image-survival"))))

    q.sort(key=lambda x: x[0])
    return [(g, f) for _, g, f in q]

def refill(seed_base):
    """réplicas extras (seeds novas) das frentes-chave, até acabar o tempo."""
    S = [seed_base, seed_base+7]
    out = []
    for sg, sd in itertools.product([0.05,0.1,0.2], S):
        out.append(("lei1_c_sigma+", lambda sg=sg, sd=sd: probe("lei1_c_sigma",
            ["--kind","c","--sigma",sg,"--ticks",800,"--comp",8,"--seed",sd])))
    for (m, dof) in [(0b111111,0),(0b111100,2),(0,6)]:
        for sd in S:
            out.append(("lei2_dof+", lambda m=m, sd=sd: probe("lei2_dof",
                ["--kind","c","--frozen",m,"--comp",1,"--pred",0.3,"--rain",0,"--ticks",2000,"--seed",sd])))
    for alvo in ["synthetic","assets/photo_texture.png"]:
        for sd in S:
            out.append(("denoise+", lambda a=alvo, sd=sd: denoise_cmd("denoise",
                ["--target",a,"--size",96,"--sigma",0.2,"--seed",sd,"--ticks",800,"--auto"])))
            out.append(("lei_populacao+", lambda a=alvo, sd=sd: probe("lei_populacao",
                ["--kind","equilibrium","--met",0.2,"--target",a,"--seed",sd,"--ticks",1500])))
    for sd in S:
        out.append(("adapt+", lambda sd=sd: probe("adapt",
            ["--kind","adapt","--met",0.5,"--seed",sd,"--ticks-a",1500,"--ticks-b",600])))
        out.append(("gradcheck+", lambda sd=sd: probe("gradcheck",
            ["--kind","gradcheck","--modality","image","--n",100,"--seed",sd])))
    return out

def main():
    random.seed(0)
    queue = build_queue()
    log(f"fila base: {len(queue)} validações · deadline {DEADLINE_H}h · workers 6")
    done = 0
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(f): g for g, f in queue}
        seed_base = 100
        while futs:
            for fut in as_completed(list(futs)):
                g = futs.pop(fut)
                done += 1
                if done % 10 == 0:
                    log(f"{done} concluídas · {len(futs)} na fila · restam {left()/3600:.2f}h")
                # refill enquanto houver tempo
                if len(futs) < 8 and left() > 1200:
                    extra = refill(seed_base)
                    seed_base += 17
                    for g2, f2 in extra:
                        futs[ex.submit(f2)] = g2
                    log(f"refill +{len(extra)} (seeds base {seed_base})")
                break
    log(f"FIM: {done} validações em {(time.time()-T0)/3600:.2f}h")
    # sumário rápido por grupo
    counts = {}
    for l in open(f"{OUT}/results.jsonl"):
        try:
            r = json.loads(l); counts[r.get("_group","?")] = counts.get(r.get("_group","?"),0)+1
        except Exception:
            pass
    log("por grupo: " + json.dumps(counts, ensure_ascii=False))

if __name__ == "__main__":
    main()
