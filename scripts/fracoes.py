#!/usr/bin/env python3
"""Testa a ideia IARA_NUMEROS_COMO_FRACOES.md nos parâmetros dos átomos vivos.

Compara, no MESMO orçamento de bits, quatro representações dos 7 parâmetros de
cada Gabor2D, e mede o efeito END-TO-END (PSNR da reconstrução após decodificar,
não só o MSE dos parâmetros — como o doc exige na §4):

  (0) float32                 — referência sem perda
  (1) fp16
  (2) uniforme de b bits      — o que o codec já faz (Q12/Q8)
  (A) alfabeto de frações     — K-means restrito a racionais de denominador ≤ D
  (B) fração contínua truncada — convergentes de Stern-Brocot, bits variáveis (Elias-gamma)

PREDIÇÃO (feita antes de rodar): os parâmetros geométricos são contínuos e não têm
estrutura racional, então (B) deve PERDER para a uniforme. (A) só ganha se a
distribuição for concentrada — o que ocorre em `sigma` (satura no teto) e `amp`
(concentrada perto de 0).
"""
import json, sys, math
from fractions import Fraction
import numpy as np

RUN = sys.argv[1]
pop = json.load(open(f"{RUN}/population.json"))
man = json.load(open(f"{RUN}/manifest.json"))
H, W = man["world"]["dims"]
DC = man["world"]["dc"]

# alvo: regenera o sintético exatamente como o Rust (targets.rs)
def synthetic_image(size):
    s = float(size)
    x = np.arange(size)[None, :].astype(float)
    y = np.arange(size)[:, None].astype(float)
    TAU = 2*np.pi
    grad = 0.25*np.sin(TAU*x/s*0.5)*np.cos(TAU*y/s*0.5)
    blob = 0.45*np.exp(-(((x-s*0.3)**2)+((y-s*0.32)**2))/(2*(s*0.12)**2))
    stripes = np.where((x > s*0.55) & (y > s*0.55), 0.22*np.sin(TAU*0.22*x), 0.0)
    square = np.where((np.abs(x-s*0.7) < s*0.09) & (np.abs(y-s*0.25) < s*0.09), 0.3, 0.0)
    return np.clip(0.45+grad+blob+stripes+square, 0, 1)

target_dc = synthetic_image(W)
target = target_dc - target_dc.mean()

YY, XX = np.mgrid[0:H, 0:W]

def render(atoms, amps):
    rec = np.zeros((H, W))
    for a, amp in zip(atoms, amps):
        cx, cy, f, th, ph, sg = a
        r = math.ceil(2.5*sg)
        x0, x1 = max(0, int(cx)-r), min(W, int(cx)+r+1)
        y0, y1 = max(0, int(cy)-r), min(H, int(cy)+r+1)
        if x1 <= x0 or y1 <= y0: continue
        dx = XX[y0:y1, x0:x1]-cx; dy = YY[y0:y1, x0:x1]-cy
        env = np.exp(-(dx*dx+dy*dy)/(2*sg*sg))
        car = np.cos(2*np.pi*f*(dx*np.cos(th)+dy*np.sin(th))+ph)
        rec[y0:y1, x0:x1] += amp*env*car
    return rec

def psnr(rec):
    mse = np.mean((target-rec)**2)
    return 99.0 if mse <= 1e-12 else 10*np.log10(1.0/mse)

A = np.array([[a["atom"]["cx"], a["atom"]["cy"], a["atom"]["freq"], a["atom"]["theta"],
               a["atom"]["phase"], a["atom"]["sigma"]] for a in pop])
AMP = np.array([a["amp"] for a in pop])
N = len(pop)

# limites por parâmetro (o decodificador conhece — vão no header, custo O(1))
LO = np.array([0.0, 0.0, 0.015, 0.0, 0.0, 1.2])
HI = np.array([float(W), float(H), 0.45, np.pi, 2*np.pi, 24.0])
AMP_R = 4.0
# BUG ENCONTRADO: amp real é ~1e-3; range fixo ±4 desperdiça todos os níveis.
# Escala adaptativa: usa o range REAL (2 floats no header, custo O(1)).
AMP_LO, AMP_HI = float(AMP.min()), float(AMP.max())

base = psnr(render(A, AMP))
print(f"população: {N} átomos · PSNR f32 = {base:.2f} dB\n")

def eval_scheme(nome, A2, AMP2, bits_total):
    p = psnr(render(A2, AMP2))
    bpa = bits_total/N/8
    print(f"{nome:<34} {bpa:>8.2f} B/átomo {bits_total/8/1024:>9.2f} KB {p:>9.2f} dB  {p-base:>+8.2f}")
    return bpa, p

print(f"{'esquema':<34} {'bytes/átomo':>8} {'total':>14} {'PSNR':>12} {'Δ vs f32':>10}")
eval_scheme("float32 (referência)", A, AMP, N*7*32)

# fp16
def f16(x): return np.float32(np.float16(x))
eval_scheme("float16", f16(A), f16(AMP), N*7*16)

# uniforme de b bits (o codec Q*)
for b in (12, 10, 8, 6, 4):
    q = np.round((A-LO)/(HI-LO)*(2**b-1)).clip(0, 2**b-1)
    Ad = LO + q/(2**b-1)*(HI-LO)
    qa = np.round((AMP+AMP_R)/(2*AMP_R)*(2**b-1)).clip(0, 2**b-1)
    AMPd = -AMP_R + qa/(2**b-1)*(2*AMP_R)
    eval_scheme(f"uniforme {b}b (amp range FIXO ±4)", Ad, AMPd, N*7*b)
    # o mesmo, com escala adaptativa da amp
    qa2 = np.round((AMP-AMP_LO)/(AMP_HI-AMP_LO)*(2**b-1)).clip(0, 2**b-1)
    AMPd2 = AMP_LO + qa2/(2**b-1)*(AMP_HI-AMP_LO)
    eval_scheme(f"uniforme {b}b (amp range ADAPTATIVO)", Ad, AMPd2, N*7*b)

# (A) alfabeto de frações: K-means restrito a racionais de denominador <= D
def alfabeto_fracoes(vals, K, D):
    lo, hi = vals.min(), vals.max()
    if hi-lo < 1e-9: return np.full_like(vals, lo), 0
    x = (vals-lo)/(hi-lo)                       # normaliza pra [0,1]
    c = np.quantile(x, np.linspace(0, 1, K))    # init por quantis
    for _ in range(25):
        idx = np.abs(x[:, None]-c[None, :]).argmin(1)
        for k in range(K):
            m = idx == k
            if m.any(): c[k] = x[m].mean()
        # PROJETA cada centróide no racional mais próximo com q <= D
        c = np.array([float(Fraction(v).limit_denominator(D)) for v in c])
    idx = np.abs(x[:, None]-c[None, :]).argmin(1)
    return lo + c[idx]*(hi-lo), idx

def alfabeto_livre(vals, K):
    """controle: K-means SEM restrição racional. Se empatar com o racional,
    o mérito é do codebook adaptativo, não das frações."""
    lo, hi = vals.min(), vals.max()
    if hi-lo < 1e-9: return np.full_like(vals, lo)
    x = (vals-lo)/(hi-lo)
    c = np.quantile(x, np.linspace(0, 1, K))
    for _ in range(25):
        idx = np.abs(x[:, None]-c[None, :]).argmin(1)
        for k in range(K):
            m = idx == k
            if m.any(): c[k] = x[m].mean()
    idx = np.abs(x[:, None]-c[None, :]).argmin(1)
    return lo + c[idx]*(hi-lo)

for K in (16, 8):
    bits = int(np.ceil(np.log2(K)))
    Ad = np.zeros_like(A)
    for j in range(6):
        Ad[:, j] = alfabeto_livre(A[:, j], K)
    AMPd = alfabeto_livre(AMP, K)
    eval_scheme(f"CONTROLE: K-means LIVRE K={K} = {bits} bits", Ad, AMPd, N*7*bits)

for K, D in ((16, 8), (16, 16), (8, 8)):
    bits = int(np.ceil(np.log2(K)))
    Ad = np.zeros_like(A)
    for j in range(6):
        Ad[:, j], _ = alfabeto_fracoes(A[:, j], K, D)
    AMPd, _ = alfabeto_fracoes(AMP, K, D)
    eval_scheme(f"alfabeto {K} frações (q≤{D}) = {bits} bits", Ad, AMPd, N*7*bits)

# (B) fração contínua truncada, bits variáveis (Elias-gamma nos termos)
def elias_gamma_bits(n):
    n = max(1, int(n))
    return 2*int(math.floor(math.log2(n)))+1

def frac_continua(x, termos):
    """x em [0,1]; devolve (aprox, bits) truncando em `termos` termos."""
    fr = Fraction(float(x)).limit_denominator(10**9)
    a, bits, cf = [], 0, []
    v = fr
    for _ in range(termos):
        ai = math.floor(v)
        cf.append(ai)
        bits += elias_gamma_bits(ai+1)
        v = v - ai
        if v == 0: break
        v = 1/v
    # reconstroi
    val = Fraction(cf[-1])
    for ai in reversed(cf[:-1]):
        if val == 0: val = Fraction(ai)
        else: val = ai + 1/val
    return float(val), bits

for termos in (2, 3, 4, 5):
    Ad = np.zeros_like(A); total_bits = 0
    for j in range(6):
        x = (A[:, j]-LO[j])/(HI[j]-LO[j])
        out = [frac_continua(v, termos) for v in x]
        Ad[:, j] = LO[j] + np.array([o[0] for o in out])*(HI[j]-LO[j])
        total_bits += sum(o[1] for o in out)
    xa = (AMP+AMP_R)/(2*AMP_R)
    outa = [frac_continua(v, termos) for v in xa]
    AMPd = -AMP_R + np.array([o[0] for o in outa])*(2*AMP_R)
    total_bits += sum(o[1] for o in outa)
    eval_scheme(f"fração contínua, {termos} termos", Ad, AMPd, total_bits)

print("\nPREDIÇÃO registrada antes de rodar: (B) fração contínua PERDE da uniforme,")
print("porque os parâmetros geométricos são contínuos e não têm estrutura racional.")
print("CONTROLE: se K-means LIVRE ≥ alfabeto de frações, o mérito é do codebook,")
print("não do fato de os valores serem racionais 'bonitos'.")
