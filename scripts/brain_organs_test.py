#!/usr/bin/env python3
"""Verifica INDEPENDENTEMENTE os achados centrais do GUIA_CEREBRO_VIVO §4.
Cérebro mínimo por órgãos: sensor→spikes→valor, aprende por plasticidade de
3 fatores (pré × dopamina/RPE × modulador). Tudo numpy, 5 seeds, held-out.

Testa 4 claims do guia:
  1. aprende e generaliza (held-out > acaso)
  2. dopamina é essencial (sem RPE → acaso)
  3. ternário-no-repouso ≈ float (precisão por propósito)
  4. contexto gateia → sem esquecimento catastrófico (2 domínios)
"""
import numpy as np

D, H, C = 8, 32, 4   # features, hidden, classes

def make_task(seed):
    r = np.random.default_rng(seed)
    W1 = r.normal(size=(D, H))                 # sensor→hidden (fixo, aleatório)
    rule = r.normal(size=(H, C))               # regra verdadeira (recompensa)
    def gen(n, rr):
        X = rr.integers(0, 2, (n, D)).astype(float)
        Hd = np.tanh(X @ W1)                    # spikes (contínuo p/ simplicidade)
        y = (Hd @ rule).argmax(1)              # classe correta = argmax da regra
        return X, Hd, y
    return W1, gen

def train_eval(seed, dopamine=True, ternary_rest=False, ctx_gate=False, two_domains=False):
    r = np.random.default_rng(seed*7+1)
    W1, gen = make_task(seed)
    n_ctx = 2 if (ctx_gate or two_domains) else 1
    w2 = np.zeros((n_ctx, H, C))
    lr = 0.05
    baseline = np.zeros(n_ctx)   # baseline de recompensa (freio/serotonina simples)
    def gen_dom(n, dom, rr):
        X, Hd, y = gen(n, rr)
        if dom == 1: y = (y + 1) % C
        return X, Hd, y
    doms = [0, 1] if two_domains else [0]
    for dom in doms:
        ctx = min(dom if (ctx_gate or two_domains) else 0, n_ctx-1)
        if two_domains and not ctx_gate: ctx = 0
        for _ in range(6000):
            X, Hd, y = gen_dom(1, dom, r)
            h = Hd[0]
            v = w2[ctx].T @ h                      # valor por classe
            # EXPLORA: amostra ação do softmax (não greedy)
            pv = np.exp(v - v.max()); pv /= pv.sum()
            a = int(r.choice(C, p=pv))
            reward = 1.0 if a == y[0] else 0.0
            # DOPAMINA = RPE (reward - baseline); sem dopamina = sinal fixo
            rpe = (reward - baseline[ctx]) if dopamine else 0.02
            baseline[ctx] = 0.99*baseline[ctx] + 0.01*reward
            # 3 fatores: pré(h) × dopamina(rpe) × (empurra a ação amostrada)
            grad = -pv.copy(); grad[a] += 1.0     # ∂logπ(a)
            w2[ctx] += lr * rpe * np.outer(h, grad)
    def acc(dom):
        ctx = min(dom if (ctx_gate or two_domains) else 0, n_ctx-1)
        Wq = w2[ctx].copy()
        if ternary_rest:
            sc = np.abs(w2[ctx]).max(0, keepdims=True) + 1e-9
            Wq = np.sign(w2[ctx]) * sc
        rr = np.random.default_rng(seed*99+dom+1)
        X, Hd, y = gen_dom(2000, dom, rr)
        return ((Hd @ Wq).argmax(1) == y).mean()
    if two_domains:
        return acc(0), acc(1)
    return acc(0)

seeds = [11,22,33,44,55]
def avg(f): 
    vs=[f(s) for s in seeds]; return np.mean(vs), np.std(vs)

print("VERIFICAÇÃO INDEPENDENTE do guia §4 (cérebro por órgãos, 5 seeds)\n")
print(f"acaso = {1/C:.2f}\n")

m,s = avg(lambda sd: train_eval(sd))
print(f"1. aprende e generaliza (held-out):  {m:.3f} ± {s:.3f}   {'✓ > acaso' if m>0.4 else '✗'}")

m2,s2 = avg(lambda sd: train_eval(sd, dopamine=False))
print(f"2. SEM dopamina (ablação):           {m2:.3f} ± {s2:.3f}   {'✓ cai pro acaso' if m2<0.35 else '✗ não caiu'}")

m3,s3 = avg(lambda sd: train_eval(sd, ternary_rest=True))
print(f"3. ternário-no-repouso:              {m3:.3f} ± {s3:.3f}   {'✓ ≈ float' if abs(m3-m)<0.1 else '✗ perde'}")

# 4. dois domínios: sem gate (esquece) vs com gate
noA,noB = avg(lambda sd: train_eval(sd, two_domains=True, ctx_gate=False)[0]), avg(lambda sd: train_eval(sd, two_domains=True, ctx_gate=False)[1])
gA = avg(lambda sd: train_eval(sd, two_domains=True, ctx_gate=True)[0])
gB = avg(lambda sd: train_eval(sd, two_domains=True, ctx_gate=True)[1])
print(f"\n4. esquecimento catastrófico (treina A, depois B):")
print(f"   SEM gate:  A={noA[0]:.3f}  B={noB[0]:.3f}   {'← A esquecido' if noA[0]<0.4 else ''}")
print(f"   COM gate:  A={gA[0]:.3f}  B={gB[0]:.3f}   {'✓ A preservado' if gA[0]>0.4 else '✗'}")
