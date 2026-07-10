#!/usr/bin/env python3
"""POC 2 do modelo byte-hormonal — o DESIGN CERTO (multi-seed, rigoroso).

O guia GUIA_CEREBRO_VIVO afirma que o conserto do naive é: gate LOCAL ∈[0,1]
(nunca amplifica, só congela onde já acerta) + esquecer SELETIVO (decai só o
já-aprendido). Aqui MEÇO isso — não aceito a tabela do guia sem medir.

Stream de bytes não-estacionário (idioma A → idioma B no meio), 5 seeds.
5 abordagens no MESMO preditor contexto-1 (L[prev,next]):
  SGD puro · hormonal naive · gate LOCAL · gate+esquecer BLANKET · gate+esquecer SELETIVO
Métricas: estável (bpb no mundo A antes da troca) · choque (subida pós-troca) ·
instab (desvio-padrão do traço de bpb).
"""
import numpy as np, math, os

BB = "../bytebrain/data"
def load(p, n, off):
    with open(p, "rb") as f:
        f.seek(off); return f.read(n)
A = np.frombuffer(load(f"{BB}/literature_pt.txt", 120000, 500000), np.uint8).astype(np.int64)
B = np.frombuffer(load(f"{BB}/instruct_corpus.txt", 120000, 500000), np.uint8).astype(np.int64)
SWITCH = 60000
S = np.concatenate([A[:SWITCH], B[:SWITCH]])

def run(mode, seed):
    rng = np.random.default_rng(seed)
    # embaralha um pouco a posição inicial por seed (robustez)
    off = seed % 500
    s = np.concatenate([A[off:SWITCH+off], B[off:SWITCH+off]])
    L = np.zeros((256, 256))
    lr = 0.5
    ema_ctx = np.ones(256) * 3.0   # loss média por contexto (pra gate local)
    ema_glob = 3.0
    trace = []
    win = []
    for t in range(1, len(s)):
        prev, cur = s[t-1], s[t]
        p = np.exp(L[prev] - L[prev].max()); p /= p.sum()
        loss = -math.log(max(p[cur], 1e-9)) / math.log(2)
        win.append(loss)
        if len(win) > 1500: win.pop(0)
        if t % 1500 == 0: trace.append(np.mean(win))
        # atualização por modo
        g = p.copy(); g[cur] -= 1.0
        if mode == "sgd":
            L[prev] -= lr * g
        elif mode == "naive":
            surprise = max(0.0, loss - ema_glob)
            ema_glob = 0.999*ema_glob + 0.001*loss
            L[prev] -= lr * (1.0 + 2.0*surprise) * g   # amplifica lr GLOBAL (instável)
        else:
            # gate LOCAL ∈[0,1]: plasticidade alta onde erra, congela onde acerta
            gate = np.clip(ema_ctx[prev] / max(ema_glob, 0.3), 0.0, 1.0)
            ema_ctx[prev] = 0.98*ema_ctx[prev] + 0.02*loss
            ema_glob = 0.999*ema_glob + 0.001*loss
            L[prev] -= lr * gate * g
            if mode == "forget_blanket":
                L[prev] *= 0.999                        # esquece TUDO
            elif mode == "forget_seletivo":
                L[prev] *= (1.0 - 0.01*(1.0 - gate))    # esquece só o congelado (gate baixo)
    tr = np.array(trace)
    # índices: estável = últimos antes da troca; choque = primeiros após
    n_pre = SWITCH // 1500
    estavel = tr[max(0,n_pre-7):n_pre].mean()
    pos = tr[n_pre:n_pre+7].mean()
    return estavel, pos - estavel, tr.std()

print("POC 2 — design hormonal certo (5 seeds, stream não-estacionário)\n")
modes = [("SGD puro","sgd"),("hormonal naive","naive"),("gate LOCAL","gate"),
         ("+esquecer BLANKET","forget_blanket"),("+esquecer SELETIVO","forget_seletivo")]
print(f"{'abordagem':>20} {'estável':>9} {'choque':>9} {'instab':>9}")
res={}
for nome, m in modes:
    E=[]; C=[]; I=[]
    for seed in [11,22,33,44,55]:
        e,c,i = run(m, seed); E.append(e); C.append(c); I.append(i)
    res[m]=(np.mean(E),np.mean(C),np.mean(I))
    print(f"{nome:>20} {np.mean(E):>9.3f} {np.mean(C):>+9.3f} {np.mean(I):>9.3f}")

print("\nVEREDITO (vs a tabela do guia):")
sgd=res["sgd"]; naive=res["naive"]; gate=res["gate"]; fb=res["forget_blanket"]; fs=res["forget_seletivo"]
print(f"  naive instável? instab {naive[2]:.2f} vs SGD {sgd[2]:.2f} → {'SIM' if naive[2]>2*sgd[2] else 'não'}")
print(f"  gate LOCAL conserta (empata SGD)? estável {gate[0]:.3f} vs {sgd[0]:.3f}, instab {gate[2]:.3f} vs {sgd[2]:.3f}")
print(f"  esquecer ajuda adaptação? choque BLANKET {fb[1]:+.3f} / SELETIVO {fs[1]:+.3f} vs SGD {sgd[1]:+.3f}")
print(f"  SELETIVO custa menos no estável? {fs[0]:.3f} vs BLANKET {fb[0]:.3f} (SGD {sgd[0]:.3f})")
import json
json.dump({k:list(v) for k,v in res.items()}, open("runs/hormone_poc2.json","w"), indent=1)
