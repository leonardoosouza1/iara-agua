#!/usr/bin/env python3
"""POC do modelo byte-hormonal: adaptação a mundo NÃO-ESTACIONÁRIO.

Preditor de bytes online (contexto-1, tabela de logits por byte anterior).
Stream MUDA de distribuição no meio (idioma A → idioma B). Compara:
  - SGD fixo: entropia cruzada, lr constante.
  - HORMONAL: lr gateada por DOPAMINA (surpresa=RPE) + decaimento METABÓLICO
    (esquece contagens velhas) + CORTISOL (sob choque de loss, sobe plasticidade).
Métrica honesta: bpb ao longo do tempo; velocidade de recuperação após a troca.
Sem torch — puro numpy, rápido, determinístico.
"""
import numpy as np, math

def stream(data, n): return np.frombuffer(data[:n], dtype=np.uint8).astype(np.int64)

# dois "mundos": dois corpora reais bem diferentes (PT literatura vs código/instruct)
def load(path, n, off=0):
    with open(path,'rb') as f: f.seek(off); return f.read(n)
A = stream(load("assets/photo_texture.png",0) or b"", 0)  # placeholder
# usa dois textos reais do bytebrain se existirem, senão sintético
import os
BB="../bytebrain/data"
if os.path.exists(f"{BB}/literature_pt.txt") and os.path.exists(f"{BB}/instruct_corpus.txt"):
    A=stream(load(f"{BB}/literature_pt.txt",120000,500000),120000)
    B=stream(load(f"{BB}/instruct_corpus.txt",120000,500000),120000)
else:
    rng=np.random.default_rng(0)
    # mundo A: markov com uma matriz; mundo B: outra bem diferente
    def synth(seed,n):
        r=np.random.default_rng(seed); P=r.dirichlet(np.ones(256)*0.1,256)
        x=[0]
        for _ in range(n): x.append(r.choice(256,p=P[x[-1]]))
        return np.array(x[1:])
    A=synth(1,120000); B=synth(2,120000)

# stream: metade A, metade B (a TROCA no meio)
half=60000
S=np.concatenate([A[:half], B[:half]])

def run(hormonal):
    L=np.zeros((256,256))          # logits[prev, next]
    cnt=np.ones((256,256))*0.1     # pseudo-contagem (a "memória")
    lr=0.5
    bpbs=[]; window=[]; ema_loss=2.0
    for t in range(1,len(S)):
        prev,cur=S[t-1],S[t]
        p=np.exp(L[prev]-L[prev].max()); p/=p.sum()
        loss=-math.log(max(p[cur],1e-9))/math.log(2)
        window.append(loss)
        if len(window)>2000: window.pop(0)
        if t%2000==0: bpbs.append((t,np.mean(window)))
        # atualização
        if hormonal:
            # DOPAMINA: surpresa = loss vs esperado (EMA). lr sobe com surpresa.
            surprise=max(0.0, loss-ema_loss)
            ema_loss=0.999*ema_loss+0.001*loss
            # CORTISOL: choque sustentado (loss média alta) → mais plasticidade
            stress=min(2.0, np.mean(window[-200:])/max(ema_loss,0.5)) if len(window)>=200 else 1.0
            eff_lr=lr*(1.0+2.0*surprise)*stress
            # METABOLISMO: esquece contagens velhas (decai a memória) → adapta
            cnt[prev]*=0.995
        else:
            eff_lr=lr
        # gradiente da cross-entropy no logit (softmax)
        g=p.copy(); g[cur]-=1.0
        L[prev]-=eff_lr*g
        cnt[prev,cur]+=1
    return bpbs

import sys
b_sgd=run(False)
b_hor=run(True)

print("bpb ao longo do stream (troca de distribuição no meio, t=60000):\n")
print(f"{'t':>8} {'SGD fixo':>10} {'hormonal':>10} {'fase':>10}")
# alinha por t
d1={t:v for t,v in b_sgd}; d2={t:v for t,v in b_hor}
for t in sorted(d1):
    if t%8000==0 or (58000<=t<=72000):
        fase="MUNDO A" if t<60000 else "MUNDO B"
        mark=" ← TROCA" if 60000<=t<=62000 else ""
        print(f"{t:>8} {d1[t]:>10.3f} {d2.get(t,float('nan')):>10.3f} {fase:>10}{mark}")

# velocidade de recuperação: bpb médio nos 10k ticks APÓS a troca
post=lambda d: np.mean([v for t,v in d.items() if 60000<t<=70000])
pre=lambda d: np.mean([v for t,v in d.items() if 50000<t<=60000])
print(f"\nbpb médio ANTES da troca (mundo A estável): SGD {pre(d1):.3f} · hormonal {pre(d2):.3f}")
print(f"bpb médio 10k APÓS a troca (adaptar a B):   SGD {post(d1):.3f} · hormonal {post(d2):.3f}")
print(f"\nchoque de adaptação (pós-pré, menor=adapta melhor):")
print(f"  SGD      +{post(d1)-pre(d1):.3f} bpb")
print(f"  hormonal +{post(d2)-pre(d2):.3f} bpb")
w = "HORMONAL adapta melhor" if (post(d2)-pre(d2)) < (post(d1)-pre(d1)) else "SGD adapta melhor (hormonal falhou)"
print(f"  → {w}")
