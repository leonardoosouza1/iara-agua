#!/usr/bin/env python3
"""A LEI DO GARIMPO — a penalidade por parâmetro é 2σ²·ln(M_alternativas).

TESE UNIFICADORA. Nossa lei do metabolismo mede c ≈ 2·ln(M_busca) e c ≈ 1 + 2·dof.
Isso não é sobre Gabors. É o RIC (Risk Inflation Criterion, Foster & George 1994):
quando você garimpa entre p candidatos, o limiar honesto é 2σ²·ln p, não 2σ² (AIC).
O AIC é o caso degenerado M ≈ e — quando você NÃO buscou.

Testes, todos fora do sistema dos seres, com predição registrada antes:

  A) GARIMPO PURO. E[max Sharpe de M estratégias sem alfa] ≈ σ_S·√(2 ln M).
     Predição: a razão → 1 com M grande. (É a raiz do Deflated Sharpe de
     Bailey & López de Prado — e a mesma equação do nosso c.)

  B) SELEÇÃO DE VARIÁVEIS (p >> n). Predição: com penalidade 2σ²·ln p (RIC) o
     suporte verdadeiro é recuperado melhor que com 2σ² (AIC), porque o AIC não
     sabe que você olhou p candidatos. BIC (σ²·ln n) fica no meio.
     Baselines externos obrigatórios: Lasso e OMP.

  C) PODA METABÓLICA DE REDE. Um neurônio "sobrevive" se a perda que ele explica
     supera 2σ²·ln(largura)·(seus parâmetros). Predição: essa poda derivada bate
     poda-por-magnitude no mesmo orçamento de neurônios.
"""
import numpy as np, math, json, sys
from collections import defaultdict

rng_global = np.random.default_rng(0)
OUT = {}

def linfit(x, y):
    a, b = np.polyfit(x, y, 1)
    pred = np.polyval((a, b), x)
    r2 = 1 - np.sum((y-pred)**2)/max(np.sum((y-np.mean(y))**2), 1e-12)
    return a, b, r2

# ══════════════ A) GARIMPO PURO: E[max Sharpe] = σ_S·√(2 ln M) ══════════════
print("="*74)
print("A) LEI DO GARIMPO — E[max Sharpe de M estratégias SEM alfa] ≈ √(2 ln M)")
print("   (Deflated Sharpe / Bailey & López de Prado; mesma raiz do nosso c)")
print("="*74)
print(f"{'M estratégias':>14} {'E[max Sharpe]':>15} {'√(2 ln M)':>12} {'razão':>8}")
rows_a = []
T = 252*2   # 2 anos de retornos diários
for M in [2, 8, 32, 128, 512, 2048, 8192]:
    reps = 400
    best = []
    for _ in range(reps):
        # M estratégias com ZERO alfa: retornos puro ruído
        r = rng_global.normal(0, 1, size=(M, T))
        sr = r.mean(axis=1)/r.std(axis=1)*math.sqrt(T)   # Sharpe anualizado-ish
        best.append(sr.max())
    e_max = float(np.mean(best))
    # σ do Sharpe amostral sob H0 ≈ 1 (por construção, sr ~ N(0,1) aprox)
    pred = math.sqrt(2*math.log(M))
    rows_a.append((M, e_max, pred))
    print(f"{M:>14} {e_max:>15.3f} {pred:>12.3f} {e_max/pred:>8.3f}")
print("\n   O Sharpe que o ACASO produz cresce com √(ln M). Um backtest que testou")
print("   1000 variantes precisa de Sharpe > 3.7 só pra empatar com o acaso.")
OUT["garimpo"] = [{"M": m, "e_max": e, "pred": p} for m, e, p in rows_a]

# ══════════════ B) SELEÇÃO DE VARIÁVEIS: AIC vs BIC vs RIC ══════════════
print("\n"+"="*74)
print("B) SELEÇÃO DE VARIÁVEIS (p >> n) — a penalidade certa é 2σ²·ln p (RIC)?")
print("="*74)

def forward_select(X, y, penalty_per_param, max_k):
    """seleção gulosa: entra a variável que mais reduz o RSS, enquanto
    ΔRSS > penalidade. É exatamente 'sobreviver = explicar mais que o custo'."""
    n, p = X.shape
    chosen, resid = [], y.copy()
    rss = float(resid @ resid)
    for _ in range(max_k):
        best_j, best_gain = -1, 0.0
        for j in range(p):
            if j in chosen: continue
            xj = X[:, j]
            nn = xj @ xj
            if nn < 1e-12: continue
            gain = (xj @ resid)**2 / nn          # redução de RSS ao entrar j
            if gain > best_gain: best_gain, best_j = gain, j
        if best_j < 0 or best_gain <= penalty_per_param:
            break
        chosen.append(best_j)
        Xs = X[:, chosen]
        beta, *_ = np.linalg.lstsq(Xs, y, rcond=None)
        resid = y - Xs @ beta
        rss = float(resid @ resid)
    return chosen

def omp(X, y, k):
    chosen, resid = [], y.copy()
    for _ in range(k):
        corr = np.abs(X.T @ resid)
        corr[chosen] = -1
        j = int(np.argmax(corr)); chosen.append(j)
        Xs = X[:, chosen]
        beta, *_ = np.linalg.lstsq(Xs, y, rcond=None)
        resid = y - Xs @ beta
    return chosen

def lasso_cd(X, y, lam, iters=300):
    n, p = X.shape
    beta = np.zeros(p)
    norms = (X**2).sum(0) + 1e-12
    r = y.copy()
    for _ in range(iters):
        for j in range(p):
            r += X[:, j]*beta[j]
            rho = X[:, j] @ r
            beta[j] = np.sign(rho)*max(abs(rho)-lam, 0)/norms[j]
            r -= X[:, j]*beta[j]
    return [j for j in range(p) if abs(beta[j]) > 1e-8]

def f1(sel, true):
    sel, true = set(sel), set(true)
    if not sel: return 0.0
    tp = len(sel & true)
    prec = tp/len(sel); rec = tp/len(true)
    return 0.0 if prec+rec == 0 else 2*prec*rec/(prec+rec)

n, k_true, sigma = 200, 10, 1.0
res = defaultdict(lambda: defaultdict(list))
for p in [50, 200, 1000]:
    for seed in range(8):
        rng = np.random.default_rng(1000+seed)
        X = rng.normal(size=(n, p))/math.sqrt(n)   # colunas ~norma 1
        true = list(rng.choice(p, k_true, replace=False))
        beta = np.zeros(p); beta[true] = rng.choice([-1, 1], k_true)*rng.uniform(0.6, 1.2, k_true)
        y = X @ beta + rng.normal(0, sigma, n)

        crit = {
            "AIC  (2σ²)":            2*sigma**2,
            "BIC  (σ²·ln n)":        sigma**2*math.log(n),
            "RIC  (2σ²·ln p)":       2*sigma**2*math.log(p),
        }
        for name, pen in crit.items():
            sel = forward_select(X, y, pen, max_k=40)
            res[p][name].append((f1(sel, true), len(sel)))
        # baselines externos
        res[p]["OMP (k=oráculo)"].append((f1(omp(X, y, k_true), true), k_true))
        lam = sigma*math.sqrt(2*math.log(p))       # limiar universal
        res[p]["Lasso (λ universal)"].append((f1(lasso_cd(X, y, lam), true), 0))

print(f"{'p':>6} {'critério':>22} {'F1 do suporte':>15} {'|selec|':>9}")
for p in sorted(res):
    for name in ["AIC  (2σ²)", "BIC  (σ²·ln n)", "RIC  (2σ²·ln p)", "OMP (k=oráculo)", "Lasso (λ universal)"]:
        v = res[p][name]
        f1m = np.mean([a for a, _ in v]); km = np.mean([b for _, b in v])
        star = "  ←" if name.startswith("RIC") else ""
        print(f"{p:>6} {name:>22} {f1m:>15.3f} {km:>9.1f}{star}")
    print()
OUT["selecao"] = {str(p): {n2: float(np.mean([a for a, _ in v])) for n2, v in d.items()} for p, d in res.items()}

print("   Leitura: se o RIC vence o AIC e a vantagem CRESCE com p, a penalidade")
print("   depende do nº de candidatos garimpados — que é exatamente a nossa lei.")

# ══════════════ C) PODA METABÓLICA DE REDE ══════════════
print("\n"+"="*74)
print("C) PODA METABÓLICA — neurônio sobrevive se explica > 2σ²·ln(largura)·dof")
print("="*74)

def train_ridge(H, y, lam=1e-6):
    A = H.T@H + lam*np.eye(H.shape[1])
    return np.linalg.solve(A, H.T@y)

rows_c = []
for seed in range(6):
    rng = np.random.default_rng(300+seed)
    n, d, W = 400, 8, 200
    Xtr = rng.normal(size=(n, d)); Xte = rng.normal(size=(n, d))
    # alvo: soma de 5 features não-lineares (esparso em "neurônios úteis")
    Wt = rng.normal(size=(d, 5))
    f = lambda X: np.tanh(X@Wt) @ np.array([1.5, -1.0, 0.8, -1.3, 0.6])
    sig = 0.5
    ytr = f(Xtr) + rng.normal(0, sig, n); yte = f(Xte)

    # rede aleatória larga (random features), treina a saída
    Wr = rng.normal(size=(d, W)); br = rng.normal(size=W)
    Htr = np.tanh(Xtr@Wr + br); Hte = np.tanh(Xte@Wr + br)
    beta_full = train_ridge(Htr, ytr)
    mse_full = float(np.mean((Hte@beta_full - yte)**2))

    # (1) PODA METABÓLICA: cada neurônio explica (h·r)²/‖h‖²; custo = 2σ²·ln W
    r = ytr.copy()
    custo = 2*sig**2*math.log(W)
    vivos = []
    for _ in range(W):
        best_j, best_g = -1, 0.0
        for j in range(W):
            if j in vivos: continue
            h = Htr[:, j]; nn = h@h
            if nn < 1e-9: continue
            g = (h@r)**2/nn
            if g > best_g: best_g, best_j = g, j
        if best_j < 0 or best_g <= custo: break
        vivos.append(best_j)
        Hs = Htr[:, vivos]; b = train_ridge(Hs, ytr); r = ytr - Hs@b
    k = len(vivos)
    b = train_ridge(Htr[:, vivos], ytr)
    mse_met = float(np.mean((Hte[:, vivos]@b - yte)**2)) if vivos else float("inf")

    # (2) PODA POR MAGNITUDE no MESMO orçamento k (o baseline padrão)
    idx_mag = list(np.argsort(-np.abs(beta_full))[:max(k, 1)])
    b2 = train_ridge(Htr[:, idx_mag], ytr)
    mse_mag = float(np.mean((Hte[:, idx_mag]@b2 - yte)**2))

    # (3) AIC ingênuo (custo 2σ², sem saber que garimpou W neurônios)
    r = ytr.copy(); vivos_aic = []
    custo_aic = 2*sig**2
    for _ in range(W):
        best_j, best_g = -1, 0.0
        for j in range(W):
            if j in vivos_aic: continue
            h = Htr[:, j]; nn = h@h
            if nn < 1e-9: continue
            g = (h@r)**2/nn
            if g > best_g: best_g, best_j = g, j
        if best_j < 0 or best_g <= custo_aic: break
        vivos_aic.append(best_j)
        Hs = Htr[:, vivos_aic]; bb = train_ridge(Hs, ytr); r = ytr - Hs@bb
    b3 = train_ridge(Htr[:, vivos_aic], ytr)
    mse_aic = float(np.mean((Hte[:, vivos_aic]@b3 - yte)**2)) if vivos_aic else float("inf")

    rows_c.append((k, len(vivos_aic), mse_full, mse_met, mse_mag, mse_aic))

k_m = np.mean([r[0] for r in rows_c]); k_a = np.mean([r[1] for r in rows_c])
print(f"   rede larga = 200 neurônios · ruído σ=0.5 · 6 sementes\n")
print(f"{'método':>28} {'neurônios':>10} {'MSE held-out':>14}")
print(f"{'rede cheia (200)':>28} {200:>10} {np.mean([r[2] for r in rows_c]):>14.4f}")
print(f"{'poda METABÓLICA (2σ²·ln W)':>28} {k_m:>10.1f} {np.mean([r[3] for r in rows_c]):>14.4f}  ←")
print(f"{'poda por magnitude (mesmo k)':>28} {k_m:>10.1f} {np.mean([r[4] for r in rows_c]):>14.4f}")
print(f"{'AIC ingênuo (2σ², sem ln W)':>28} {k_a:>10.1f} {np.mean([r[5] for r in rows_c]):>14.4f}")
wins = sum(1 for r in rows_c if r[3] < r[4])
print(f"\n   metabólica < magnitude em {wins}/{len(rows_c)} sementes")
OUT["poda"] = {"k_met": float(k_m), "k_aic": float(k_a),
               "mse_full": float(np.mean([r[2] for r in rows_c])),
               "mse_met": float(np.mean([r[3] for r in rows_c])),
               "mse_mag": float(np.mean([r[4] for r in rows_c])),
               "mse_aic": float(np.mean([r[5] for r in rows_c])),
               "wins_vs_mag": wins, "n": len(rows_c)}

json.dump(OUT, open("runs/_battery8h/lei_garimpo.json", "w"), indent=1)
print("\n→ runs/_battery8h/lei_garimpo.json")
