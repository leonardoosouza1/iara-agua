#!/usr/bin/env python3
"""Fila COMPLEMENTAR da bateria 8h — o pesado e o diverso (roda em paralelo).
Resultados no MESMO results.jsonl (append, 1 linha JSON por validação)."""
import json, os, subprocess, time, itertools
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = f"{ROOT}/target/release/agua-lab"
OUT = f"{ROOT}/runs/_battery8h"
RESULTS = open(f"{OUT}/results.jsonl", "a", buffering=1)
LOG = open(f"{OUT}/progress_extra.log", "a", buffering=1)
DEADLINE_H = float(os.environ.get("BATTERY_HOURS", "6.8"))
T0 = time.time()
def left(): return DEADLINE_H*3600 - (time.time()-T0)
def log(m):
    stamp = f"[extra {(time.time()-T0)/3600:5.2f}h]"
    print(stamp, m); LOG.write(f"{stamp} {m}\n")

def emit(group, row, wall):
    row = dict(row); row["_group"] = group; row["_wall_s"] = round(wall, 2)
    row["_extra"] = True
    RESULTS.write(json.dumps(row) + "\n")

def probe(group, args):
    t = time.time()
    try:
        r = subprocess.run([BIN, "probe"] + [str(a) for a in args],
                           capture_output=True, text=True, timeout=7200)
        line = [l for l in r.stdout.strip().split("\n") if l.startswith("{")]
        if not line:
            emit(group, {"error": r.stderr[-300:] or "sem saida", "args": [str(a) for a in args]}, time.time()-t); return
        emit(group, json.loads(line[-1]), time.time()-t)
    except Exception as e:
        emit(group, {"error": str(e)[:300]}, time.time()-t)

def denoise_cmd(group, args):
    t = time.time()
    try:
        r = subprocess.run([BIN, "denoise"] + [str(a) for a in args],
                           capture_output=True, text=True, timeout=7200)
        line = [l for l in r.stdout.strip().split("\n") if l.startswith("[{")]
        if not line:
            emit(group, {"error": r.stderr[-300:], "args": [str(a) for a in args]}, time.time()-t); return
        rows = json.loads(line[-1]); d = {x["metodo"]: x for x in rows}
        emit(group, {"kind": "denoise", "args": [str(a) for a in args],
                     "gauss_or": d.get("gauss@oráculo", {}).get("psnr_limpo"),
                     "vida": d.get("SOBREVIVÊNCIA", {}).get("psnr_limpo"),
                     "vida_n": d.get("SOBREVIVÊNCIA", {}).get("atomos"),
                     "mp_or": d.get("MP@N-oráculo", {}).get("psnr_limpo")}, time.time()-t)
    except Exception as e:
        emit(group, {"error": str(e)[:300]}, time.time()-t)

Q = []
S3 = [11, 22, 33]

# denoise em ESCALA (a lei diz que met_c independe da resolução — teste aplicado!)
for size, sd in itertools.product([128, 192, 256], S3):
    Q.append(("denoise_escala", lambda s=size, sd=sd: denoise_cmd("denoise_escala",
        ["--target", "assets/photo_texture.png", "--size", s, "--sigma", 0.15,
         "--seed", sd, "--ticks", 1000, "--auto"])))

# áudio REAL (fala e trompete) — denoise + morozov
for alvo, sg, sd in itertools.product(["assets/speech.wav", "assets/trumpet.wav"], [0.1, 0.2], S3):
    Q.append(("denoise_audio_real", lambda a=alvo, sg=sg, sd=sd: denoise_cmd("denoise_audio_real",
        ["--modality", "audio", "--target", a, "--size", 8192, "--sigma", sg,
         "--seed", sd, "--ticks", 800, "--auto", "--k", 25])))

# γ(SNR) fino nas 2 modalidades (grade densa de k)
for k, sg in itertools.product([12, 16, 20, 25, 30, 35, 45, 55, 70], [0.08, 0.12, 0.18, 0.25]):
    Q.append(("gamma_fino_img", lambda k=k, sg=sg: probe("gamma_fino_img",
        ["--kind", "morozov", "--k", k, "--sigma", sg, "--seed", 11, "--ticks", 800,
         "--target", "assets/photo_texture.png"])))

# deriva LONGA (memória/estabilidade): 100k ticks
for sd in [11, 22]:
    Q.append(("drift_100k", lambda sd=sd: probe("drift_100k",
        ["--kind", "drift", "--ticks", 100000, "--seed", sd, "--met", 0.3])))

# Lei da População em VÁRIAS resoluções (N = E/met deve valer em qualquer tamanho)
for size, met, sd in itertools.product([48, 96, 128], [0.1, 0.3], S3):
    Q.append(("lei_populacao_escala", lambda s=size, m=met, sd=sd: probe("lei_populacao_escala",
        ["--kind", "equilibrium", "--size", s, "--met", m, "--seed", sd, "--ticks", 1500])))

# Espelho espectral com MAIS alvos (fotos reais em 2 tamanhos)
for alvo, size, sd in itertools.product(
        ["assets/photo_texture.png", "assets/photo_smooth.png"], [96, 128], S3):
    Q.append(("lei_espelho_real", lambda a=alvo, s=size, sd=sd: probe("lei_espelho_real",
        ["--kind", "spectrum", "--target", a, "--size", s, "--met", 0.2, "--seed", sd, "--ticks", 1500])))

# adaptação com magnitudes de mudança diferentes (A→B sintético jitter? usa fotos)
for met, sd in itertools.product([0.3, 0.8], S3):
    Q.append(("adapt_extra", lambda m=met, sd=sd: probe("adapt_extra",
        ["--kind", "adapt", "--met", m, "--seed", sd, "--ticks-a", 2500, "--ticks-b", 1000])))

# codec em populações grandes (met baixo = milhares de átomos)
for met, sd in itertools.product([0.03, 0.1], S3):
    Q.append(("codec_grande", lambda m=met, sd=sd: probe("codec_grande",
        ["--kind", "codec", "--target", "assets/photo_texture.png", "--met", m,
         "--seed", sd, "--ticks", 2000])))

# perf com tudo grande (o teto da máquina)
for size in [192, 256, 384]:
    Q.append(("perf_teto", lambda s=size: probe("perf_teto",
        ["--kind", "perf", "--size", s, "--ticks", 800, "--cap", 6000, "--met", 0.05, "--seed", 11])))

log(f"fila extra: {len(Q)} validações pesadas · deadline {DEADLINE_H}h · workers 4")
done = 0
with ThreadPoolExecutor(max_workers=4) as ex:
    futs = {ex.submit(f): g for g, f in Q if left() > 0}
    for fut in as_completed(futs):
        done += 1
        if done % 5 == 0:
            log(f"{done}/{len(Q)} · restam {left()/3600:.2f}h")
log(f"FIM extra: {done} validações em {(time.time()-T0)/3600:.2f}h")
