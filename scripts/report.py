#!/usr/bin/env python3
"""F6 — gerador de report headless. Lê o summary.json de uma bateria (produzido
por analise8h.py) e emite um painel HTML autocontido, sem GUI.

O Bevy NATIVO (o "mundo vivo" que o Leonardo quer ver escoar) precisa da máquina
dele — o sandbox mata processo GUI (exit 144) e ele já rejeitou UIs de browser.
Este report é a alternativa headless que o roadmap lista: visibilidade sem GUI.

uso: python3 scripts/report.py [runs/_battery8h/summary.json] [saida.html]
"""
import json, sys, os, html

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/runs/_battery8h/summary.json"
out = sys.argv[2] if len(sys.argv) > 2 else f"{ROOT}/runs/_battery8h/report.html"
S = json.load(open(src))

def card(tag, eq, ev, verdict, vclass="ok"):
    return f'''<div class="c"><div class="tag">{html.escape(tag)}</div>
      <div class="eq">{eq}</div><div class="ev">{ev}</div>
      <div class="v v-{vclass}">{html.escape(verdict)}</div></div>'''

cards = []
if "lei1" in S:
    l = S["lei1"]
    cards.append(card("Lei 1 · metabolismo crítico", "met<sub>c</sub> = c·σ²",
        f"c = {l['c_mean']:.1f} · CV {l['cv_pct']:.1f}% · n={l['n']}", "confirmada"))
if "lei2" in S:
    l = S["lei2"]
    cards.append(card("Lei 2 · custo por parâmetro (AIC)", f"c = {l['a']:.2f}·dof + {l['b']:.2f}",
        f"R² = {l['r2']:.3f} · n={l['n']}", "confirmada"))
if "denoise" in S:
    d = S["denoise"]
    cards.append(card("Denoise zero-shot", "met = k·σ̂²",
        f"vs MP@oráculo {d['vs_mp'][0]}/{d['vs_mp'][1]} (+{d['d_mp']:.2f} dB) · vs wavelet {d['vs_wavelet'][0]}/{d['vs_wavelet'][1]}", "191+/191+", "ok"))
if "espelho" in S:
    e = S["espelho"]
    cards.append(card("Lei do Espelho", "assinatura ↔ sinal",
        f"corr bruta {e['circular']:.3f} (trivial) → discriminação honesta", "só como assinatura", "op"))
if "adapt" in S:
    a = S["adapt"]
    cards.append(card("Adaptação online", "sinal troca → reescoa",
        f"recupera {100*a['frac']:.0f}% do dano sem reiniciar", "confirmada"))
if "codec" in S:
    c = S["codec"]
    q8 = c.get("q8", {})
    cards.append(card("Codec vs JPEG", "params → bytes",
        f"Q8 {q8.get('bpa',0):.1f} B/átomo — JPEG esmaga 5–18×", "compressão: jogo perdido", "no"))

flow = [
    ("Fundação", "Física agnóstica de modalidade (imagem/áudio/vídeo)", "ok"),
    ("Compressão", "Bater o matching pursuit — jogo perdido (−4 a −27 dB)", "no"),
    ("Gradiente", "Adam nos contínuos, seleção no discreto", "ok"),
    ("A lei", "met_c = σ²(1 + 2·dof) = penalidade do AIC", "ok"),
    ("Denoise", "sem hiperparâmetro, bate MP com N-oráculo", "ok"),
    ("Garimpo", "a lei sai do solver: finanças, estatística, poda de rede", "ok"),
    ("Não-linear", "oclusão — onde o guloso não pode jogar", "ok"),
    ("Geração (F4)", "prior importado → textura nova, sem alvo", "ok"),
]
steps = "".join(f'''<div class="s"><div class="sw">{html.escape(w)}</div>
  <div class="sd">{html.escape(d)} <span class="v v-{c}">{'✓' if c=='ok' else ('×' if c=='no' else '~')}</span></div></div>''' for w, d, c in flow)

doc = f'''<!doctype html><html><head><meta charset="utf-8">
<title>IARA-Água — report</title><style>
:root{{--bg:#f2f5f4;--fg:#14262a;--dim:#52696d;--rule:#c3d0cd;--sur:#fff;--ok:#0f8f7d;--no:#b23b32;--op:#9a6a12;
--serif:"Iowan Old Style",Palatino,Georgia,serif;--mono:ui-monospace,Menlo,monospace}}
@media(prefers-color-scheme:dark){{:root{{--bg:#0d1a1e;--fg:#dfe9e7;--dim:#90a6a8;--rule:#2b444a;--sur:#16292e;--ok:#4fd1b8;--no:#e77c72;--op:#dcae5c}}}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--fg);font-family:var(--serif);line-height:1.5}}
.w{{max-width:60rem;margin:0 auto;padding:2.5rem 1.5rem 4rem}}
h1{{font-size:2rem;font-weight:500;margin:.2rem 0}}
.eyebrow{{font-family:var(--mono);font-size:.72rem;letter-spacing:.12em;text-transform:uppercase;color:var(--dim)}}
.dek{{color:var(--dim);max-width:52ch}}
h2{{font-size:1.3rem;font-weight:500;margin:2.5rem 0 1rem;border-bottom:1px solid var(--rule);padding-bottom:.3rem}}
.flow{{display:grid;gap:0}}.s{{display:grid;grid-template-columns:8rem 1fr;gap:1rem;padding:.6rem 0;border-top:1px solid var(--rule)}}
.sw{{font-family:var(--mono);font-size:.78rem;color:var(--dim)}}.sd{{font-size:.95rem}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(16rem,1fr));gap:1px;background:var(--rule);border:1px solid var(--rule)}}
.c{{background:var(--sur);padding:1.1rem}}.tag{{font-family:var(--mono);font-size:.68rem;letter-spacing:.08em;text-transform:uppercase;color:var(--dim)}}
.eq{{font-style:italic;font-size:1.2rem;padding:.4rem 0;border-bottom:1px solid var(--rule);margin:.4rem 0}}
.ev{{font-family:var(--mono);font-size:.8rem;color:var(--fg)}}
.v{{display:inline-block;font-family:var(--mono);font-size:.66rem;letter-spacing:.06em;text-transform:uppercase;padding:.14rem .45rem;border-radius:2px;margin-top:.5rem}}
.v-ok{{background:color-mix(in srgb,var(--ok) 18%,transparent);color:var(--ok)}}
.v-no{{background:color-mix(in srgb,var(--no) 18%,transparent);color:var(--no)}}
.v-op{{background:color-mix(in srgb,var(--op) 18%,transparent);color:var(--op)}}
footer{{margin-top:3rem;padding-top:1rem;border-top:1px solid var(--rule);color:var(--dim);font-size:.82rem;font-family:var(--mono)}}
</style></head><body><div class="w">
<div class="eyebrow">iara-água · report headless · gerado de {html.escape(os.path.basename(src))}</div>
<h1>O que a bateria provou</h1>
<p class="dek">Report autocontido gerado sem GUI. O observatório Bevy nativo (a água escoando ao vivo) precisa da máquina do Leonardo.</p>
<h2>O caminho</h2><div class="flow">{steps}</div>
<h2>As leis</h2><div class="cards">{"".join(cards)}</div>
<footer>gerado por scripts/report.py · dados em runs/_battery8h/ · leis em docs/</footer>
</div></body></html>'''

open(out, "w").write(doc)
print(f"report → {out} ({len(doc)} bytes)")
print("Bevy nativo (mundo vivo) = máquina do Leonardo; este é o headless.")
