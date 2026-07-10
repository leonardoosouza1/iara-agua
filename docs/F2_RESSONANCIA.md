# F2 — Órgão temporal CfC + a ressonância como equação

## Paridade do CfC (exata)

Portei `CfCCell` — a fórmula do Universe, idêntica:
`state ← state·exp(−decay·dt) + input·gain·dt`, saída `tanh(state)`.
Teste de paridade contra a fórmula manual: **max|diff| = 0.0 (exata)**. O
ecossistema agora tem o CfC em Rust dentro do agua-lab.

## A ressonância como equação (previ errado, medi certo)

Um CfC rastreando um sinal senoidal ruidoso: qual `decay` minimiza o erro?

| freq | decay* | período |
|---|---|---|
| 0.004 | 0.185 | 250 |
| 0.016 | 0.319 | 62.5 |
| 0.064 | 0.552 | 15.6 |
| 0.128 | 0.795 | 7.8 |

**Ajuste: decay* = 1.78·f^0.406, R² = 0.992.**

⚠️ **Previ α ≈ 1 (linear), deu α ≈ 0.4.** A ressonância é uma **lei de potência
sublinear**, não linear. Monótona — confirma qualitativamente as baterias 10/12 do
Universe (o decay ótimo cresce com a frequência) — mas o expoente ~0.4 é novo e
não trivial.

Não é bias-variância simples: o expoente fica estável (0.39–0.43) para σ ∈ {0.1,
0.3, 0.6} e só degenera em σ = 1.0 (ruído domina, medida perde sentido). Então a
constante de tempo do filtro NÃO casa o período linearmente — casa a raiz ~0.4.
A explicação exata fica em aberto (provavelmente ligada ao ganho do integrador
leaky num sinal oscilante, que atenua ∝ 1/√(decay²+ω²)).

## Veredito de F2

- ✅ CfC portado com paridade exata.
- ✅ Ressonância confirmada como equação medida (não só monotonia): potência 0.4.
- ⚠️ Minha predição de linearidade foi refutada — registrado.
- A pergunta da harmonia (visão §7) já foi respondida na Lei 6: a harmonia NÃO
  emerge sozinha, precisa de predação (transição de fase θ*≈0.5). F2 não repete.
