# Aula 1 — Como achar a equação de uma lei (o movimento-mestre)

> Pra Leonardo. Objetivo: você sair daqui sabendo criar uma lei/fórmula sozinho. Uma coisa só, funda.

## A ideia única

Uma **lei** responde a pergunta *"qual é o MELHOR?"*. E a matemática acha o melhor de um jeito só:

> Escreve um **custo** (o quão errado você está) → acha onde ele é **mínimo** → esse ponto é a lei.

E tem um truque pra achar o mínimo: **no fundo do vale, a inclinação é zero.** Inclinação = derivada.
Então a receita inteira é:

**custo → derivada = 0 → resolve → LEI.**

É isso. Metade dos nossos ouros é exatamente esse movimento. Vou te mostrar um nascendo.

## Exemplo trabalhado (do nosso próprio ouro): a amp do ser-água

Pergunta real do solver: *"tenho um resíduo `r` (o que falta explicar) e um átomo com uma forma `ψ`.
Quanto de ψ — chama de `a` — eu adiciono pra melhor cancelar `r`?"*

**Passo 1 — escreve o custo.** "Errado" = o que sobra depois de adicionar `a·ψ`, medido em soma de
quadrados (o tamanho do vetor):

```
C(a) = ‖ r − a·ψ ‖²          (soma de (r_i − a·ψ_i)² em cada pixel i)
```

**Passo 2 — abre o quadrado.** Usando ‖u‖² = ⟨u,u⟩ (produto interno de u consigo):

```
C(a) = ⟨r,r⟩  −  2a·⟨r,ψ⟩  +  a²·⟨ψ,ψ⟩
     =  ‖r‖²  −  2a·⟨r,ψ⟩  +  a²·‖ψ‖²
```

Olha o formato: é `const − (algo)·a + (algo)·a²`. **É uma parábola em `a`** (abre pra cima porque ‖ψ‖²>0).
Toda parábola pra cima tem UM fundo. É lá que mora o melhor `a`.

**Passo 3 — derivada = 0.** Deriva C em relação a `a` (⟨r,r⟩ é constante → some; −2a⟨r,ψ⟩ → −2⟨r,ψ⟩;
a²‖ψ‖² → 2a‖ψ‖²):

```
dC/da = −2⟨r,ψ⟩ + 2a·‖ψ‖² = 0
```

**Passo 4 — resolve pra `a`:**

```
a = ⟨r,ψ⟩ / ‖ψ‖²        ← A LEI  (a nossa "amp em forma-fechada")
```

**A intuição (guarda isso):** `⟨r,ψ⟩` mede o quanto `r` e `ψ` apontam pro mesmo lado (alinhamento);
`‖ψ‖²` é o tamanho do próprio átomo. Então a lei diz: **"bebe o quanto você se alinha, dividido pelo seu
tamanho".** Isso é *projeção* — a sombra de r na direção de ψ. É por isso que o solver não explode: cada
ser bebe EXATAMENTE o ótimo, calculado, não chutado.

Você acabou de ver uma lei nascer de um princípio, em 4 passos.

## O padrão: TODOS os ouros são esse mesmo movimento

- **Metabolismo (AIC)** `met = σ²(1+2·dof)`: mesmo custo dos quadrados, mais uma conta de *"cada parâmetro
  custa σ²"* (σ² = variância do ruído). Um ser só vale a pena se a energia que ele explica passa do que
  ele custa. O ponto de **empate** É a lei.
- **Garimpo (RIC)** `pen = 2σ²·ln(M)`: quando você escolhe o MELHOR de M candidatos, o ruído sozinho às
  vezes parece bom por sorte. O "mais sortudo de M" cresce como `σ·√(2 ln M)` (isso vem de estatística de
  extremos). Elevando ao quadrado → o limiar pra ser REAL é `2σ²·ln M`. Por isso escolher entre muitos
  exige um limiar maior.
- **Água÷grau**: o "melhor recall" quando o sinal flui como um passeio aleatório. O passeio empilha nos
  nós de grau alto (∝ grau) → divide pelo grau pra recuperar o sinal verdadeiro.

Mesmo esqueleto sempre: **define o melhor → otimiza (derivada=0) ou tira a esperança → lei.**

## A SUA receita (pra achar uma lei nova)

1. **Nomeia** a grandeza que você quer que seja a melhor (erro? plausibilidade? quantos seres?).
2. **Escreve** como CUSTO (pra minimizar) ou como PROBABILIDADE (pra maximizar). Soma-de-quadrados é o
   cavalo de batalha.
3. **Acha o ponto plano**: derivada = 0. (Se for "melhor de muitos", usa esperança / extremos.)
4. **Resolve** → sua lei.
5. **VALIDA** nos dados (regra da casa: honestidade brutal). Lei que não prevê é só álgebra bonita.

## Seu primeiro exercício (a gente resolve junto no próximo passo)

Na pirâmide de imagem: *"quando adicionar o grau ℓ+1 PAGA?"* Ele paga se o detalhe que explica bate o que
custa (mais seres). Escreve:

```
ganho(ℓ+1) = redução do resíduo²  −  (nº de seres novos)·σ²
```

Acha o **empate** (`ganho = 0`) → é a lei de "quando parar de descer a pirâmide". É um metabolismo por
grau. Traz sua tentativa que eu corrijo — é assim que você aprende a achar equação: fazendo uma.
