# IARA-Imagem — arquitetura nativa (pirâmide de grafos + água, no latente)

> Princípio: a imagem NÃO nasce de um passo. Ela EMERGE descendo uma pirâmide de grafos, com ÁGUA
> trazendo detalhe grau a grau, tudo na MATRIZ LATENTE, os seres comendo o resíduo. Os "passos" do
> modelo de imagem (denoising grosso→fino) viram **água descendo o grafo** — o análogo IARA-nativo.

## Fluxo
1. **imagem → VAE.encode → latente Z `[4,64,64]`** (a matriz). O VAE do SD = córtex herdado latente↔pixel.
2. (opcional) **Z → bytes (int8)**: o grafo opera em byte (tese-substrato; float→byte preserva estrutura).
3. **Pirâmide**: Z em L níveis grosso→fino (Z₀ subamostrado … Z_L cheio). Cada nível = um **grafo** (células latentes = nós).
4. **PASSOS = descer o grau:**
   - **grau 0**: seres comem o resíduo de Z₀ (estrutura global). Metabolismo decide a população. Amp forma-fechada. Gauss-Seidel harmoniza.
   - **água÷grau**: propaga nível ℓ→ℓ+1 (pai→filho). O grosso semeia o fino; a água escorre trazendo detalhe; o resíduo na escala fina é a comida dos seres finos.
   - repete até a resolução cheia.
5. **Cara a cara**: em cada nível, resíduo = latente-alvo − latente-gerado = a comida. Estabilizar = latentes iguais = a água virou a foto.
6. **latente-gerado → VAE.decode → imagem.**
7. **Verificador (Morozov)**: sabe quando um grau fechou (resíduo ~ piso de ruído) → quando descer / parar.

## Onde entra cada OURO
| ouro | papel aqui |
|---|---|
| solver-água (`Field` N-dim agnóstico) | um latente É um sinal → os seres o devoram; a pirâmide emerge |
| água÷grau | propagação entre graus (pai→filho), o multi-hop |
| metabolismo (AIC) | quantos seres por grau |
| garimpo (RIC) | podar ser/grau que não paga |
| amp forma-fechada | cada ser bebe o ótimo (não explode) |
| bytes | latente quantizado, grafo em byte |
| SD (VAE) | pesado que semeia (latente↔pixel + talvez o grau-0); o refino é da IARA |
| verificador | sei/não-sei quando o grau fechou |

## Leis a ACHAR (matemática — ver AULA_MATEMATICA_LEIS.md)
- **Quantos seres no grau ℓ?** (metabolismo por nível)
- **Quando o grau ℓ+1 PAGA?** (break-even: detalhe explicado vs custo dos novos seres)
- **Quanto de água passa pai→filho?** (condução ótima entre graus)
