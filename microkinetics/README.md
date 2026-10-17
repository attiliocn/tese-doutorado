# Microcinética: conferir os resultados da tese

Programa para rodar a cinética da troca de ligantes com propagação de incerteza por
Monte Carlo e **conferir o resultado contra as tabelas de dados deste repositório**.
O motor é o mesmo do pipeline da tese; o que muda é o que está em volta dele: aqui as
condições são explícitas, os caminhos são argumentos, e a saída é comparável com os
arquivos publicados.

## Como rodar

```bash
cd microkinetics
```

**Requisitos:** Python 3.12, `numpy` 2.2, `pandas` 2.3 e `scipy` 1.16 — as versões em que
o programa foi testado. Não há mais nada além da biblioteca padrão.

**Uma reação, as quatro condições** — é o caminho prático:

```bash
python run_reaction.py im0-ROH28                        # ~10 min nesta máquina
python run_reaction.py im1-ROH3 --condicoes A1 B2       # só duas condições
python run_reaction.py im1-ROH3 --n-samples 200         # conferência rápida
python run_reaction.py --listar                         # reações boas para testar
```

**Um conjunto de reações, pelo notebook:** abra o `run.ipynb`, aponte o
`ENERGIES_FILE` para o seu arquivo de energias e rode. Depois o `validate.ipynb` lê a
pasta gravada e dá o veredito.

O `imports/README.md` descreve o formato do arquivo de energias — inclusive a
pegadinha do estado padrão.

## O que o veredito quer dizer

O ESI é a fonte da verdade. A conferência separa duas perguntas:

| veredito | o que significa |
|---|---|
| `confere` | o produto previsto é o mesmo e os números estão dentro da tolerância |
| `produto confere, números divergem` | o produto majoritário bate, mas as frações não |
| `produto divergente` | a rodada prevê o outro produto |

A tolerância padrão é 0,05 na fração `p_MAP` e 1 ponto percentual nas medianas — mas
**partindo das energias do ESI, com a mesma semente, a diferença tem de ser zero**: se
não for, alguma coisa está errada no caminho. A tolerância existe para quem traz
energias próprias, e aí o que se compara é o **produto previsto** (o lado do corte de
0,5) e a **margem** (a distância a 0,5), não o valor de `p_MAP`.

`p_MAP` é a fração das amostras em que o MAP é o majoritário — o lado do corte de 0,5,
que é a definição de majoritário, não uma convenção. Ela mede o erro do cálculo, não o
acerto experimental: não é probabilidade de acerto nem confiança.

## O teste de aceite

```bash
python check_against_thesis.py      # ~4 s, não roda cinética
```

Confere, contra as tabelas deste repositório e as trajetórias que o acompanham:

1. `tables/ch2.libraryBestConformers.csv` → as tabelas de energias relativas e de barreiras da biblioteca (6 reações);
2. `tables/ch2.RealSystemBestConformers.csv` (energias a 1 atm) → a tabela das reações publicadas calculadas em DFT (8 reações);
3. barreiras → Eyring → as EDOs → as trajetórias sem perturbação de `reference/`, nas quatro condições padrão (16 casos).

Se este script não disser `CONFERE`, nada depois dele vale.

## Custo

Uma integração custa entre 10 e 200 ms, e **as condições forçadas (70 °C) custam uma
ordem de grandeza mais que as brandas**: o sistema fica rígido e o integrador trabalha.
Medido nesta máquina (6 núcleos reais, 12 processos):

| escala | tempo |
|---|---|
| 1 reação × 4 condições × 1000 amostras | ~10 min |
| 3 reações × 2 condições × 300 amostras | ~3 min |
| a biblioteca inteira (136 × 4 × 1000) | dezenas de horas |

É por isso que o caminho é **uma reação por vez**. Amostras cuja integração passa de
60 s são descartadas, como na tese; o resumo de cada condição traz a coluna
`n_timed_out` dizendo quantas.

## Estrutura

```
microkinetics/
├── check_against_thesis.py   o teste de aceite (energias -> tabelas -> simulação)
├── run_reaction.py           uma reação nas quatro condições + veredito (linha de comando)
├── run.ipynb                 um conjunto de reações (notebook)
├── validate.ipynb            lê uma rodada gravada e compara com o ESI
├── src/
│   ├── mechanism.py          energias -> perfil -> barreiras (e a correção de estado padrão)
│   ├── kinetics.py           Eyring
│   ├── ode_system.py         as duas trocas, integradas com BDF
│   ├── montecarlo.py         o esquema de perturbação correlacionada
│   ├── runner.py             o laço reação × condição, os arquivos e o resumo
│   ├── compare.py            a conferência contra as tabelas publicadas
│   ├── inputs.py             monta o arquivo de energias a partir das tabelas
│   ├── colors.py             as cores do produto da troca (MAP e BIS), para as figuras
│   └── config.py             as quatro condições padrão e os parâmetros do Monte Carlo
├── imports/                  as energias do usuário (formato no README de lá)
├── reference/                as trajetórias sem perturbação das 16 simulações da tese,
│                             para o teste de aceite (4 arquivos, ~90 KB)
└── exports/                  as rodadas
```

## O que a rodada grava

Por condição (`exports/<rodada>/<condição>/`):

| arquivo | o que é |
|---|---|
| `aggregate_metrics.csv` | uma linha por reação, com as **35 colunas** do `aggregate_metrics.csv` da tese mais `n_samples` e `n_timed_out` — é o que o `validate.ipynb` compara |
| `<reação>_samples.csv` | as 1000 amostras: o desvio `A`, os `x_ts1`/`x_ts2`, se estourou o tempo, e as grandezas finais |
| `<reação>_baseline.csv` | a trajetória sem perturbação |
| `run_options.json` | tudo o que foi rodado: temperatura, tempo, concentrações, σ, semente, solver |

As trajetórias completas (os 1001 pontos de tempo de cada amostra, ~28 MB por reação
por condição) **não são gravadas por padrão** — o resumo e as amostras bastam para a
conferência. Ligue `SAVE_TRAJECTORIES = True` no `run.ipynb` se precisar delas.

## Procedência

O motor vem do pipeline da tese (fora deste repositório), copiado para cá:
`mechanism.py` do parser de perfis e `kinetics.py`, `ode_system.py` e a sistemática de
Monte Carlo do `src/` daquele pipeline. As cópias são fiéis nas contas — inclusive onde
o original tem um valor arredondado, como `R = 8.314` — porque é isso que faz o
resultado bater com os arquivos da tese. Não "melhore" uma constante sem reconferir com
`check_against_thesis.py`.

As tabelas da pasta `tables/` e as trajetórias de `reference/` acompanham a tese e são
a fonte da verdade da conferência.

## Conferências já feitas

| o que | resultado |
|---|---|
| `check_against_thesis.py` (perfis, barreiras e simulação) | CONFERE |
| `im0-ROH28` nas quatro condições, 1000 amostras (3,3 min) | 4/4 conferem, `p_MAP` e medianas idênticos |
| `im1-ROH20` (publicada) nas condições padrão, 1000 amostras (1,3 min) | 2/2 conferem |

A `im0-ROH28` foi escolhida por ser a reação com a **menor margem ao corte** entre as
136 (0,014 em A1, com `p_MAP` = 0,486) — é o caso em que uma diferença de convenção
apareceria. A saída está em `exports/im0-ROH28/`.
