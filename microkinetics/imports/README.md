# Formato do arquivo de energias

O `run.ipynb` lê um CSV com **uma linha por espécie**. As energias são energias livres
absolutas, em hartree, do confôrmero escolhido de cada espécie.

```csv
reaction,species,G_hartree
im1-ROH3,bispyrrolide,-1555.20841600
im1-ROH3,roh,-307.50666930
im1-ROH3,pyrrole,-288.81575970
im1-ROH3,a1,-1862.70962900
im1-ROH3,ts1,-1862.69869300
im1-ROH3,map,-1573.91440700
im1-ROH3,a2,-1881.41953800
im1-ROH3,ts2,-1881.40639000
im1-ROH3,bis,-1592.61934400
```

## Colunas

| coluna | o que é |
|---|---|
| `reaction` | rótulo da reação, no formato `im<imido>-<fenol>` — por exemplo `im1-ROH3` |
| `species` | a espécie (ver abaixo) |
| `G_hartree` | energia livre **absoluta**, em hartree |

O nome da coluna de energia pode ser outro: passe `value_column="..."` para
`mechanism.read_energies`. As linhas podem vir em qualquer ordem, e as espécies
opcionais podem faltar.

## Espécies

| espécie | obrigatória? | o que é |
|---|---|---|
| `bispyrrolide` | sim | o complexo de partida |
| `roh` | sim | o fenol |
| `pyrrole` | sim | o pirrol liberado (2,5-dimetilpirrol) |
| `ts1` | sim | estado de transição da primeira troca |
| `map` | sim | monoarilóxido |
| `ts2` | sim | estado de transição da segunda troca |
| `bis` | sim | bisarilóxido |
| `a1` | não | aduto antes do TS1 |
| `a2` | não | aduto entre o MAP e o TS2 |

**Os adutos são opcionais de propósito.** Eles não entram no modelo cinético: as
barreiras são medidas das espécies separadas (reagentes → TS1, e MAP + fenol → TS2),
porque os adutos são endergônicos no capítulo 1. Se você os tiver, inclua — eles
aparecem no perfil, mas não mudam a cinética.

Se faltar alguma espécie obrigatória, a reação é pulada com uma mensagem dizendo qual
faltou, e a rodada continua com as outras.

## Estado padrão: a pegadinha

O programa precisa saber se as suas energias já estão no estado padrão de **1 mol/L**,
que é o que a equação de Eyring usada aqui pressupõe. Isso é a célula 1 do
`run.ipynb`:

| `ENERGIES_AT_STANDARD_STATE` | quando usar |
|---|---|
| `False` | as energias **já estão** a 1 mol/L. É o caso de quem parte do `tables/ch2.libraryBestConformers.csv`, que o pipeline do autor já converteu |
| `True` | as energias estão como saíram do cálculo, a 1 atm. O programa soma 1,894198 kcal/mol **a cada espécie** antes de montar o perfil |

A conversão é a do parser da tese: soma-se RT ln 24.46 à energia de cada espécie, e só
depois o perfil é montado. Como o perfil é uma soma de espécies contra um referencial
de três, cada ponto fica com a correção uma vez por espécie somada — os intermediários
e os estados de transição (`a1`, `ts1`, `a2`, `ts2`, duas espécies cada) ficam
1,894 kcal/mol abaixo do que ficariam se a correção fosse somada ao ponto, e `map` e
`bis` (três espécies cada, como o referencial) não mudam.

**É esta convenção que faz o resultado bater com a tese.** Verifique com
`python check_against_thesis.py` antes de confiar numa rodada.

## De onde tirar as energias

**Biblioteca do capítulo 3 (136 reações).** O `tables/ch2.libraryBestConformers.csv` já
tem o confôrmero de menor energia livre de cada espécie, no estado padrão de 1 mol/L.
Para montar o arquivo:

```python
from src import inputs
energies = inputs.library_energies("../tables", reactions=["im1-ROH3", "im2-ROH4"])
inputs.write_energies_file(energies, "imports/minhas-energias.csv")
```

**Reações publicadas (16 em DFT).** As energias absolutas estão na tabela do sistema
real do Apêndice D (`tables/appendix/RealSystemBestConformers.tex`, coluna `G`, em
hartree, a 1 atm) — então `ENERGIES_AT_STANDARD_STATE = True`. Os rótulos das espécies
nessa tabela seguem o código antigo do fenol (`im1-ROH20_a1`), que é o mesmo rótulo da
reação.

**Suas próprias energias.** Mesmo formato, e decida o estado padrão conforme a tabela
acima.
