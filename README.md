# Dados da tese de doutorado

Dados produzidos durante o meu doutorado e que acompanham a tese **"Ferramentas estatísticas
aplicadas no design e predição de processos catalíticos"** (Attilio Chiavegatti Neto, Programa de
Pós-graduação em Química, UFSCar, 2026; orientador Prof. Dr. Marco Antonio Barbosa Ferreira).

O trabalho trata da **troca de ligantes em complexos de molibdênio-alquilideno imido**: o
mecanismo da protonólise pelo fenol (Capítulo 2 da tese) e o escopo computacional de uma
biblioteca de reações de troca (Capítulo 3), com modelos microcinéticos e modelos estatísticos
para prever qual produto se forma, o monoarilóxido ou o bisarilóxido.

Os arquivos são os mesmos do material suplementar da tese. A numeração das tabelas e dos
apêndices citada aqui é a do PDF da tese.

## Como os arquivos estão organizados

O prefixo do nome indica o capítulo de resultados:

- `ch1.` — Capítulo 2, "Estudo mecanístico da protonólise de complexos Mo-alquilidenos" (sistema modelo);
- `ch2.` — Capítulo 3, "Explorando computacionalmente o escopo de troca de ligantes de Mo-alquilidenos" (biblioteca de reações);
- `esi.` — material geral.

Os `.xyz` contêm várias estruturas cada um, em coordenadas cartesianas e em Å. Cada estrutura é
formada pelo número de átomos, por uma linha de comentário com o rótulo da estrutura e pelas
coordenadas. As estruturas aparecem na mesma ordem e com o mesmo rótulo da tabela ou do `.csv`
correspondente.

**Rótulos das espécies do mecanismo** (`ch2.`, e também nos `.csv`): `a1`, aduto da primeira
troca; `ts1`, primeiro estado de transição; `map`, monoarilóxido; `a2`, aduto da segunda troca;
`ts2`, segundo estado de transição; `bis`, bisarilóxido. As espécies comuns a várias reações têm
rótulos próprios: `bispyrrolide_im<n>_<n>`, `roh_<ROH>_<n>` e `pyrrole_dimethylpyrrole_<n>`. Os
códigos de imido (`im0` a `im5`) e de fenol (`ROH<n>`) são os do índice de ligantes da tese
(Apêndice A), com a correspondência com os códigos antigos da biblioteca em
`ch2.phenolCodes.csv`.

## Geometrias (`.xyz`)

| Arquivo | Conteúdo | Tabela na tese |
|---|---|---|
| `ch1.modelMecOmegaTZ.xyz` | 30 estruturas do sistema modelo (mínimos e estados de transição do mecanismo) otimizadas em ωB97X-D4/def2-TZVP CPCM (benzeno) | Tabela B.1, Apêndice B |
| `ch1.modelMecOmegaDZ.xyz` | 18 estruturas do sistema modelo otimizadas em ωB97X-D4/def2-SVP CPCM (benzeno) | Tabela B.2, Apêndice B |
| `ch1.modelMecXTB1.xyz` | 18 estruturas do sistema modelo otimizadas em GFN1-xTB, em vácuo | Tabela B.3, Apêndice B |
| `ch1.modelMecXTB1speDZ.xyz` | 18 geometrias GFN1-xTB: o confôrmero de cada espécie escolhido pelo protocolo composto ωB97X-D4/def2-SVP SMD (benzeno)//GFN1-xTB | Tabela B.4, Apêndice B |
| `ch1.modelMecXTB1speTZ.xyz` | 18 geometrias GFN1-xTB: o confôrmero de cada espécie escolhido pelo protocolo composto ωB97X-D4/def2-TZVP SMD (benzeno)//GFN1-xTB | Tabela B.5, Apêndice B |
| `ch2.RealSystemBestConformers.xyz` | As 108 estruturas do sistema real das reações publicadas, com geometrias e frequências em ωB97X-D4/def2-SVP CPCM (benzeno): os 6 bispirrolidos, o 2,5-dimetilpirrol, os 11 ROH e as 90 espécies do mecanismo das 16 reações. O rótulo é o da própria tabela, sem o número do confôrmero | Tabela D.1, Apêndice D |
| `ch2.libraryEnergies.xyz` | Geometrias GFN1-xTB dos 7033 confôrmeros das 873 espécies da biblioteca de reações | Acompanha `ch2.libraryEnergies.csv` (Apêndice D) |
| `ch2.libraryBestConformers.xyz` | Geometrias GFN1-xTB do confôrmero de menor energia livre de cada uma das 873 espécies da biblioteca | Acompanha `ch2.libraryBestConformers.csv`; são as estruturas das Tabelas D.2 a D.4 |

Observações:

- os arquivos `ch1.modelMecXTB1speDZ.xyz` e `ch1.modelMecXTB1speTZ.xyz` podem repetir estruturas
  de `ch1.modelMecXTB1.xyz`. Os três usam geometrias GFN1-xTB, e cada arquivo traz o confôrmero
  escolhido na sua tabela;
- em `ch2.libraryEnergies.xyz`, o rótulo da linha de comentário é o da coluna `geometry label` do
  `.csv`. Ele difere do rótulo da energia em 163 confôrmeros, que recebem o sufixo `_g<N>` para
  distinguir confôrmeros de projetos de cálculo diferentes que tinham o mesmo rótulo;
- `ch2.RealSystemBestConformers.xyz` é o único arquivo com geometrias de nível DFT; os demais
  trazem geometrias GFN1-xTB.

## Energias (`ch2.libraryEnergies.csv`, `ch2.libraryBestConformers.csv` e `esi.literatureSearch.tsv`)

| Arquivo | Conteúdo | Tabela na tese |
|---|---|---|
| `ch2.libraryEnergies.csv` | Energias dos 7033 confôrmeros das 873 espécies da biblioteca, no protocolo composto TPSS0-D3(BJ)/def2-TZVP SMD (benzeno)//GFN1-xTB. Colunas: índice; `standardized label` (rótulo da espécie e do confôrmero); `E(el)`, a energia eletrônica TPSS0, em E<sub>h</sub>; `G`, a energia livre, igual a E(el) mais a correção térmica GFN1-xTB, em E<sub>h</sub>, a 298.15 K e 1 mol/L; `geometry label` (rótulo da estrutura no `.xyz`) | Apêndice D |
| `ch2.libraryBestConformers.csv` | O confôrmero de menor energia livre de cada uma das 873 espécies. É o confôrmero usado nas barreiras e nas energias de reação do Capítulo 3. Colunas: índice; `standardized label`; `E(el)`; `G` (mesmas unidades e estado padrão do arquivo anterior) | Fonte das Tabelas D.2, D.3 e D.4 |
| `esi.literatureSearch.tsv` | Levantamento bibliográfico de 72 reações de troca de ligantes com complexos Mo-alquilideno (preparo *in situ* ou com isolamento por cristalização). Colunas: entrada; referência (chave BibTeX); composto; imido; ROH; pirrolido; equivalentes de ROH; solvente; temperatura (°C); tempo (min); preparo; produtos BP, MAP e BIS (%); conversão ou rendimento (%); notas | Tabela C.1, Apêndice C |

Abreviações: **BP** = bispirrolido; **MAP** = monoarilóxido pirrolido; **BIS** = bisarilóxido.

## Microcinética (`.csv`)

Resultados da microcinética com propagação de incerteza por Monte Carlo (1000 amostras, semente
42), com energias livres a 1 mol/L. Para cada reação e condição:

- `BP`, `MAP` e `BIS`: medianas das 1000 amostras da quantidade final de bispirrolido,
  monoarilóxido e bisarilóxido, em porcentagem do complexo de partida. As medianas são calculadas
  separadamente e não somam necessariamente 100;
- `p_MAP`: fração das amostras em que o MAP é o produto majoritário.

As quatro condições padrão são:

- **A1**: 25 °C, 6 h, 0.01 mol/L, 1 equivalente de ROH;
- **A2**: 70 °C, 12 h, 0.1 mol/L, 1 equivalente;
- **B1**: 25 °C, 6 h, 0.01 mol/L, 2 equivalentes;
- **B2**: 70 °C, 12 h, 0.1 mol/L, 2 equivalentes.

| Arquivo | Conteúdo | Tabela na tese |
|---|---|---|
| `ch2.microkineticsLibrary.csv` | As 136 reações da biblioteca nas condições A1, A2, B1 e B2 (colunas `<condição> BP`, `<condição> MAP` etc.), com a classe (`class`: GM, GE ou GB; `excluded` nas duas reações retiradas antes da classificação, com o 2,6-di-*terc*-butilfenol) e o padrão do produto majoritário (`pattern`: M quando p_MAP > 0.5 e B caso contrário, nas condições A1, A2, B1 e B2, nessa ordem), com as energias do protocolo composto TPSS0-D3(BJ)/def2-TZVP SMD (benzeno)//GFN1-xTB e incerteza de 3 kcal/mol. Rótulos `im<imido>-<ROH>` como em `ch2.libraryBestConformers.csv` | Base das classes GM, GE e GB e da seção "Biblioteca de reações" do Capítulo 3 |
| `ch2.microkineticsReported.csv` | As 17 reações publicadas calculadas, nas condições de preparo de cada complexo, com as energias DFT ωB97X-D4/def2-TZVP CPCM (benzeno)//ωB97X-D4/def2-SVP CPCM (benzeno) (`DFT`, incerteza de 2 kcal/mol) e do protocolo composto (`composite`, incerteza de 3 kcal/mol). `entry` é a numeração da tabela da tese; o sufixo `-1eq`/`-2eq` separa as reações feitas com quantidades diferentes de ROH | Tabela 3.6 |
| `ch2.microkineticsReportedStandard.csv` | As 15 reações publicadas nas condições padrão A1 a B2, com as energias do protocolo DFT (incerteza de 2 kcal/mol). Nas condições padrão, a quantidade de ROH é fixa, por isso as entradas 8/9 e 14/15 viram uma reação só | Validação das classes (Capítulo 3) |

As classes: **GM**, o monoarilóxido é o produto (o bisarilóxido só se forçar mais); **GE**, o
produto segue a estequiometria (monoarilóxido com 1 equivalente, bisarilóxido com 2); **GB**, o
bisarilóxido é o produto.

Nos rótulos, `ROH<n>` segue o Apêndice A. Os fenóis sintéticos e os ROH das reações publicadas,
que também entram na biblioteca, têm o código antigo (`d<n>m`, `w<n>p`, `oc<n>`, `ad1`, `bitet<n>`)
em `ch2.phenolCodes.csv`. O `oc3` é duplicata do fenol de código `35` (mesmo ROH70) e aparece como
`ROH70b` nos arquivos da biblioteca.

## Modelos estatísticos (`ch2.modelsData.csv`)

As 125 reações da biblioteca usadas nos modelos estatísticos do Capítulo 3 (reações classificadas
com todos os descritores calculados, inclusive as dos fenóis sintéticos). Colunas: `reaction`
(rótulo `im<imido>-ROH<n>`); `ROH` (código do fenol no índice de ligantes, Apêndice A); `imido`;
`phenol_code` (código do fenol na biblioteca); `synthetic_phenol`; `barrier1` e `barrier2`
(barreiras de ativação da primeira e da segunda troca, kcal/mol, protocolo composto, 1 mol/L);
`class` (GM, GE ou GB); `pattern` (produto majoritário em A1, A2, B1 e B2, como em
`ch2.microkineticsLibrary.csv`); `threshold_single`, `threshold_double`, `mlr_barrier2` e
`mlr_barrier1` (`train`, `validation` ou vazio, quando a reação não entra no modelo); e os 322
descritores, com os nomes definidos no Apêndice E. Correspondem às Tabelas E.1 a E.5.

## Cinética

*(em preparação)* — os arquivos da cinética das reações (perfis de energia, agregados das
simulações de Monte Carlo e as saídas por condição) serão acrescentados aqui.

## Métodos

Os cálculos foram feitos com ORCA (DFT e semiempíricos) e CREST/GFN-FF (busca conformacional); a
microcinética e os modelos estatísticos, em Python. Os níveis de teoria de cada arquivo estão nas
tabelas acima, e a descrição completa dos procedimentos está no capítulo de Métodos da tese.

## Código

O código usado nas análises (busca conformacional, microcinética, Monte Carlo e modelos
estatísticos) ainda não está publicado aqui.

## Licença e citação

**MIT** — os dados e o código podem ser usados, modificados e redistribuídos livremente, inclusive
comercialmente, mantendo o aviso de copyright (veja o [LICENSE](LICENSE)).

Se você usar estes dados num trabalho acadêmico, cite a tese:

> CHIAVEGATTI NETO, A. *Ferramentas estatísticas aplicadas no design e predição de processos
> catalíticos*. 2026. Tese (Doutorado em Química) — Universidade Federal de São Carlos, São Carlos,
> 2026.

Este trabalho teve apoio da Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP),
processo 2020/13563-7.
