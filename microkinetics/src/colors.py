"""Cores das figuras da tese (padrão de 14/10/2026, a partir do Colour Schemes de Paul Tol, SRON/EPS/TN/09-002).

Vermelho e azul são só do produto da troca (MAP e BIS) e das classes GM e GB; o laranja e o
verde-azulado são só a comparação de equivalentes de fenol.
"""

# Produto da troca (conjunto high-contrast); é a mesma cor dos quadrados M e B das tabelas de padrão
MAP = "#BB5566"
BIS = "#004488"
PRODUCT = {"MAP": MAP, "BIS": BIS, "M": MAP, "B": BIS}

# Fundos do p_MAP: incerto (pale) e robusto (medium-contrast)
UNCERTAIN = {"MAP": "#FFCCCC", "BIS": "#BBCCEE"}
ROBUST = {"MAP": "#EE99AA", "BIS": "#6699CC"}

# Classes: GM monoarilóxido (antes G1), GE estequiométrica (antes Z3), GB bisarilóxido (antes G2); trio claro do medium-contrast
CLASS = {"GM": "#EE99AA", "GE": "#EECC66", "GB": "#6699CC"}

# Metátese de Ferreira 2019 (vibrant)
EQUIV1 = "#EE7733"
EQUIV2 = "#009988"
NO_PHENOL = "#BBBBBB"

# EDA e energia estérica (cap. 1)
EDA = {"Esteric": "#BB5566", "Eelect": "#004488", "Eorb": "#6699CC", "dDisp": "#BBCCEE"}

# Imidos (muted, sem os tons de vermelho e azul reservados ao produto)
IMIDO = {"im0": "#BBBBBB", "im1": "#332288", "im2": "#44AA99", "im3": "#AA4499", "im4": "#999933", "im5": "#88CCEE"}

# Família do fenol: com e sem arila em orto (fig:mlrBarrier1All, 17/10). Verde e roxo retomam o sentido do
# diagrama retirado (verde = com arila), em tons do conjunto bright do Tol, fora das cores reservadas ao
# produto, às classes e aos imidos. A regra química está no script (SMARTS, a mesma do make_effect_diagram.py).
PHENOL_FAMILY = {"com arila em orto": "#228833", "sem arila em orto": "#AA3377"}

# Origem das reações da biblioteca do cap. 2 (fig:libraryComposition): conjunto de Ferreira,
# conjunto nosso (as séries desenhadas e os pais) e literatura. Tons do Tol que não gastam o
# vermelho e o azul do produto nem as cores dos imidos.
ORIGIN = {"Ferreira": "#882255", "nosso": "#117733", "literatura": "#DDAA33"}


# --- Perfis de energia livre do cap. 2 (aprovados pelo autor em 14/10/2026) -----------------------------------------
# Cor = ROH (o fenol); estilo da linha = imido; assim a mesma entrada tem a mesma curva em todas as figuras.
# Vermelho e azul (produto) e laranja/verde-azulado (equivalentes) ficam fora das curvas das entradas.
PROFILE_ROH = {"ROH57": "#AA4499", "ROH62": "#44AA99", "ROH70": "#999933", "ROH73": "#332288",
               "ROH28": "#33BBEE", "ROH47": "#117733", "ROH61": "#DDAA33", "ROH41": "#555555",
               "ROH20": "#EE6677"}  # o ROH20 (bifenil-2-ol, entrada 12) entrou em 17/10, com o set4 novo
PROFILE_IMIDO_STYLE = {"adamantila": "-", "2,6-Me2": "--", "2-CF3": ":", "2,6-iPr2": "-.", "C6F5": "-"}
# entrada da tab:reportedReactionData -> (ROH, imido)
PROFILE_ENTRY = {1: ("ROH57", "adamantila"), 2: ("ROH62", "adamantila"), 3: ("ROH70", "adamantila"),
                 4: ("ROH73", "adamantila"), 6: ("ROH57", "2,6-Me2"), 8: ("ROH57", "2-CF3"),
                 # 12 e 14 entraram em 17/10, com o set4 novo (6, 8, 12 e 14)
                 12: ("ROH20", "2,6-iPr2"), 14: ("ROH47", "C6F5"),
                 11: ("ROH70", "2,6-iPr2"), 13: ("ROH28", "C6F5"), 15: ("ROH47", "C6F5"),
                 16: ("ROH61", "C6F5"), 17: ("ROH41", "C6F5")}
# entrada da tab:reportedReactionData -> código antigo da reação (imido-fenol da biblioteca, como em tables/ch2.libraryEnergies.csv)
PROFILE_ENTRY_OLD = {1: "im5-oc1", 2: "im5-oc2", 3: "im5-oc3", 4: "im5-oc4", 6: "im2-oc1", 8: "im4-oc1", 11: "im1-oc3",
                     12: "im1-ad1", 13: "im3-15", 14: "im3-29", 15: "im3-29", 16: "im3-31", 17: "im3-26"}  # a 12 é do set4 novo; a 14 e a 15 são a mesma combinação, com 1 e 2 equivalentes
# comparação de protocolos e de níveis (sistema modelo e entradas 2 e 3): cor = protocolo, tracejado = só a energia xTB
# A entrada 4 (Ibrahem 2009) é o ROH73 (oc4), não o ROH72: o ROH72 é C26H30 e não pertence à série
# octa-hidro (C26H34) das entradas 1 a 4 (correção de 17/10, a partir dos SMILES do índice de ligantes).
PROFILE_DFT = "#EE7733"        # protocolo DFT (ωB97X-D4)
PROFILE_COMPOSITE = "#33BBEE"  # protocolo composto (TPSS0-D3(BJ)//GFN1-xTB)
PROFILE_REFERENCE = "#000000"  # DLPNO-CCSD(T1)/CBS
PROFILE_LEVELS = {"DLPNO-CCSD(T1)/CBS": (PROFILE_REFERENCE, "-"),
                  "ωB97X-D4/def2-TZVP": (PROFILE_DFT, "-"),
                  "TPSS0-D3(BJ)/def2-TZVP SMD (geometria xTB)": (PROFILE_COMPOSITE, "-"),
                  "GFN1-xTB": (PROFILE_COMPOSITE, "--")}
