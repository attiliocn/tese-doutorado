"""Constantes do sistema microcinético.

Tudo o que antes era variável global mutável (`configure_standard_condition`
reescrevia o módulo em tempo de execução) está aqui como dado: as condições são
dicionários, e quem roda passa a condição explicitamente. Sem estado global, a
rodada fica reproduzível e o mesmo processo pode rodar várias condições em sequência.
"""

# Espécies do sistema de EDOs, na ordem usada pelo integrador.
ENTITIES = ["bispyr", "roh", "me2pyr", "map", "bis"]

# Quem recebe o desvio sorteado. Os adutos (a1, a2) não entram no modelo cinético,
# mas ficam na lista porque é o vocabulário do mecanismo.
INTERMEDIATES = ["map", "bis"]
TRANSITION_STATES = ["ts1", "ts2"]

# As quatro condições padrão da tese (Métodos, sec:methodsMonteCarlo).
# `sigma` não é usado por `configure_standard_condition` do pipeline antigo: o desvio
# depende do nível de teoria (3 kcal/mol para o protocolo composto, 2 para DFT), não
# da condição. Está aqui só como referência dos valores usados na tese.
STANDARD_CONDITIONS = {
    "A1": {"temperature": 25, "simulation_time": 3600 * 6,
           "bispyr": 0.01, "roh": 0.01, "sigma": 3.0},
    "B1": {"temperature": 25, "simulation_time": 3600 * 6,
           "bispyr": 0.01, "roh": 0.02, "sigma": 3.0},
    "A2": {"temperature": 70, "simulation_time": 3600 * 12,
           "bispyr": 0.1, "roh": 0.1, "sigma": 3.0},
    "B2": {"temperature": 70, "simulation_time": 3600 * 12,
           "bispyr": 0.1, "roh": 0.2, "sigma": 3.0},
}

# Concentrações iniciais: o complexo de partida e o fenol, com os produtos em zero.
def initial_concentrations(condition):
    """Concentrações iniciais de uma condição padrão, no vocabulário das EDOs."""
    settings = STANDARD_CONDITIONS[condition]
    return {"bispyr": settings["bispyr"], "roh": settings["roh"],
            "me2pyr": 0.0, "map": 0.0, "bis": 0.0}


# Integrador: sistema rígido (as constantes diferem em ordens de grandeza).
ODE_SOLVER_OPTIONS = {"method": "BDF", "rtol": 1e-12, "atol": 1e-12}

# Pontos de tempo, mais densos no início (t_i = t_final * s_i^3).
TRAJECTORY_POINTS_DEFAULT = 1000
TRAJECTORY_EXPONENT_DEFAULT = 3

# Monte Carlo: os valores das rodadas da tese.
PERTURBATION_SIGMA_DEFAULT = 3.0      # protocolo composto; 2.0 para DFT
PERTURBATION_SEED_DEFAULT = 42
PERTURBATION_N_SAMPLES_DEFAULT = 1000

# Limite de tempo por amostra, em segundos. O pipeline do autor tinha 60 na
# assinatura e 300 na docstring; aqui é uma constante única. É parte do desenho:
# em B2, dezenas de reações da biblioteca estouram o limite.
JOB_TIMEOUT_SECONDS = 60

# Amostras enviadas por vez a cada processo de trabalho. Cada amostra devolve um
# dicionário pequeno, e o custo de despachar uma a uma é comparável ao da integração:
# lotes de 25 cortam o vaivém entre processos sem atrasar o balanceamento.
POOL_CHUNKSIZE = 25

# Mapeamento barreira -> constante de velocidade e a ordem esperada pelas EDOs.
BARRIER_COLUMN_MAPPING = {"b1": "k1d", "b2": "k2d", "b3": "k1r", "b4": "k2r"}
RATE_COLUMN_ORDER = ["k1d", "k1r", "k2d", "k2r"]