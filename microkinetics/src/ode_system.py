"""As duas trocas consecutivas de ligantes e a integração das EDOs.

A rede é a do capítulo 2 (Equações `eq:kineticsBP` a `eq:kineticsPyrH`), com lei de
ação das massas e etapas elementares reversíveis:

    bispirrolido + ROH  ⇌  MAP + pirrol        (k1d, k1r)
    MAP + ROH           ⇌  bisarilóxido + pirrol (k2d, k2r)

Os adutos (a1, a2) não entram: as barreiras são medidas das espécies separadas.
"""

import numpy as np
from scipy.integrate import solve_ivp

from . import config


def derivative(t, concentrations, k1d, k1r, k2d, k2r):
    """Derivadas das cinco concentrações, na ordem de `config.ENTITIES`."""
    bispyr, roh, me2pyr, map_, bis = concentrations

    first_forward = k1d * bispyr * roh
    first_reverse = k1r * map_ * me2pyr
    second_forward = k2d * map_ * roh
    second_reverse = k2r * bis * me2pyr

    return [
        -first_forward + first_reverse,
        -first_forward + first_reverse - second_forward + second_reverse,
        first_forward - first_reverse + second_forward - second_reverse,
        first_forward - first_reverse - second_forward + second_reverse,
        second_forward - second_reverse,
    ]


def trajectory_points(total_time, num_points=None, exponent=None):
    """Pontos de tempo, mais densos no início (`t_i = t_final * s_i^expoente`).

    A distribuição cúbica existe para descrever a variação rápida das concentrações
    nos instantes iniciais (Métodos, `sec:methodsMicrokinetics`).
    """
    num_points = config.TRAJECTORY_POINTS_DEFAULT if num_points is None else num_points
    exponent = config.TRAJECTORY_EXPONENT_DEFAULT if exponent is None else exponent
    s = np.linspace(0, 1, num_points)
    return (s ** exponent) * total_time


def integrate(rates, concentrations, times):
    """Integra o sistema com o BDF e as tolerâncias da tese. Devolve o `OdeResult`."""
    initial = [concentrations[entity] for entity in config.ENTITIES]
    arguments = tuple(rates[key] for key in config.RATE_COLUMN_ORDER)

    return solve_ivp(
        fun=derivative,
        t_span=[0, times[-1]],
        y0=initial,
        args=arguments,
        t_eval=times,
        **config.ODE_SOLVER_OPTIONS,
    )


def simulate(rates, concentrations, total_time, num_points=None, exponent=None):
    """Roda uma simulação e devolve as porcentagens finais do complexo de partida.

    As grandezas são as do `analysis.ipynb` do autor: tudo em % do bispirrolido
    inicial, e a conversão é `100 - % de bispirrolido` restante.
    """
    times = trajectory_points(total_time, num_points, exponent)
    solution = integrate(rates, concentrations, times)
    if not solution.success:
        raise RuntimeError(f"a integração falhou: {solution.message}")

    values = dict(zip(config.ENTITIES, solution.y))
    reference = values["bispyr"][0]
    return {
        "bispyr": values["bispyr"][-1] / reference * 100,
        "map": values["map"][-1] / reference * 100,
        "bis": values["bis"][-1] / reference * 100,
        "roh": values["roh"][-1] / values["roh"][0] * 100,
        "conversion": 100 - values["bispyr"][-1] / reference * 100,
    }