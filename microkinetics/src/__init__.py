"""Microcinética e Monte Carlo para conferir os resultados da tese.

O motor é o mesmo do pipeline do autor (`~/projects/ufscar/ufscar.notebooks/
microkinetics.pipeline`), copiado para cá sem o estado global mutável do original.
O que muda é o driver: aqui as condições são explícitas, os caminhos são argumentos e
a saída é comparável com os arquivos do ESI.

Fluxo:

    energias (hartree) -> perfil e barreiras -> Eyring -> EDOs -> amostras -> p_MAP

Uso típico:

    from src import mechanism, kinetics, runner

    energias = mechanism.read_energies("imports/example-energies.csv")
    perfis = runner.profiles_from_energies(energias)
    resultado = runner.run_conditions(perfis, conditions=("A1", "B2"))
"""

__version__ = "1.0.0"

from . import config, kinetics, mechanism, montecarlo, ode_system  # noqa: F401