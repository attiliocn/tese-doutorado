"""Constante de velocidade pela equação de Eyring, com coeficiente de transmissão 1.

Cópia fiel de `src/kinetics.py` do pipeline do autor (`microkinetics.pipeline`):
muda só a assinatura, que passa a receber as constantes por argumento em vez de as
trazer de um módulo mutável. A expressão e as constantes são as mesmas, e é isso que
faz a rodada bater com os arquivos da tese — não "melhore" a constante dos gases nem
a conversão de unidades sem reconferir contra `check_against_thesis.py`.
"""

import numpy as np

# Constante dos gases em J mol^-1 K^-1, como no pipeline (valor arredondado de
# propósito: é o que gerou os números da tese).
GAS_CONSTANT = 8.314
BOLTZMANN = 1.380649e-23   # J/K
PLANCK = 6.62607015e-34    # J s


def calculate_rate(reaction_barrier, temperature):
    """Barreira (kcal/mol) e temperatura (°C) -> constante de velocidade.

    O resultado sai em s^-1 e é usado sem outra correção como L mol^-1 s^-1 nas
    etapas bimoleculares, o que pressupõe energias no estado padrão de 1 mol/L.
    """
    kelvin = temperature + 273.15
    rt = (GAS_CONSTANT * kelvin) / 1000          # kJ/mol
    prefactor = BOLTZMANN * kelvin / PLANCK      # 1/s
    return prefactor * np.exp(-reaction_barrier / (rt / 4.184))


def rates_from_barriers(barriers, temperature):
    """As quatro barreiras de uma reação -> `{k1d, k1r, k2d, k2r}`."""
    return {
        "k1d": calculate_rate(barriers["b1"], temperature),
        "k1r": calculate_rate(barriers["b3"], temperature),
        "k2d": calculate_rate(barriers["b2"], temperature),
        "k2r": calculate_rate(barriers["b4"], temperature),
    }