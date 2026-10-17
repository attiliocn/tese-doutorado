"""Energias -> perfil do mecanismo -> barreiras.

Porta, para fora das células dos notebooks, o que o pipeline do autor fazia em
`exchange-parser/main.ipynb` (célula 16) e `dft-optimized.ipynb` (célula 3) mais a
correção de estado padrão (célula 3 do primeiro, célula 8 do segundo).

Entrada: energias livres **absolutas**, em hartree, uma por espécie:
`bispyrrolide`, `roh`, `pyrrole`, `map`, `bis`, `ts1`, `ts2` e, quando houver,
`a1` e `a2`.

Saída: o perfil relativo (kcal/mol, reagentes isolados = 0) e as quatro barreiras.

Regra da correção de estado padrão (1 atm -> 1 mol/L), como está no parser: **soma-se
RT ln 24.46 = 1.894198 kcal/mol a cada espécie, antes de montar o perfil**. Como o
perfil é montado somando espécies, cada ponto herda a correção uma vez por espécie
somada; o referencial (bispirrolido + 2 ROH) tem três espécies, então:

    a1, ts1, a2, ts2  -> 2 espécies  -> o relativo perde 1 correção  (-1.894198)
    map, bis          -> 3 espécies  -> o relativo não muda

Equivale, em kcal/mol, à regra do parser de perfis do capítulo 2
(`reaction[3:] -= 1.894/627.5096`), que só subtrai nos pontos bimoleculares.

A correção é opcional (`standard_state=True`), para aceitar energias que já venham
convertidas — por exemplo as do `tables/ch2.libraryBestConformers.csv`, que o pipeline
já gravou a 1 mol/L.
"""

from scipy import constants
import numpy as np
import pandas as pd

# Hartree -> kcal/mol, pela mesma expressão do parser (não um literal arredondado).
HARTREE_TO_KCAL_MOL = (
    constants.physical_constants["Hartree energy"][0] / constants.calorie
) * constants.N_A / 1e3

# RT ln 24.46 a 298.15 K: a conversão de 1 atm para 1 mol/L.
STANDARD_STATE_CORRECTION_KCAL = (
    constants.R * 298.15 * np.log(24.46) / constants.calorie / 1e3
)
STANDARD_STATE_CORRECTION_EH = STANDARD_STATE_CORRECTION_KCAL / HARTREE_TO_KCAL_MOL

# Espécies de cada ponto do perfil. O referencial são os reagentes isolados.
PROFILE_COMPOSITION = {
    "reactants": ("bispyrrolide", "roh", "roh"),
    "a1": ("a1", "roh"),
    "ts1": ("ts1", "roh"),
    "map": ("map", "roh", "pyrrole"),
    "a2": ("a2", "pyrrole"),
    "ts2": ("ts2", "pyrrole"),
    "bis": ("bis", "pyrrole", "pyrrole"),
}

# Pontos que o modelo cinético usa. Os adutos (a1, a2) ficam fora: as barreiras são
# medidas das espécies separadas (ver microkinetics-methods-notes.md §1).
REQUIRED_SPECIES = ("bispyrrolide", "roh", "pyrrole", "ts1", "map", "ts2", "bis")
OPTIONAL_SPECIES = ("a1", "a2")

BARRIER_DEFINITION = {
    "b1": ("ts1", "reactants"),   # primeira troca, sentido direto  -> k1d
    "b2": ("ts2", "map"),         # segunda troca, sentido direto    -> k2d
    "b3": ("ts1", "map"),         # primeira troca, sentido inverso  -> k1r
    "b4": ("ts2", "bis"),         # segunda troca, sentido inverso   -> k2r
}


def standard_state_correction_kcal(n_species):
    """Correção líquida que sobra num ponto do perfil, por espécie somada.

    Um ponto com `n_species` espécies contra um referencial de três perde
    `3 - n_species` correções.
    """
    return (3 - n_species) * STANDARD_STATE_CORRECTION_KCAL


def calculate_mechanism(energies, standard_state=True):
    """Perfil relativo do mecanismo, em kcal/mol, com os reagentes isolados em 0.

    `energies` é um mapeamento espécie -> energia livre absoluta (hartree). Espécies
    ausentes do mapeamento ficam ausentes do perfil (`NaN`), como no parser.
    """
    # `None` entra como ausente: o parser usa `np.nan` para o que não existe.
    values = {}
    for key, value in energies.items():
        values[key] = np.nan if value is None else float(value)
    if standard_state:
        correction = STANDARD_STATE_CORRECTION_EH
        values = {k: (v + correction if v == v else v) for k, v in values.items()}

    profile = {}
    for point, composition in PROFILE_COMPOSITION.items():
        profile[point] = sum(values.get(sp, np.nan) for sp in composition)

    profile = pd.Series(profile)
    profile = (profile - profile["reactants"]) * HARTREE_TO_KCAL_MOL
    return profile


def calculate_barriers(profile):
    """As quatro barreiras, em kcal/mol, a partir do perfil relativo."""
    return pd.Series(
        {name: profile[products] - profile[reference]
         for name, (products, reference) in BARRIER_DEFINITION.items()}
    )


def profile_from_energies(energies, standard_state=True):
    """Perfil e barreiras de uma reação. Devolve `(profile, barriers)`."""
    profile = calculate_mechanism(energies, standard_state=standard_state)
    return profile, calculate_barriers(profile)


def check_species(energies, label=""):
    """Confere se o conjunto de espécies dá para simular. Devolve a lista de faltantes."""
    missing = [sp for sp in REQUIRED_SPECIES if sp not in energies or energies[sp] != energies[sp]]
    if missing:
        raise ValueError(
            f"{label or 'reação'}: faltam as espécies obrigatórias {missing}. "
            f"Obrigatórias: {list(REQUIRED_SPECIES)}; opcionais: {list(OPTIONAL_SPECIES)}."
        )
    return [sp for sp in OPTIONAL_SPECIES if sp not in energies or energies[sp] != energies[sp]]


def read_energies(path, standard_state=True, value_column="G_hartree"):
    """Lê o arquivo de energias do usuário: `reaction, species, <value_column>`.

    Devolve `{reaction: {species: energia}}`. O arquivo pode trazer as espécies em
    qualquer ordem e sem os adutos.
    """
    table = pd.read_csv(path)
    for column in ("reaction", "species", value_column):
        if column not in table.columns:
            raise ValueError(
                f"{path}: falta a coluna '{column}'. "
                f"Colunas encontradas: {list(table.columns)}."
            )

    energies = {}
    for reaction, block in table.groupby("reaction", sort=False):
        energies[str(reaction)] = dict(zip(block["species"].astype(str), block[value_column].astype(float)))

    if standard_state:
        # A correção entra aqui, uma vez por espécie, e não no perfil: é o que o
        # parser faz, e é o que permite aceitar arquivos com ou sem a conversão.
        energies = {
            reaction: {sp: value + STANDARD_STATE_CORRECTION_EH for sp, value in block.items()}
            for reaction, block in energies.items()
        }
    return energies