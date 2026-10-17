"""Confere o motor contra os números da tese, sem rodar a cinética.

Serve de teste de aceite da pasta: se os perfis e as barreiras que saem das energias
do ESI não forem os das tabelas do capítulo, nada depois disso vale.

Três conferências:
1. `tables/ch2.libraryBestConformers.csv` (energias a 1 mol/L) -> `tab:mechanism` e `tab:barriers`;
2. `tables/appendix/RealSystemBestConformers.tex` (energias a 1 atm) -> `tab:mechanismsReactionsDFT`;
3. barreiras -> Eyring -> EDOs -> o `unperturbated.csv` das rodadas do autor,
   nas quatro condições padrão (o teste forte do motor, porque passa por tudo).

Uso: python check_against_thesis.py
"""

from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from src import mechanism as mech  # noqa: E402
from src import config, inputs, kinetics, ode_system  # noqa: E402

ROOT = HERE.parents[0]
ESI = ROOT / "tables"   # as tabelas do material suplementar
TABLES = ROOT / "tables"
RUNS = HERE / "reference"       # as referências do autor para o teste de aceite

# tabela do apêndice D: perfil relativo (kcal/mol) de seis reações da biblioteca
MECHANISM_TABLE = {
    "im1-ROH3": (3.4, 10.3, -9.5, -8.5, -0.2, -18.3),
    "im1-ROH4": (4.4, 10.8, -8.3, -6.1, 1.8, -19.0),
    "im1-ROH5": (6.3, 15.0, -8.2, -5.9, 4.5, -20.1),
    "im1-ROH8": (7.0, 13.4, -7.8, -7.1, 3.5, -19.2),
    "im1-ROH9": (3.3, 9.8, -9.4, -9.1, -0.1, -19.2),
    "im1-ROH14": (6.5, 14.3, -9.5, -7.6, -0.8, -21.3),
}
# tabela do apêndice D: as quatro barreiras das mesmas reações
BARRIERS_TABLE = {
    "im1-ROH3": (10.3, 9.2, 19.8, 18.0),
    "im1-ROH4": (10.8, 10.0, 19.1, 20.7),
    "im1-ROH5": (15.0, 12.7, 23.1, 24.6),
    "im1-ROH8": (13.4, 11.4, 21.3, 22.8),
    "im1-ROH9": (9.8, 9.3, 19.2, 19.1),
    "im1-ROH14": (14.3, 8.7, 23.8, 20.5),
}
PROFILE_POINTS = ("a1", "ts1", "map", "a2", "ts2", "bis")


def load_best_conformers():
    """Energias livres (hartree, 1 mol/L) do confôrmero de menor G de cada espécie."""
    table = pd.read_csv(ESI / "ch2.libraryBestConformers.csv")
    table.columns = [c.strip() for c in table.columns]
    return dict(zip(table["standardized label"], table["G"]))


BEST = load_best_conformers()


def _find(energies, prefix):
    """Energia da espécie `prefix`, aceitando o rótulo exato ou `prefix_<índice>`.

    Os dois arquivos de origem rotulam de formas diferentes: a tabela do sistema real
    usa o rótulo exato (`im1-ROH20_map`, `dimethylpyrrole`), e o arquivo do ESI traz o
    índice do confôrmero (`im0-ROH3_map_33`, `pyrrole_dimethylpyrrole_2`). O sufixo
    `_` no segundo caso é obrigatório: sem ele, `roh_ROH3` casaria com `roh_ROH30_1`.
    """
    candidates = [v for k, v in energies.items()
                  if k == prefix or k.startswith(prefix + "_")]
    return min(candidates) if candidates else None


def library_energies(imido, roh):
    """Energias da biblioteca para uma reação, a partir do arquivo do ESI."""
    return {
        "bispyrrolide": _find(BEST, f"bispyrrolide_{imido}"),
        "roh": _find(BEST, f"roh_{roh}"),
        "pyrrole": _find(BEST, "pyrrole_dimethylpyrrole"),
        **{s: _find(BEST, f"{imido}-{roh}_{s}") for s in PROFILE_POINTS},
    }


def check_library():
    print("=== 1. biblioteca: tables/ch2.libraryBestConformers.csv -> tabelas do apêndice D ===")
    ok = 0
    for reaction, reference in MECHANISM_TABLE.items():
        imido, roh = reaction.split("-")
        profile = mech.calculate_mechanism(library_energies(imido, roh), standard_state=False)
        barriers = mech.calculate_barriers(profile)
        got = tuple(round(float(profile[p]), 1) for p in PROFILE_POINTS)
        got_barriers = tuple(round(float(barriers[b]), 1) for b in ("b1", "b2", "b3", "b4"))
        same = got == reference and got_barriers == BARRIERS_TABLE[reaction]
        ok += same
        print(f"  {'ok  ' if same else 'DIFF'} {reaction:<11} perfil={got}")
        print(f"       {'':<11} esperado={reference}")
        if got_barriers != BARRIERS_TABLE[reaction]:
            print(f"       barreiras={got_barriers} esperado={BARRIERS_TABLE[reaction]}")
    print(f"  -> {ok}/{len(MECHANISM_TABLE)} reações com perfil e barreiras exatos\n")
    return ok == len(MECHANISM_TABLE)


def check_dft():
    print("=== 2. publicadas: tables/ch2.RealSystemBestConformers.csv (G a 1 atm) "
          "-> tab:mechanismsReactionsDFT ===")
    energies = inputs.published_energies(ESI)
    published = {  # rótulo -> perfil esperado (a1, ts1, map, a2, ts2, bis)
        "im1-ROH20": (2.3, 12.4, -17.0, -16.8, -5.5, -29.7),
        "im1-ROH57": (4.4, 11.6, -16.1, None, 0.0, -26.1),
        "im1-ROH70": (8.3, 14.6, -14.6, None, 8.9, -21.8),
        "im2-ROH49": (2.9, 12.0, -13.8, -5.7, 1.7, -24.9),
        "im2-ROH57": (3.6, 15.2, -16.1, -15.6, -2.9, -26.5),
        "im3-ROH28": (9.3, 16.9, -13.8, -10.8, 2.8, -30.8),
        "im3-ROH41": (12.7, 20.5, -15.9, None, 21.1, -20.4),
        "im3-ROH47": (9.6, 17.9, -15.5, None, 4.4, -29.7),
    }
    ok = 0
    for reaction, reference in published.items():
        profile = mech.calculate_mechanism(energies[reaction], standard_state=True)
        got = tuple(None if profile[p] != profile[p] else round(float(profile[p]), 1)
                    for p in PROFILE_POINTS)
        same = all(b is None or abs(a - b) < 0.2 for a, b in zip(got, list(reference)))
        ok += same
        print(f"  {'ok  ' if same else 'DIFF'} {reaction:<11} calculado={got}")
        print(f"       {'':<11} esperado ={reference}")
    print(f"  -> {ok}/{len(published)} reações conferem; "
          f"{len(energies)} perfis publicados no total\n")
    return ok == len(published)


def check_simulation():
    print("=== 3. barreiras -> Eyring -> EDOs -> unperturbated.csv do autor ===")
    if not RUNS.exists():
        print("  pulado: as referências do autor não estão em microkinetics/reference/\n")
        return True
    rows = pd.read_csv(ESI / "ch2.microkineticsLibrary.csv").set_index("reaction")
    codes = pd.read_csv(ESI / "ch2.phenolCodes.csv")
    lib2roh = dict(zip(codes["library code"].astype(str), codes["ROH code in ESI labels"]))
    ok = total = 0

    for condition, folder in [("A1", "lib_a1"), ("A2", "lib_a2"),
                              ("B1", "lib_b1"), ("B2", "lib_b2")]:
        path = RUNS / f"{folder}-unperturbated.csv"
        if not path.exists():
            continue
        table = pd.read_csv(path)
        settings = config.STANDARD_CONDITIONS[condition]
        for _, row in table.head(4).iterrows():
            label = row["trajectory_alias"]
            imido, code = label.split("-")
            reaction = f"{imido}-{lib2roh.get(code, code)}"
            if reaction not in rows.index:
                continue
            profile = mech.calculate_mechanism(library_energies(imido, lib2roh.get(code, code)),
                                               standard_state=False)
            barriers = mech.calculate_barriers(profile)
            rates = kinetics.rates_from_barriers(barriers, settings["temperature"])
            outcome = ode_system.simulate(
                rates, config.initial_concentrations(condition), settings["simulation_time"])
            got = (outcome["bispyr"], outcome["map"], outcome["bis"])
            reference = (row["bispyr"], row["map"], row["bis"])
            same = all(abs(a - b) < 0.01 for a, b in zip(got, reference))
            ok += same
            total += 1
            print(f"  {'ok  ' if same else 'DIFF'} {condition} {label:<11} "
                  f"BP/MAP/BIS={tuple(round(v, 2) for v in got)} esperado={tuple(round(v, 2) for v in reference)}")
    if total:
        print(f"  -> {ok}/{total} simulações idênticas ao arquivo do autor\n")
    return ok == total


def main():
    results = [check_library(), check_dft(), check_simulation()]
    print("=" * 70)
    print("CONFERE com a tese" if all(results) else "DIVERGE: ver as linhas marcadas acima")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())