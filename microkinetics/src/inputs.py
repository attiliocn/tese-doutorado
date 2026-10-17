"""Entradas: montar o arquivo de energias a partir dos arquivos do ESI da tese.

O `run.ipynb` aceita um arquivo de energias do usuário (`reaction, species, G_hartree`).
Para quem quer partir do que a tese já tem, este módulo monta esse arquivo a partir de
`tables/ch2.libraryBestConformers.csv` (a biblioteca, 1 mol/L) ou da tabela de energias
absolutas das reações publicadas (1 atm).

A montagem é a mesma que a do `check_against_thesis.py`: pega a espécie de menor
energia livre de cada papel. Na biblioteca isso é o próprio conteúdo do arquivo (uma
linha por espécie); nas reações publicadas, a tabela já traz um confôrmero por espécie.

As energias saem **como estão no arquivo de origem**: quem chama diz se elas já estão
no estado padrão de 1 mol/L (`standard_state=False`, o caso da biblioteca) ou não
(`standard_state=True`, o caso da tabela a 1 atm).
"""

from pathlib import Path
import re

import pandas as pd

# Espécies na ordem em que o arquivo de energias é escrito.
SPECIES_ORDER = ("bispyrrolide", "roh", "pyrrole", "a1", "ts1", "map", "a2", "ts2", "bis")

PROFILE_POINTS = ("a1", "ts1", "map", "a2", "ts2", "bis")


def _lowest(energies, prefix):
    """Energia da espécie `prefix`, aceitando o rótulo exato ou `prefix_<índice>`.

    O sufixo `_` é obrigatório no segundo caso: sem ele `roh_ROH3` casaria também com
    `roh_ROH30_1`, e o arquivo passaria a pegar o fenol errado.
    """
    candidates = [v for k, v in energies.items()
                  if k == prefix or k.startswith(prefix + "_")]
    return min(candidates) if candidates else None


def library_energies(esi_dir, reactions=None, phenol_codes=None):
    """Energias da biblioteca a partir do ESI. Devolve `{reação: {espécie: G}}`.

    `reactions` limita o conjunto (rótulos `im<imido>-ROH<n>`); sem ele, sai a
    biblioteca inteira. Os rótulos já são os da tese, então não há código antigo
    para traduzir aqui.
    """
    esi_dir = Path(esi_dir)
    table = pd.read_csv(esi_dir / "ch2.libraryBestConformers.csv")
    table.columns = [c.strip() for c in table.columns]
    energies = dict(zip(table["standardized label"], table["G"]))

    labels = [r for r in (reactions or _reactions_from(energies))]
    result = {}
    for reaction in labels:
        imido, roh = reaction.split("-", 1)
        block = {
            "bispyrrolide": _lowest(energies, f"bispyrrolide_{imido}"),
            "roh": _lowest(energies, f"roh_{roh}"),
            "pyrrole": _lowest(energies, "pyrrole_dimethylpyrrole"),
            **{s: _lowest(energies, f"{imido}-{roh}_{s}") for s in PROFILE_POINTS},
        }
        if block["bispyrrolide"] is None or block["roh"] is None:
            continue
        result[reaction] = block
    return result


# Rótulo de reação + ponto do mecanismo. O índice do confôrmero é opcional: os rótulos
# da biblioteca trazem (`im0-ROH3_map_33`) e os das reações publicadas não
# (`im1-ROH20_map`), e os dois arquivos passam por aqui.
REACTION_LABEL = re.compile(r"^(im\d+-[^_]+)_(?:a1|ts1|map|a2|ts2|bis)(?:_\d+)?$")


def _reactions_from(energies):
    """Rótulos de reação presentes nas chaves do arquivo, na ordem de aparição."""
    seen = []
    for key in energies:
        match = REACTION_LABEL.match(key)
        if match and match.group(1) not in seen:
            seen.append(match.group(1))
    return seen


def published_energies(esi_dir, reactions=None):
    """Energias das reações publicadas, a partir do `tables/ch2.RealSystemBestConformers.csv`.

    As energias livres saem **a 1 atm** (é a coluna `G` da tabela do Apêndice D, como
    saiu do ORCA), então quem chama usa `standard_state=True` para converter.

    No arquivo convivem duas convenções de rótulo: as espécies do mecanismo usam o
    código do ESI (`im1-ROH20_a1`) e as compartilhadas, o código antigo (`roh_ad1`,
    `bispyrrolide_im1`). A ponte é o `ch2.phenolCodes.csv`.
    """
    esi_dir = Path(esi_dir)
    table = pd.read_csv(esi_dir / "ch2.RealSystemBestConformers.csv")
    energy = dict(zip(table["custom_label"], table["G_hartree"]))

    bridge = {}
    for row in pd.read_csv(esi_dir / "ch2.phenolCodes.csv").to_dict("records"):
        esi_code = str(row["ROH code in ESI labels"]).strip()
        bridge.setdefault(esi_code, []).append(str(row["library code"]).strip())

    labels = list(reactions) if reactions else _reactions_from(energy)
    result = {}
    for reaction in labels:
        imido, roh = reaction.split("-", 1)
        old_codes = bridge.get(roh, []) + bridge.get(roh + "b", [])
        roh_energy = next((energy[f"roh_{old}"] for old in old_codes
                           if f"roh_{old}" in energy), None)
        block = {
            "bispyrrolide": energy.get(f"bispyrrolide_{imido}"),
            "roh": roh_energy,
            "pyrrole": energy.get("dimethylpyrrole"),
            **{s: energy.get(f"{imido}-{roh}_{s}") for s in PROFILE_POINTS},
        }
        if block["bispyrrolide"] is None or block["roh"] is None:
            continue
        result[reaction] = {k: v for k, v in block.items() if v is not None}
    return result


def write_energies_file(energies, path, value_column="G_hartree"):
    """Grava o arquivo de energias no formato que o `run.ipynb` lê."""
    rows = []
    for reaction, block in energies.items():
        for species in SPECIES_ORDER:
            if species in block and block[species] is not None:
                rows.append({"reaction": reaction, "species": species,
                             value_column: block[species]})
    frame = pd.DataFrame(rows)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, float_format="%.8f")
    return len(frame), frame["reaction"].nunique()