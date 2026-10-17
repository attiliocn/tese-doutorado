"""Conferência de uma rodada contra os arquivos do ESI da tese.

O ESI é a fonte da verdade: `tables/ch2.microkineticsLibrary.csv` para as 136 reações da
biblioteca nas quatro condições padrão, `tables/ch2.microkineticsReported.csv` para as
reações publicadas nas condições de preparo, e
`tables/ch2.microkineticsReportedStandard.csv` para as mesmas nas condições padrão.

Duas perguntas diferentes, e o relatório separa as duas:

1. **reprodução** — o usuário rodou com as energias do próprio ESI? Então os números
   têm de sair iguais, porque a semente e o esquema de perturbação são os mesmos.
   Compara `p_MAP`, as medianas e o padrão do produto com tolerância apertada.
2. **concordância** — o usuário trouxe energias próprias? Então o que se compara é o
   **produto previsto** (o lado do corte de 0,5) e a **margem** (a distância a 0,5),
   que é o que a tese usa para decidir se a previsão depende do erro das energias.
   O p_MAP mede o erro do cálculo, não o acerto experimental.

Definição da p_MAP e do corte: README §4. O corte de 0,5 é a definição de majoritário;
não há categoria binária de confiança (decisão do autor, 16/10/2026).
"""

from pathlib import Path

import numpy as np
import pandas as pd

CONDITIONS = ("A1", "A2", "B1", "B2")
PRODUCT_CUTOFF = 0.5


def load_library_reference(esi_dir):
    """`tables/ch2.microkineticsLibrary.csv`, indexado pelo rótulo da reação."""
    frame = pd.read_csv(Path(esi_dir) / "ch2.microkineticsLibrary.csv")
    return frame.set_index("reaction")


def load_reported_reference(esi_dir):
    """As reações publicadas nas condições de preparo (DFT e protocolo composto)."""
    frame = pd.read_csv(Path(esi_dir) / "ch2.microkineticsReported.csv")
    return frame.set_index("reaction")


def load_reported_standard_reference(esi_dir):
    """As reações publicadas nas quatro condições padrão (só o protocolo DFT)."""
    frame = pd.read_csv(Path(esi_dir) / "ch2.microkineticsReportedStandard.csv")
    return frame.set_index("reaction")


def product_of(p_map):
    """Produto previsto: do lado do corte de 0,5. MAP acima, bisarilóxido abaixo."""
    if p_map != p_map:                      # NaN
        return None
    return "MAP" if p_map > PRODUCT_CUTOFF else "BIS"


def margin_of(p_map):
    """Margem da previsão: a distância ao corte. 0 é empate; 0,5 é o máximo."""
    if p_map != p_map:
        return float("nan")
    return abs(p_map - PRODUCT_CUTOFF)


def expected_pattern(reference_row):
    """O `pattern` do ESI como quatro letras, na ordem A1, A2, B1, B2."""
    pattern = reference_row.get("pattern")
    if isinstance(pattern, str) and pattern:
        return pattern
    letters = []
    for condition in CONDITIONS:
        value = reference_row.get(f"{condition} p_MAP", np.nan)
        letters.append("?" if value != value else ("M" if value > PRODUCT_CUTOFF else "B"))
    return "".join(letters)


def compare_library(run_frames, esi_dir, p_map_tolerance=0.05, median_tolerance=1.0,
                    conditions=CONDITIONS):
    """Compara a rodada da biblioteca com o ESI. Devolve um `DataFrame` por condição."""
    reference = load_library_reference(esi_dir)
    rows = []

    for condition in conditions:
        if condition not in run_frames:
            continue
        run = run_frames[condition]

        for reaction, run_row in run.iterrows():
            entry = {"reaction": reaction, "condition": condition}
            if reaction not in reference.index:
                entry.update(verdict="sem referência no ESI",
                             note="a reação não está no arquivo do ESI")
                rows.append(entry)
                continue

            reference_row = reference.loc[reaction]
            expected_map = reference_row.get(f"{condition} p_MAP", np.nan)
            entry["p_MAP_rodada"] = run_row["p_map"]
            entry["p_MAP_ESI"] = expected_map
            entry["p_MAP_diferenca"] = run_row["p_map"] - expected_map
            entry["produto_rodada"] = product_of(run_row["p_map"])
            entry["produto_ESI"] = product_of(expected_map)
            entry["margem_rodada"] = margin_of(run_row["p_map"])
            entry["margem_ESI"] = margin_of(expected_map)

            for quantity, column in (("MAP", "map_median"), ("BIS", "bis_median"),
                                     ("BP", "bispyr_median")):
                got = run_row.get(column, np.nan)
                want = reference_row.get(f"{condition} {quantity}", np.nan)
                entry[f"{quantity}_rodada"] = got
                entry[f"{quantity}_ESI"] = want
                entry[f"{quantity}_diferenca"] = got - want

            same_product = entry["produto_rodada"] == entry["produto_ESI"]
            close_map = abs(entry["p_MAP_diferenca"]) <= p_map_tolerance
            close_medians = all(
                abs(entry[f"{q}_diferenca"]) <= median_tolerance
                for q in ("MAP", "BIS", "BP")
                if entry[f"{q}_diferenca"] == entry[f"{q}_diferenca"])

            if same_product and close_map and close_medians:
                entry["verdict"] = "confere"
            elif same_product:
                entry["verdict"] = "produto confere, números divergem"
            else:
                entry["verdict"] = "produto divergente"
                entry["note"] = (f"rodada prevê {entry['produto_rodada']}; "
                                 f"ESI prevê {entry['produto_ESI']}")
            rows.append(entry)

    return pd.DataFrame(rows)


def compare_reported(run_frames, esi_dir, conditions=CONDITIONS):
    """Compara as reações publicadas nas condições padrão com o ESI.

    O `pattern` do arquivo é por condição; as entradas compostas (`8/9`, `14/15`)
    aparecem aqui pelo rótulo da reação, que nas condições padrão não tem sufixo de
    equivalentes.
    """
    reference = load_reported_standard_reference(esi_dir)
    rows = []

    for condition in conditions:
        if condition not in run_frames:
            continue
        run = run_frames[condition]
        for reaction, run_row in run.iterrows():
            entry = {"reaction": reaction, "condition": condition}
            if reaction not in reference.index:
                entry.update(verdict="sem referência no ESI")
                rows.append(entry)
                continue
            reference_row = reference.loc[reaction]
            expected_map = reference_row.get(f"{condition} p_MAP", np.nan)
            entry.update(p_MAP_rodada=run_row["p_map"], p_MAP_ESI=expected_map,
                         p_MAP_diferenca=run_row["p_map"] - expected_map,
                         produto_rodada=product_of(run_row["p_map"]),
                         produto_ESI=product_of(expected_map),
                         entry=reference_row.get("entry"))
            entry["verdict"] = ("confere" if entry["produto_rodada"] == entry["produto_ESI"]
                                else "produto divergente")
            rows.append(entry)

    return pd.DataFrame(rows)


def summary_of(comparison):
    """Contagens por veredito, para o cabeçalho do relatório."""
    if comparison.empty:
        return pd.Series(dtype=int)
    return comparison["verdict"].value_counts()


def report(comparisons, title="Conferência contra o ESI"):
    """Imprime o relatório e devolve o resumo. O usuário lê isto e nada mais."""
    print("=" * 78)
    print(title)
    print("=" * 78)
    overall = []
    for name, comparison in comparisons.items():
        if comparison is None or comparison.empty:
            continue
        counts = summary_of(comparison)
        total = len(comparison)
        agrees = int(counts.get("confere", 0))
        print(f"\n{name}: {agrees}/{total} linhas conferem")
        for verdict, count in counts.items():
            print(f"    {verdict:<34} {count}")
        disagreements = comparison[comparison["verdict"] != "confere"]
        if not disagreements.empty:
            print(f"    divergências:")
            for _, row in disagreements.head(15).iterrows():
                detail = f"p_MAP {row.get('p_MAP_rodada', float('nan')):.3f} x "
                detail += f"{row.get('p_MAP_ESI', float('nan')):.3f}"
                print(f"      {row['reaction']:<14} {row['condition']:<3} "
                      f"{row['verdict']:<32} {detail}")
            if len(disagreements) > 15:
                print(f"      ... e mais {len(disagreements) - 15}")
        overall.append((name, agrees, total))

    print("\n" + "=" * 78)
    agrees = sum(a for _, a, _ in overall)
    totals = sum(t for _, _, t in overall)
    print(f"TOTAL: {agrees}/{totals} linhas conferem com o ESI")
    print("=" * 78)
    return overall