"""Roda UMA reação nas quatro condições padrão e confere contra o ESI.

É o caminho prático de conferência: uma reação × quatro condições × 1000 amostras
leva da ordem de 15 min nesta máquina, e é o que basta para saber se o programa
reproduz a tese. A biblioteca inteira são 136 reações × 4 condições, e as condições
forçadas (A2 e B2) custam ~10× mais que as brandas — rodar tudo é trabalho de muitas
horas, e não é o que este script faz.

Uso:
    python run_reaction.py im1-ROH3                    # as quatro condições
    python run_reaction.py im1-ROH3 --condicoes A1 B2
    python run_reaction.py im1-ROH3 --n-samples 200    # conferência rápida
    python run_reaction.py --listar                    # boas candidatas

"""

from argparse import ArgumentParser
from pathlib import Path
import sys
import time

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from src import compare, inputs, runner  # noqa: E402

ROOT = HERE.parents[0]
ESI = ROOT / "tables"   # as tabelas do material suplementar


def listar_candidates(esi_dir=ESI):
    """Reações da biblioteca com p_MAP perto do corte em alguma condição.

    São as que dizem alguma coisa: onde o produto previsto é óbvio (p_MAP 0 ou 1) a
    conferência passa por qualquer caminho.
    """
    reference = compare.load_library_reference(esi_dir)
    rows = []
    for reaction, row in reference.iterrows():
        values = [row[f"{c} p_MAP"] for c in compare.CONDITIONS]
        if any(value != value for value in values):
            continue
        margens = [compare.margin_of(value) for value in values]
        rows.append({"reação": reaction,
                     "p_MAP A1/A2/B1/B2": "/".join(f"{v:.2f}" for v in values),
                     "margem mínima": min(margens),
                     "classe": row.get("class")})
    table = pd.DataFrame(rows).sort_values("margem mínima")
    return table


def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=__doc__ and __import__("argparse").RawDescriptionHelpFormatter)
    parser.add_argument("reaction", nargs="?", help="rótulo da reação, por exemplo im1-ROH3")
    parser.add_argument("--condicoes", nargs="+", default=list(compare.CONDITIONS),
                        choices=list(compare.CONDITIONS), help="condições a rodar")
    parser.add_argument("--n-samples", type=int, default=1000, help="amostras por condição")
    parser.add_argument("--sigma", type=float, default=3.0, help="kcal/mol (3 composto, 2 DFT)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=-1, help="-1 = todos os núcleos")
    parser.add_argument("--saida", default=None, help="pasta da rodada")
    parser.add_argument("--listar", action="store_true", help="listar reações candidatas e sair")
    args = parser.parse_args()

    if args.listar:
        print(listar_candidates().head(20).to_string(index=False))
        return 0

    if not args.reaction:
        parser.error("informe a reação (ou use --listar)")

    reaction = args.reaction
    output_dir = Path(args.saida) if args.saida else HERE / f"exports/{reaction}"

    energies = inputs.library_energies(ESI, reactions=[reaction])
    if reaction not in energies:
        print(f"{reaction}: não achei o perfil completo no ESI")
        return 1

    # As energias do ESI já estão no estado padrão de 1 mol/L: nada de correção.
    profiles, barriers, _ = runner.profiles_from_energies(energies, standard_state=False)

    print(f"\n{reaction} — barreiras (kcal/mol):")
    print(pd.DataFrame(barriers).T.round(2).to_string())
    print(f"\n{len(args.condicoes)} condições × {args.n_samples} amostras, "
          f"semente {args.seed}, processos {'todos' if args.n_jobs == -1 else args.n_jobs}")

    started = time.time()
    frames = runner.run_conditions(
        profiles, conditions=tuple(args.condicoes), output_dir=output_dir,
        sigma=args.sigma, n_samples=args.n_samples, random_seed=args.seed,
        n_jobs=args.n_jobs, verbose=True)
    print(f"\ntempo total: {(time.time() - started) / 60:.1f} min")

    comparison = compare.compare_library(frames, ESI, conditions=tuple(args.condicoes))
    compare.report({f"{reaction} (biblioteca)": comparison},
                   title=f"Conferência de {reaction} contra o ESI")

    destino = output_dir / "conferencia.csv"
    comparison.to_csv(destino, index=False)
    print(f"\nrelatório: {destino}")
    return 0 if (comparison["verdict"] == "confere").all() else 1


if __name__ == "__main__":
    sys.exit(main())