"""Driver: energias -> perfis -> rodada nas condições padrão -> agregados por condição.

Substitui o que no pipeline do autor era feito por duas células de notebook
(`main.ipynb` células 6 e 12 a 15) e por um segundo notebook (`analysis.ipynb`) que só
sabia ler uma pasta de trajetórias. Aqui:

- as condições são explícitas, não globais mutáveis;
- os caminhos são argumentos;
- a agregação sai no **mesmo formato** do `aggregate_metrics.csv` do autor (as 35
  colunas), para que os dois conjuntos possam ser comparados célula a célula;
- cada condição grava um `run_options.json` com tudo o que é preciso para reproduzir.

O `p_MAP` é a fração das amostras em que o MAP é o majoritário, isto é, o lado do
corte de 0,5 — que é a definição de majoritário, não uma convenção (README §4).
As amostras empatadas não entram no denominador, como no `analysis.ipynb`.
"""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json
import time

import numpy as np
import pandas as pd

from . import config, kinetics, mechanism, montecarlo, ode_system

# Colunas do `aggregate_metrics.csv` do autor, na ordem: as estatísticas de cada
# grandeza, depois os contadores e as duas probabilidades.
AGGREGATED_QUANTITIES = ("conversion", "roh", "bispyr", "map", "bis", "f_map")
STATISTICS = ("avg", "std", "median", "CI_low", "CI_up")
CI_PERCENTILES = (2.5, 97.5)


def profiles_from_energies(energies, standard_state=True, verbose=True):
    """`{reação: energias}` -> `{reação: perfil}` e `{reação: barreiras}`.

    A montagem é a do parser do autor (`calculate_mechanism` do `exchange-parser`),
    inclusive a correção de estado padrão, que é somada **à energia de cada espécie**
    antes de montar o perfil.
    """
    profiles, barriers, skipped = {}, {}, {}
    for reaction, block in energies.items():
        try:
            mechanism.check_species(block, label=reaction)
        except ValueError as error:
            skipped[reaction] = str(error)
            continue
        profile, barrier = mechanism.profile_from_energies(block, standard_state=standard_state)
        profiles[reaction] = profile
        barriers[reaction] = barrier

    if verbose and skipped:
        print(f"  {len(skipped)} reação(ões) fora da rodada, por falta de espécie:")
        for reaction, why in list(skipped.items())[:5]:
            print(f"    - {why}")
    if verbose:
        print(f"  {len(profiles)} reações com perfil completo")
    return profiles, barriers, skipped


def settings_for(condition, temperature=None, simulation_time=None, concentrations=None):
    """Os parâmetros de uma condição padrão, permitindo trocar o que vier de fora.

    As reações publicadas usam as condições de preparo de cada experimento, que não
    são as quatro padrão — daí os três parâmetros opcionais.
    """
    base = config.STANDARD_CONDITIONS[condition]
    return {
        "condition": condition,
        "temperature": base["temperature"] if temperature is None else temperature,
        "simulation_time": base["simulation_time"] if simulation_time is None else simulation_time,
        "concentrations": config.initial_concentrations(condition) if concentrations is None else concentrations,
    }


def _profile_to_barriers(profile):
    return mechanism.calculate_barriers(profile)


def _barriers_to_rates(barriers, temperature):
    return kinetics.rates_from_barriers(barriers, temperature)


def summarise(samples, baseline, quantity_columns=("conversion", "map", "bis")):
    """Média, desvio, mediana e intervalo de 95 % de cada grandeza.

    As amostras descartadas por estouro de tempo ficam fora de tudo, inclusive dos
    contadores — é o que a tese faz ao descartá-las. O número de descartadas vai no
    resumo, porque ele diz se a rodada daquela reação é confiável.
    """
    if "timed_out" in samples.columns:
        discarded = int(samples["timed_out"].sum())
        frame = samples[~samples["timed_out"]].copy()
    else:
        discarded, frame = 0, samples.copy()

    frame["f_map"] = frame["map"] / (frame["map"] + frame["bis"])

    summary = {"n_samples": float(len(frame)), "n_timed_out": float(discarded)}
    for quantity in AGGREGATED_QUANTITIES:
        if quantity not in frame.columns:
            continue
        values = frame[quantity]
        for name in STATISTICS:
            key = f"{quantity}_{name}"
            if name == "avg":
                summary[key] = float(values.mean())
            elif name == "std":
                summary[key] = float(values.std())
            elif name == "median":
                summary[key] = float(values.median())
            else:
                summary[key] = float(np.percentile(values, CI_PERCENTILES[0 if name == "CI_low" else 1]))

    map_majority = int((frame["map"] > frame["bis"]).sum())
    bis_majority = int((frame["bis"] > frame["map"]).sum())
    decided = map_majority + bis_majority
    summary["counts_map_major_product"] = float(map_majority)
    summary["counts_bis_major_product"] = float(bis_majority)
    summary["p_map"] = float(map_majority / decided) if decided else float("nan")
    summary["p_bis"] = float(bis_majority / decided) if decided else float("nan")
    summary["_baseline"] = baseline
    return summary


def run_reaction(profile, settings, sigma=None, n_samples=None, random_seed=None,
                 n_jobs=-1, quantities=("conversion", "map", "bis", "roh", "bispyr"),
                 pool=None):
    """Uma reação numa condição. Devolve `(resumo, linha_de_base, amostras)`."""
    baseline, samples = montecarlo.run_monte_carlo(
        profile, _profile_to_barriers, _barriers_to_rates, settings,
        sigma=sigma, n_samples=n_samples, random_seed=random_seed, n_jobs=n_jobs,
        quantities=quantities, pool=pool)
    return summarise(samples, baseline), baseline, samples


MODE_FULL = "full"
MODE_FAST = "fast"

# No modo rápido a linha de base entra no denominador: a amostra não perturbada é o
# limite do conjunto, e usá-la dá um p_MAP estável com poucas amostras. Uma integração
# custa ~130 ms, então 1000 amostras × 4 condições × 136 reações são ~19 h seriais.
def run_conditions(profiles, conditions=("A1", "A2", "B1", "B2"), output_dir=None,
                   settings_by_reaction=None, sigma=None, n_samples=None,
                   random_seed=None, n_jobs=-1, save_trajectories=False,
                   save_samples=True, mode=MODE_FULL, verbose=True):
    """Roda um conjunto de perfis nas condições pedidas e grava os resultados.

    `settings_by_reaction` permite dar a cada reação uma condição própria (as reações
    publicadas): é um mapeamento `{reação: {temperature, simulation_time, concentrations}}`.
    Sem ele, todas as reações usam a condição padrão.

    `mode`: `full` é o da tese — `n_samples` amostras perturbadas, com a linha de base
    fora do denominador do p_MAP; `fast` inclui a linha de base entre as amostras, para
    descoberta e comparação rápida. O padrão da tese é o `full`.

    Devolve `{condição: DataFrame}` com uma linha por reação, no formato do
    `aggregate_metrics.csv` do autor.
    """
    output_dir = Path(output_dir) if output_dir else None
    sigma = config.PERTURBATION_SIGMA_DEFAULT if sigma is None else sigma
    n_samples = config.PERTURBATION_N_SAMPLES_DEFAULT if n_samples is None else n_samples
    random_seed = config.PERTURBATION_SEED_DEFAULT if random_seed is None else random_seed

    # Um pool para a rodada inteira. Abrir um processo por reação custava mais que as
    # próprias integrações: são 544 chamadas numa rodada da biblioteca.
    pool = None
    if n_jobs != 1 and _pool_needed(len(profiles), n_samples):
        pool = ProcessPoolExecutor(max_workers=None if n_jobs in (-1, None) else n_jobs)
        if verbose:
            print(f"processos de trabalho: {pool._max_workers}")

    try:
        return _run_conditions(profiles, conditions, output_dir, settings_by_reaction,
                               sigma, n_samples, random_seed, n_jobs, save_trajectories,
                               save_samples, mode, verbose, pool)
    finally:
        if pool is not None:
            pool.shutdown()


def _pool_needed(n_reactions, n_samples):
    """Vale abrir processos? Só quando há amostras suficientes para amortizar o custo."""
    return n_reactions * n_samples > 20


def _run_conditions(profiles, conditions, output_dir, settings_by_reaction, sigma,
                    n_samples, random_seed, n_jobs, save_trajectories, save_samples,
                    mode, verbose, pool):
    results = {}
    for condition in conditions:
        condition_dir = None
        if output_dir:
            condition_dir = output_dir / condition
            condition_dir.mkdir(parents=True, exist_ok=True)
            (condition_dir / "trajectories").mkdir(exist_ok=True)

        rows, started = [], time.time()
        if verbose:
            print(f"\n--- {condition}: {len(profiles)} reações, "
                  f"{n_samples} amostras, sigma = {sigma} kcal/mol, semente {random_seed}")

        for index, (reaction, profile) in enumerate(profiles.items(), start=1):
            overrides = (settings_by_reaction or {}).get(reaction, {})
            settings = settings_for(condition, **overrides)

            summary, baseline, samples = run_reaction(
                profile, settings, sigma=sigma, n_samples=n_samples,
                random_seed=random_seed, n_jobs=n_jobs, pool=pool)
            if mode == MODE_FAST:
                samples = _add_baseline_sample(samples, baseline)
                summary = summarise(samples, baseline)

            rows.append({"trajectory_alias": reaction, **{k: v for k, v in summary.items() if k != "_baseline"}})

            if condition_dir and save_samples:
                samples.to_csv(condition_dir / f"{reaction}_samples.csv", index=False)
                pd.DataFrame([{**{"sample": -1, "is_baseline": 1}, **baseline}]).to_csv(
                    condition_dir / f"{reaction}_baseline.csv", index=False)
            if condition_dir and save_trajectories:
                _save_trajectory(condition_dir / "trajectories", reaction, settings, profile,
                                 sigma, n_samples, random_seed, n_jobs)

            if verbose and (index % 10 == 0 or index == len(profiles)):
                elapsed = time.time() - started
                remaining = elapsed / index * (len(profiles) - index)
                print(f"    {index}/{len(profiles)}  ({elapsed:.0f} s decorridos, "
                      f"~{remaining:.0f} s restantes)")

        frame = pd.DataFrame(rows).set_index("trajectory_alias")
        results[condition] = frame

        if condition_dir:
            frame.to_csv(condition_dir / "aggregate_metrics.csv")
            _write_run_options(condition_dir, condition, settings_for(condition),
                               sigma, n_samples, random_seed, list(profiles),
                               save_trajectories, mode)
            if verbose:
                print(f"    gravado em {condition_dir}/aggregate_metrics.csv")

    return results


def load_run_frames(directory, conditions=("A1", "A2", "B1", "B2")):
    """Lê uma rodada gravada: `{condição: DataFrame}` a partir dos `aggregate_metrics.csv`.

    É o que o `validate.ipynb` usa — as duas pontas (rodar e conferir) se comunicam
    pelo disco, e não por estado em memória.
    """
    directory = Path(directory)
    frames = {}
    for condition in conditions:
        path = directory / condition / "aggregate_metrics.csv"
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        if "trajectory_alias" in frame.columns:
            frame = frame.set_index("trajectory_alias")
        frames[condition] = frame
    if not frames:
        raise FileNotFoundError(
            f"{directory}: nenhum aggregate_metrics.csv encontrado. "
            f"Rode o run.ipynb primeiro, ou confira o nome da pasta da rodada."
        )
    return frames


def estimate_runtime(n_reactions, conditions=("A1", "A2", "B1", "B2"), n_samples=None,
                     n_jobs=-1, seconds_per_sample=0.13):
    """Estimativa grosseira do custo de uma rodada, em minutos.

    Uma integração custa ~0,13 s nesta máquina; o paralelismo divide por
    `min(n_jobs, núcleos)`. Serve para o usuário decidir antes de soltar o laço.
    """
    import os
    n_samples = config.PERTURBATION_N_SAMPLES_DEFAULT if n_samples is None else n_samples
    workers = os.cpu_count() if n_jobs in (-1, None) else min(n_jobs, os.cpu_count())
    workers = max(workers, 1)
    total_samples = n_reactions * len(conditions) * (n_samples + 1)
    return total_samples * seconds_per_sample / workers / 60


def _add_baseline_sample(samples, baseline):
    """Acrescenta a linha de base ao conjunto de amostras (modo rápido)."""
    row = {key: np.nan for key in samples.columns}
    row.update({"timed_out": False, **baseline})
    return pd.concat([samples, pd.DataFrame([row])], ignore_index=True)


def _write_run_options(directory, condition, settings, sigma, n_samples, random_seed,
                       reactions, save_trajectories, mode=MODE_FULL):
    """O registro do que foi rodado. Sem ele a rodada não é reproduzível."""
    options = {
        "condition": condition,
        "temperature_C": settings["temperature"],
        "simulation_time_s": settings["simulation_time"],
        "initial_concentrations": settings["concentrations"],
        "sigma_kcal_mol": sigma,
        "n_samples": n_samples,
        "random_seed": random_seed,
        "mode": mode,
        "n_reactions": len(reactions),
        "reactions": reactions,
        "save_trajectories": bool(save_trajectories),
        "integration": {"solver": config.ODE_SOLVER_OPTIONS["method"],
                        "rtol": config.ODE_SOLVER_OPTIONS["rtol"],
                        "atol": config.ODE_SOLVER_OPTIONS["atol"],
                        "points": config.TRAJECTORY_POINTS_DEFAULT,
                        "exponent": config.TRAJECTORY_EXPONENT_DEFAULT},
        "perturbation": {
            "intermediates": list(config.INTERMEDIATES),
            "transition_states": list(config.TRANSITION_STATES),
            "scheme": "A ~ N(0, sigma) nos intermediarios; x_j A, x_j ~ U(0,1), em cada estado de transicao",
        },
    }
    (directory / "run_options.json").write_text(json.dumps(options, indent=1, ensure_ascii=False))


def _save_trajectory(directory, reaction, settings, profile, sigma, n_samples,
                     random_seed, n_jobs):
    """Grava os pontos finais das amostras. A trajetória completa só se pedida.

    O pipeline do autor gravava os 1001 pontos de tempo de todas as amostras
    (~28 MB por reação por condição; 36 GB no conjunto das rodadas dele). Aqui o
    padrão é não gravar: o resumo e as amostras bastam para a conferência.
    """
    _, samples = montecarlo.run_monte_carlo(
        profile, _profile_to_barriers, _barriers_to_rates, settings,
        sigma=sigma, n_samples=n_samples, random_seed=random_seed, n_jobs=n_jobs,
        quantities=("conversion", "map", "bis", "roh", "bispyr"))
    samples.insert(0, "reaction", reaction)
    samples.to_pickle(directory / f"{reaction}_samples.pkl.gz", compression="gzip")