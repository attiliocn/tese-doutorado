"""Propagação da incerteza das energias por Monte Carlo (Métodos, `sec:methodsMonteCarlo`).

Esquema correlacionado, como está nas equações `eq:mcIntermediates` e
`eq:mcTransitionStates` da tese:

    G'_I     = G_I + A,        A ~ N(0, sigma^2)      (MAP e bisarilóxido)
    G'_TS_j  = G_TS_j + x_j A, x_j ~ U(0, 1)          (um x por estado de transição)

Um único desvio A por amostra, comum aos dois intermediários; cada estado de transição
recebe uma fração aleatória dele; os reagentes são a referência e não são perturbados.
Como MAP e bis recebem o mesmo A, a energia de reação da segunda troca não muda entre
as amostras — mudam a energia da primeira troca e as quatro barreiras.

O que é perturbado é o **perfil relativo**, em kcal/mol, e não as energias absolutas:
é assim no pipeline do autor (as colunas do `*_perturbed_energies.csv` são
`reactants, a1, ts1, map, a2, ts2, bis`, com `reactants` sempre 0,0). Isso importa
porque o desvio é somado depois da conversão de hartree para kcal/mol, e porque σ
está em kcal/mol.

O sorteio tem de seguir esta ordem exata (A, depois x_ts1, depois x_ts2) e a semente
tem de ser a mesma: é o que torna a rodada comparável com os arquivos da tese.
"""

from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd

from . import config, ode_system


def draw_samples(sigma, n_samples, random_seed):
    """Sorteia as amostras de uma vez, na ordem exata do pipeline do autor.

    Separar o sorteio da simulação é o que permite paralelizar as integrações sem
    alterar o resultado: a sequência de números aleatórios não depende do paralelismo.
    """
    rng = np.random.default_rng(random_seed)
    samples = []
    for _ in range(n_samples):
        a_sample = float(rng.normal(0.0, sigma))
        x_ts1 = float(rng.uniform(0.0, 1.0))
        x_ts2 = float(rng.uniform(0.0, 1.0))
        samples.append({"A": a_sample, "x_ts1": x_ts1, "x_ts2": x_ts2,
                        "B_ts1": a_sample * x_ts1, "B_ts2": a_sample * x_ts2})
    return samples


def perturb_profile(profile, sample):
    """Aplica um desvio a um perfil relativo (kcal/mol). Devolve o perfil perturbado."""
    perturbed = dict(profile)
    for label in config.INTERMEDIATES:
        value = perturbed.get(label, np.nan)
        if value == value:
            perturbed[label] = value + sample["A"]
    for label in config.TRANSITION_STATES:
        value = perturbed.get(label, np.nan)
        if value == value:
            perturbed[label] = value + sample[f"B_{label}"]
    return perturbed


def _run_sample(job):
    """Uma amostra: perturba o perfil, recalcula barreiras e taxas, integra.

    Função de nível de módulo, para poder ser enviada a um processo de trabalho.
    Devolve `None` quando a integração passa do limite de tempo — a amostra é
    descartada, como na tese ("amostras cuja integração não terminou em 60 s foram
    descartadas"), e as descartadas ficam marcadas no `DataFrame` de saída.
    """
    profile, sample, settings, profile_to_barriers, barriers_to_rates = job
    barriers = profile_to_barriers(perturb_profile(profile, sample))
    rates = barriers_to_rates(barriers, settings["temperature"])
    return _with_timeout(ode_system.simulate, (rates, settings["concentrations"],
                                               settings["simulation_time"]),
                         config.JOB_TIMEOUT_SECONDS)


def _with_timeout(function, arguments, seconds):
    """Roda `function` com um limite de relógio de parede. `None` se estourar.

    Usa SIGALRM, então só vale em Linux e macOS — é o mesmo mecanismo do pipeline do
    autor, e é o que permite descartar uma amostra patológica sem perder a rodada.
    Precisa estar no processo de trabalho: o alarme não atravessa processos.
    """
    import signal

    def _raise(signum, frame):
        raise TimeoutError

    previous = signal.signal(signal.SIGALRM, _raise)
    signal.alarm(seconds)
    try:
        return function(*arguments)
    except TimeoutError:
        return None
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def run_monte_carlo(profile, profile_to_barriers, barriers_to_rates, settings,
                    sigma=None, n_samples=None, random_seed=None, n_jobs=-1,
                    quantities=("conversion", "map", "bis"), pool=None):
    """Linha de base mais as amostras perturbadas de uma reação.

    `profile` é o perfil relativo em kcal/mol (reagentes isolados = 0).
    `profile_to_barriers(profile)` devolve as barreiras b1..b4 e
    `barriers_to_rates(barriers, temperatura)` devolve k1d, k1r, k2d, k2r — passados
    de fora para que este módulo não dependa do formato das energias.

    `pool` é um `ProcessPoolExecutor` já aberto. Uma rodada de 136 reações × 4
    condições chama esta função 544 vezes: abrir um processo por chamada custava mais
    que as próprias integrações, por isso quem roda em lote abre o pool uma vez e
    passa aqui.

    Devolve `(linha_de_base, amostras)`: a linha de base é o perfil sem perturbação
    (equivalente ao `unperturbated.csv` do autor); as amostras são um `DataFrame` com
    o desvio de cada amostra e as grandezas pedidas.
    """
    sigma = config.PERTURBATION_SIGMA_DEFAULT if sigma is None else sigma
    n_samples = config.PERTURBATION_N_SAMPLES_DEFAULT if n_samples is None else n_samples
    random_seed = config.PERTURBATION_SEED_DEFAULT if random_seed is None else random_seed

    baseline_rates = barriers_to_rates(profile_to_barriers(profile), settings["temperature"])
    baseline = ode_system.simulate(baseline_rates, settings["concentrations"],
                                   settings["simulation_time"])

    if n_samples == 0:
        frame = pd.DataFrame(columns=["A", "x_ts1", "x_ts2", "B_ts1", "B_ts2", "timed_out", *quantities])
        frame.attrs["timed_out"] = []
        return baseline, frame

    samples = draw_samples(sigma, n_samples, random_seed)
    jobs = [(profile, sample, settings, profile_to_barriers, barriers_to_rates)
            for sample in samples]

    if n_jobs == 1:
        outcomes = [_run_sample(job) for job in jobs]
    elif pool is not None:
        outcomes = list(pool.map(_run_sample, jobs, chunksize=config.POOL_CHUNKSIZE))
    else:
        with ProcessPoolExecutor(max_workers=None if n_jobs in (-1, None) else n_jobs) as own_pool:
            outcomes = list(own_pool.map(_run_sample, jobs, chunksize=config.POOL_CHUNKSIZE))

    timed_out = [i for i, outcome in enumerate(outcomes) if outcome is None]
    frame = pd.DataFrame([
        {**{key: sample[key] for key in ("A", "x_ts1", "x_ts2", "B_ts1", "B_ts2")},
         "timed_out": outcome is None,
         **({} if outcome is None else {name: outcome[name] for name in quantities})}
        for sample, outcome in zip(samples, outcomes)
    ])
    frame.attrs["timed_out"] = timed_out
    return baseline, frame