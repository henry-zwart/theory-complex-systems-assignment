from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import seaborn as sns
from scipy.optimize import minimize

from tcs_assignment import configure_mpl

RANDOM_SEED = 42
DATA_PATH = Path("data")

configure_mpl()


@dataclass
class SupremeCourtData:
    votes: npt.NDArray[np.int64]
    fit_h: npt.NDArray[np.float64]
    fit_j: npt.NDArray[np.float64]

    @classmethod
    def load(cls) -> SupremeCourtData:
        votes_path = DATA_PATH / "supreme_court_votes.txt"
        fit_h_path = DATA_PATH / "fit_parameters_h.txt"
        fit_j_path = DATA_PATH / "fit_parameters_j.txt"

        with votes_path.open("r") as f:
            votes = np.asarray(
                [[int(vote) for vote in row.strip()] for row in f.readlines()]
            )
            votes = (votes * 2) - 1

        with fit_h_path.open("r") as f:
            fit_h = np.asarray([float(h.strip()) for h in f.readlines()])

        n = fit_h.size
        fit_j = np.zeros((n, n), dtype=np.float64)
        with fit_j_path.open("r") as f:
            upper_diag_j_flat = np.asarray([float(j.strip()) for j in f.readlines()])
            preceding_elts = 0
            for i in range(n - 1):
                row_elts = (n - i) - 1
                fit_j[i, i + 1 :] = upper_diag_j_flat[
                    preceding_elts : preceding_elts + row_elts
                ]
                preceding_elts += row_elts

        return SupremeCourtData(votes, fit_h, fit_j)


def unpack_param_vec(
    g: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    for n in range(g.size):
        if n * (n + 1) / 2 == g.size:
            break
    else:
        raise RuntimeError("Invalid length for parameter vector 'g'.")

    h = g[:n]
    flat_j = g[n:]
    j = np.zeros((n, n), dtype=np.float64)
    preceding_elts = 0
    for i in range(n - 1):
        row_elts = (n - i) - 1
        j[i, i + 1 :] = flat_j[preceding_elts : preceding_elts + row_elts]
        preceding_elts += row_elts

    return h, j


def flatten_params(
    h: npt.NDArray[np.float64], j: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    n = h.size
    g = np.empty(n * (n + 1) // 2, dtype=np.float64)
    g[:n] = h
    g[n:] = extract_flat_above_diag(j)
    return g


def thermal_average(
    fn: Callable[[npt.NDArray[np.int64]], float | int],
    h: npt.NDArray[np.float64],
    j: npt.NDArray[np.float64],
    partition: float | None = None,
) -> float:
    """Compute the thermal average of a quantity for a fit model.

    The thermal average is calculated as the expected value of the
    quantity applied to each possible microstate.

    Args:
        fn: Function which accepts a microstate and returns a numeric value.
        h: 1D numpy array of fit h_i parameters.
        j: 2D upper-triangular matrix of fit J_(i j) parameters.
        partition: Optional pre-computed normalisation constant for microstate
            distribution.

    Returns:
        Thermal average of the input quantity (fn) as a float.
    """
    if partition is None:
        partition = calculate_partition(h, j)

    result = 0
    for s in enumerate_microstates(h.size):
        value = fn(s)
        if np.isclose(value, 0):
            continue

        state_probability = boltzmann_probability(s, h, j, partition)
        result += state_probability * value

    return result


def empirical_avg_local_magnetisation(
    votes: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    """Compute the local average magnetisation per spin for a dataset.

    The empirical local average magnetisation for a given spin is
    calculated as the average value of that spin across the rows in
    the data.

    Args:
        votes: Numpy array of shape (court issue, judge) where the rows
            correspond to distinct hearings, and the columns correspond
            to each judge's votes.

    Returns:
        1D numpy array with length equal to the number of judges (spins),
        where each element is the empirical local average magnetisation of
        the corresponding spin in the dataset.
    """
    return votes.mean(axis=0)


def empirical_avg_local_correlation(
    votes: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    """Compute the local average correlation between spin pairs in a dataset.

    The empirical local average correlation for a pair of spins is
    calculated as the average value of s_i * s_j across the rows in
    the data.

    This is calculated as 1/N X.T @ X, where X is the dataset.

    Args:
        votes: Numpy array of shape (court issue, judge) where the rows
            correspond to distinct hearings, and the columns correspond
            to each judge's votes.

    Returns:
        2D numpy array with shape (n judges, n judges) where entry (i,j)
        is the empirical local average correlation between the i'th and j'th judges.
    """
    n = votes.shape[0]
    return ((1 / n) * (votes.T @ votes)).astype(np.float64)


def empirical_probability_conservative_vote(
    votes: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    """Calculate probability (per judge) of voting conservative in a dataset.

    Args:
        votes: 2D numpy array of votes (+1 conservative, -1 liberal),
            where rows correspond to issues, and columns correspond to
            the vote of a given justice on a particular issue.

    Returns:
        1D numpy array with empirical probabilities of voting conservative for each
        justice.
    """
    return (votes > 0).sum(axis=0) / votes.shape[0]


def model_avg_local_magnetisation(
    h: npt.NDArray[np.float64],
    j: npt.NDArray[np.float64],
    partition: float | None = None,
) -> npt.NDArray[np.float64]:
    """Compute the model local average magnetisation per spin.

    The local average magnetisation for a given spin is calculated as
    the thermal average value of that spin over all possible microstates,
    weighted by the probability of observing each microstate.

    Args:
        h: 1D numpy array of fit h_i parameters.
        j: 2D upper-triangular matrix of fit J_(i j) parameters.
        partition: Optional pre-computed normalisation constant for microstate
            distribution.

    Returns:
        1D numpy array with length equal to the number of judges (spins),
        where each element is the local average magnetisation of the corresponding
        spin in the fit model.
    """
    if partition is None:
        partition = calculate_partition(h, j)

    s_i = np.zeros_like(h)
    for s in enumerate_microstates(h.size):
        s_i += s * boltzmann_probability(s, h, j, partition)

    return s_i


def model_avg_local_correlation(
    h: npt.NDArray[np.float64],
    j: npt.NDArray[np.float64],
    partition: float | None = None,
) -> npt.NDArray[np.floating[Any]]:
    """Compute the local average correlation between spin pairs in a model.

    The model local average correlation for a pair of spins is
    calculated as the average value of s_i * s_j across all microstates,
    weighted by the microstate probability in the fit model.

    Args:
        h: 1D numpy array of fit h_i parameters.
        j: 2D upper-triangular matrix of fit J_(i j) parameters.
        partition: Optional pre-computed normalisation constant for microstate
            distribution.

    Returns:
        2D numpy array with shape (n judges, n judges) where entry (i,j)
        is the local average correlation between the i'th and j'th judges in the
        fit model.
    """
    if partition is None:
        partition = calculate_partition(h, j)

    si_sj = np.zeros_like(j)
    for s in enumerate_microstates(h.size):
        outer = np.outer(s, s)
        si_sj += outer * boltzmann_probability(s, h, j, partition)

    return si_sj


def extract_flat_above_diag(a):
    n = a.shape[0]
    above_diag_idxs = np.triu_indices(n, k=1)
    return a[above_diag_idxs]


def calculate_microstate_energy(
    s: npt.NDArray[np.int64],
    h: npt.NDArray[np.float64],
    j: npt.NDArray[np.float64],
) -> float:
    field_energy = -np.dot(s, h)
    interaction_energy = -(j * np.outer(s, s)).sum()
    return field_energy + interaction_energy


def enumerate_microstates(n: int):
    for i in range(2**n):
        binary_repr = np.binary_repr(i, width=n)
        microstate = np.asarray(list(map(int, binary_repr)), dtype=np.int64)
        microstate = microstate * 2 - 1
        yield microstate


def calculate_partition(
    h: npt.NDArray[np.float64],
    j: npt.NDArray[np.float64],
) -> float:
    z = 0
    # for microstate in enumerate_microstates(n):
    n = h.size
    for i in range(2**n):
        binary_repr = np.binary_repr(i, width=n)
        microstate = np.asarray(list(map(int, binary_repr)), dtype=np.int64)
        microstate = microstate * 2 - 1
        z += np.exp(
            -calculate_microstate_energy(
                microstate,
                h,
                j,
            )
        )

    return z


def boltzmann_probability(
    s: npt.NDArray[np.int64],
    h: npt.NDArray[np.float64],
    j: npt.NDArray[np.float64],
    partition: float | None = None,
) -> float:
    if partition is None:
        partition = calculate_partition(h, j)

    energy = calculate_microstate_energy(s, h, j)
    log_prob = -energy - np.log(partition)
    prob = np.exp(log_prob)
    return prob


def empirical_probability(
    s: npt.NDArray[np.int64],
    dataset: npt.NDArray[np.int64],
) -> float:
    occurrences = ((dataset * s).sum(axis=1) == s.size).sum()
    return occurrences / dataset.shape[0]


def probability_k_conservative_independent_recursive(
    k: int, p: npt.NDArray[np.float64], n: int
) -> float:
    i = len(p) - n

    if k > n:  # Require more conservative votes than remaining judges
        return 0
    elif k == 0:  # Rest of judges vote liberal
        return np.exp(np.log(1 - p[i:]).sum())
    elif k == n:
        return np.exp(np.log(p[i:]).sum())
    elif k < n:
        case1 = p[i] * probability_k_conservative_independent_recursive(k - 1, p, n - 1)
        case2 = (1 - p[i]) * probability_k_conservative_independent_recursive(
            k, p, n - 1
        )
        return case1 + case2
    else:
        raise RuntimeError("Huh?")


def probability_k_conservative_independent(k: int, p: npt.NDArray[np.float64]) -> float:
    return probability_k_conservative_independent_recursive(k, p, n=len(p))


def probability_k_conservative_data(k: int, votes: npt.NDArray[np.int64]) -> float:
    net_vote = votes.sum(axis=1)
    k_conversvative_count = (net_vote == 2 * k - votes.shape[1]).sum()
    return k_conversvative_count / votes.shape[0]


def probability_k_conservative_model(
    k: int,
    h: npt.NDArray[np.float64],
    j: npt.NDArray[np.float64],
) -> float:
    return thermal_average(lambda s: (s.sum() == 2 * k - h.size).astype(int), h, j)


def negative_log_likelihood(
    votes: npt.NDArray[np.int64],
    h: npt.NDArray[np.float64],
    j: npt.NDArray[np.float64],
    partition: float,
    data_distribution: npt.NDArray[np.float64] | None = None,
) -> float:
    # Compute p_d if not provided
    if data_distribution is None:
        data_distribution = np.array(
            [
                empirical_probability(s, votes)
                for s in enumerate_microstates(votes.shape[1])
            ]
        )
    log_likelihood = 0
    for i, s in enumerate(enumerate_microstates(votes.shape[1])):
        p_d = data_distribution[i]
        p_g = boltzmann_probability(s, h, j, partition)
        log_likelihood += p_d * np.log(p_g)
    nll = -votes.shape[0] * log_likelihood
    return nll


def fit_model(
    votes: npt.NDArray[np.int64],
) -> tuple[
    tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]], npt.NDArray[np.float64]
]:
    nll_values = []

    # Pre-compute distribution of microstates in dataset
    p_d = np.array(
        [empirical_probability(s, votes) for s in enumerate_microstates(votes.shape[1])]
    )

    def objective_fn(g) -> float:
        h, j = unpack_param_vec(g)
        partition = calculate_partition(h, j)
        nll = negative_log_likelihood(votes, h, j, partition, p_d)
        nll_values.append(nll)
        return nll

    h_init = empirical_avg_local_magnetisation(votes)
    j_init = empirical_avg_local_correlation(votes)
    g_init = flatten_params(h_init, j_init)

    result = minimize(objective_fn, x0=g_init)

    return unpack_param_vec(result.x), np.asarray(nll_values)


def main():
    data = SupremeCourtData.load()
    votes = data.votes

    # Compute avg local magnetisation <s_i>_D and  correlation <s_i s_j>_D
    si_data = empirical_avg_local_magnetisation(votes)
    si_data_sort_idx = np.argsort(si_data)
    si_sj_data = empirical_avg_local_correlation(votes)

    fig, ax = plt.subplots(figsize=(3.5, 2), constrained_layout=True)
    ax.scatter(np.arange(votes.shape[1]), si_data[si_data_sort_idx], marker="x")
    ax.axhline(y=0, linestyle="dashed", linewidth=0.5, color="grey")
    ax.set_xticks(np.arange(votes.shape[1]))
    ax.set_yticks(np.linspace(-0.6, 0.6, 7))
    ax.set_ylim(-0.6, 0.6)
    ax.set_xlabel("Judge $i$")
    ax.set_ylabel(r"$\langle s_i \rangle_D$")
    ax.set_title("Average vote per judge (increasing order)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(
        "results/figures/empirical_avg_local_mag_sorted.svg",
        bbox_inches="tight",
        transparent=True,
    )

    fig, ax = plt.subplots(figsize=(4, 3.2), constrained_layout=True)
    sns.heatmap(
        si_sj_data[si_data_sort_idx][:, si_data_sort_idx],
        square=True,
        vmin=0,
        vmax=1,
        cmap="flare",
        cbar_kws={"label": r"$\langle s_i s_j \rangle_D$"},
        ax=ax,
    )
    ax.set_xlabel(r"Judge (increasing $\langle s_i \rangle_D$)")
    ax.set_ylabel(r"Judge (increasing $\langle s_i \rangle_D$)")
    ax.set_title(r"Empirical local correlation between judges")
    fig.savefig(
        "results/figures/empirical_avg_local_corr_sorted_by_mag.svg",
        bbox_inches="tight",
        transparent=True,
    )

    fig, ax = plt.subplots(figsize=(3.5, 2), constrained_layout=True)
    ax.scatter(np.arange(data.fit_h.size), data.fit_h[si_data_sort_idx], marker="x")
    ax.axhline(y=0, linestyle="dashed", linewidth=0.5, color="grey")
    ax.set_xticks(np.arange(votes.shape[1]))
    ax.set_yticks(np.linspace(-0.6, 0.6, 7))
    ax.set_ylim(-0.6, 0.6)
    ax.set_xlabel("Judge $i$")
    ax.set_ylabel(r"$h_i$")
    ax.set_title(r"Fit $h_i$ parameters (increasing $\langle s_i \rangle_D$)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(
        "results/figures/fit_parameter_h_scatterplot.svg",
        bbox_inches="tight",
        transparent=True,
    )

    fig, ax = plt.subplots(figsize=(4, 3.2), constrained_layout=True)
    fit_j = data.fit_j + data.fit_j.T
    sns.heatmap(
        fit_j[si_data_sort_idx][:, si_data_sort_idx],
        square=True,
        vmin=-1,
        vmax=1,
        center=0,
        cmap="coolwarm",
        # cmap=sns.diverging_palette(220, 20, as_cmap=True),
        # cmap="vlag",
        cbar_kws={"label": r"$J_{ij}$"},
        ax=ax,
    )
    ax.set_xlabel(r"Judge (increasing $\langle s_i \rangle_D$)")
    ax.set_ylabel(r"Judge (increasing $\langle s_i \rangle_D$)")
    ax.set_title(r"Fit pairwise interation parameters $J_{ij}$")
    fig.savefig(
        "results/figures/fit_parameter_j_heatmap.svg",
        bbox_inches="tight",
        transparent=True,
    )

    unique_microstate_obs = np.unique(votes, axis=0)
    empirical_prob = np.empty(unique_microstate_obs.shape[0], dtype=np.float64)
    model_prob = np.empty(unique_microstate_obs.shape[0], dtype=np.float64)
    partition = calculate_partition(data.fit_h, data.fit_j)
    for i, microstate in enumerate(unique_microstate_obs):
        empirical_prob[i] = empirical_probability(microstate, votes)
        model_prob[i] = boltzmann_probability(
            microstate, data.fit_h, data.fit_j, partition
        )

    fig, ax = plt.subplots(figsize=(3, 3), constrained_layout=True)
    ax.scatter(model_prob, empirical_prob, marker="x")
    ax.plot(
        np.linspace(0, 0.3, 10),
        np.linspace(0, 0.3, 10),
        linestyle="dashed",
        linewidth=0.5,
        color="grey",
    )
    # ax.axis("equal")
    ax.set_xlim(0, 0.3)
    ax.set_ylim(0, 0.3)
    ax.set_xlabel(r"$p_\boldsymbol{g}(\boldsymbol{s})$")
    ax.set_ylabel(r"$p_D (\boldsymbol{s})$")
    ax.set_title(
        r"Empirical probability vs. model probability (observed states)", pad=15
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(
        "results/figures/empirical_vs_model_probability.svg",
        bbox_inches="tight",
        transparent=True,
    )

    si_model = model_avg_local_magnetisation(data.fit_h, data.fit_j)
    min_val = min(si_model.min(), si_data.min()) - 0.1
    max_val = max(si_model.max(), si_data.max()) + 0.1
    linear = np.linspace(min_val, max_val, 10)
    fig, ax = plt.subplots(figsize=(3, 3), constrained_layout=True)
    ax.scatter(si_model, si_data, marker="+")
    ax.plot(linear, linear, linestyle="dashed", linewidth=0.5, color="grey")
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.set_xlabel(r"$\langle s_i \rangle$")
    ax.set_ylabel(r"$\langle s_i \rangle_D$")
    ax.set_title("Average local magnetisation: data vs. model", pad=15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(
        "results/figures/avg_local_magnetisation_data_vs_model.svg",
        bbox_inches="tight",
        transparent=True,
    )

    si_sj_model = model_avg_local_correlation(data.fit_h, data.fit_j)
    min_val = 0
    max_val = 1.03
    linear = np.linspace(min_val, max_val, 10)
    fig, ax = plt.subplots(figsize=(3, 3), constrained_layout=True)
    ax.scatter(si_sj_model, si_sj_data, marker="+")
    ax.plot(linear, linear, linestyle="dashed", linewidth=0.5, color="grey")
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.set_xlabel(r"$\langle s_i s_j \rangle$")
    ax.set_ylabel(r"$\langle s_i s_j \rangle_D$")
    ax.set_title("Average local correlation: data vs. model", pad=15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(
        "results/figures/avg_local_correlation_data_vs_model.svg",
        bbox_inches="tight",
        transparent=True,
    )

    conservative_vote_prob = empirical_probability_conservative_vote(votes)
    ks = np.arange(votes.shape[1] + 1, dtype=int)
    k_conservative_votes_prob_independent = np.empty_like(ks, dtype=np.float64)
    for k in ks:
        k_conservative_votes_prob_independent[k] = (
            probability_k_conservative_independent(k, conservative_vote_prob)
        )
    fig, ax = plt.subplots(figsize=(4.5, 2), constrained_layout=True)
    ax.plot(ks, k_conservative_votes_prob_independent, label=r"$P_I (k)$")
    ax.scatter(ks, k_conservative_votes_prob_independent, marker="x")
    ax.set_xticks(ks)
    ax.set_xlabel(r"Number of conservative votes ($k$)")
    ax.set_ylabel(r"$P(k)$")
    ax.set_title(r"$P(k\text{ conservative})$: Independent judges model")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        loc="center left", bbox_to_anchor=(1.01, 0.5), labelspacing=1, frameon=False
    )
    fig.savefig(
        "results/figures/conservative_vote_count_distribution.svg",
        bbox_inches="tight",
        transparent=True,
    )

    k_conservative_votes_prob_data = np.empty_like(ks, dtype=np.float64)
    for k in ks:
        k_conservative_votes_prob_data[k] = probability_k_conservative_data(k, votes)

    fig, ax = plt.subplots(figsize=(4.5, 2), constrained_layout=True)
    ax.plot(ks, k_conservative_votes_prob_independent, label=r"$P_I (k)$")
    ax.scatter(ks, k_conservative_votes_prob_independent, marker="x")
    ax.plot(ks, k_conservative_votes_prob_data, label=r"$P_D (k)$")
    ax.scatter(ks, k_conservative_votes_prob_data, marker="x")
    ax.set_xticks(ks)
    ax.set_xlabel(r"Number of conservative votes ($k$)")
    ax.set_ylabel(r"$P(k)$")
    ax.set_title(r"$P(k\text{ conservative})$: Empirical vs. independent")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        loc="center left", bbox_to_anchor=(1.01, 0.5), labelspacing=1, frameon=False
    )
    fig.savefig(
        "results/figures/conservative_vote_count_distribution_2.svg",
        bbox_inches="tight",
        transparent=True,
    )

    k_conservative_votes_prob_model = np.empty_like(ks, dtype=np.float64)
    for k in ks:
        k_conservative_votes_prob_model[k] = probability_k_conservative_model(
            k, data.fit_h, data.fit_j
        )

    fig, ax = plt.subplots(figsize=(4.5, 2), constrained_layout=True)
    ax.plot(ks, k_conservative_votes_prob_independent, label=r"$P_I (k)$")
    ax.scatter(ks, k_conservative_votes_prob_independent, marker="x")
    ax.plot(ks, k_conservative_votes_prob_data, label=r"$P_D (k)$")
    ax.scatter(ks, k_conservative_votes_prob_data, marker="x")
    ax.plot(ks, k_conservative_votes_prob_model, label=r"$P_P (k)$")
    ax.scatter(ks, k_conservative_votes_prob_model, marker="x")
    ax.set_xticks(ks)
    ax.set_xlabel(r"Number of conservative votes ($k$)")
    ax.set_ylabel(r"$P(k)$")
    ax.set_title(
        r"$P(k\text{ conservative})$: Ising vs. empirical vs. independent", pad=10
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        loc="center left", bbox_to_anchor=(1.01, 0.5), labelspacing=1, frameon=False
    )
    fig.savefig(
        "results/figures/conservative_vote_count_distribution_3.svg",
        bbox_inches="tight",
        transparent=True,
    )

    _, nll = fit_model(votes)
    nll_fit_params = negative_log_likelihood(votes, data.fit_h, data.fit_j, partition)
    fig, ax = plt.subplots(figsize=(4, 2), constrained_layout=True)
    ax.plot(np.arange(len(nll)), nll, label=r"$-\mathcal{L}^{(i)}(\boldsymbol{g})$")
    ax.axhline(
        y=nll_fit_params,
        linestyle="dashed",
        linewidth=0.5,
        color="grey",
        label=r"$-\mathcal{L}(\boldsymbol{g}^\star)$",
    )
    ax.set_xlim(0, None)
    ax.set_xlabel("Optimisation iteration ($i$)")
    ax.set_ylabel(r"$-\mathcal{L}(\boldsymbol{g})$")
    ax.set_title("Negative log-likelihood during optimisation")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend()
    fig.savefig(
        "results/figures/neg_log_likelihood_model_fitting.svg",
        bbox_inches="tight",
        transparent=True,
    )

    metadata = {
        "n_data": int(votes.shape[0]),
        "n_spins": int(votes.shape[1]),
        "n_unique_obs": int(np.unique(votes, axis=0).shape[0]),
    }

    with Path("results/supreme_court.json").open("w") as f:
        json.dump(metadata, f)


if __name__ == "__main__":
    main()
