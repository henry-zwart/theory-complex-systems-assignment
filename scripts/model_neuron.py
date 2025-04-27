"""Experiment code for Q1. Modelling the activity of a single neuron."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import seaborn as sns
from matplotlib.axes import Axes
from scipy.stats import linregress

from tcs_assignment import configure_mpl

RANDOM_SEED = 42

configure_mpl()


NEURON_DATA_PATH = Path("data/neuron.txt")


def determine_refractory(interspike_times: npt.NDArray[np.float64]) -> float:
    """Determines the refractory period from an array on neuron inter-spike times.

    The refractory period is determined to be the minimum possible duration
    between spikes. It is bounded above by the minimum observed duration.
    In reality, the duration of the refractory period is itself probabilistic,
    but with sufficient data it is approximated well by the minimum observed
    inter-spike time.

    Args:
        interspike_times: 1D numpy array of durations between consecutive spikes,
            ordered by event time.

    Returns:
        Estimate for the refractory period duration, calculated as the minimum
        observed duration between spikes.
    """
    return interspike_times.min()


def fit_exponential_function(interspike_times: npt.NDArray[np.float64]):
    """Calculate decay parameter for exponential function.

    Args:
        interspike_times: 1D numpy array of durations between consecutive spikes,
            ordered by event time.

    Returns:
        (lambda, P(tau_0)), calculated from linear least squares
        (rvalue, pvalue, stderr)
    """
    # lmbda_inv = (interspike_times - determine_refractory(interspike_times)).mean()
    # lmbda = 1 / lmbda_inv
    # return lmbda

    tau_0 = determine_refractory(interspike_times)
    interspike_times = interspike_times - tau_0
    prob_tau, tau_bins = np.histogram(interspike_times, bins=50, density=True)

    # Calculate each tau as the center of each bin
    binwidth = tau_bins[1] - tau_bins[0]
    taus = tau_bins[:-1] + binwidth / 2

    # Disregard taus with zero-probability
    nonzero = prob_tau > 0
    taus = taus[nonzero]
    prob_tau = prob_tau[nonzero]

    # Inputs for linear regression: tau - tau_0, and -log prob (tau - tau_0)
    neg_log_prob_tau = -np.log(prob_tau)

    # Fit linear model
    result = linregress(taus, neg_log_prob_tau)

    lmbda = result.slope
    prob_tau_0 = np.exp(-result.intercept)
    return (
        (lmbda, prob_tau_0),
        (taus, neg_log_prob_tau),
        (result.rvalue, result.pvalue, result.stderr),
    )


def plot_distribution_between_spikes(
    interspike_times: npt.NDArray[np.float64],
    alpha: float = 1.0,
    label: str | None = None,
    show_tau_0: bool = True,
    ax: Axes | None = None,
) -> Axes:
    """Plot histogram of durations between consecutive neuron spikes.

    Args:
        interspike_times: 1D numpy array of durations between consecutive spikes,
            ordered by event time.
        ax: Optional Matplotlib Axes object. If not provided, uses the current artist.

    Returns:
        Matplotlib Axes object with histogram of inter-spike duration times.
    """
    if ax is None:
        ax = plt.gca()

    tau_0 = determine_refractory(interspike_times)
    sns.histplot(interspike_times, ax=ax, stat="density", alpha=alpha, label=label)
    if show_tau_0:
        ax.axvline(
            x=tau_0,
            color="red",
            linestyle="dashed",
            linewidth=1,
            label=f"$\\tau_0 = {tau_0:.2f}$ms",
        )

    ax.set_xlabel(r"Inter-spike duration $\tau$ (ms)")
    ax.set_ylabel(r"$P(\tau)$")
    ax.set_title("Distribution of durations between neuron spikes")

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.legend(frameon=False)
    return ax


def plot_linear_fit(x, y, theoretical_slope, theoretical_int, ax: Axes | None = None):
    if ax is None:
        ax = plt.gca()

    # Empirical results
    ax.plot(x, y, linewidth=1)

    # Theoretical results
    x_theoretical = np.linspace(x.min(), x.max(), 1000)
    y_theoretical = theoretical_int + theoretical_slope * x_theoretical

    ax.plot(x_theoretical, y_theoretical, linestyle="dashed", linewidth=1)

    ax.set_xlim(0, None)
    ax.set_ylim(0, None)

    ax.set_xlabel(r"$\tau - \tau_0$")
    ax.set_ylabel(r"$-\log{P(\tau - \tau_0)}$")
    ax.set_title("Log-transformed distribution of inter-spike durations")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return ax


def sample_delayed_exponential(
    lmbda: float, tau_0: float, n: int, rng: np.random.Generator
) -> npt.NDArray[np.float64]:
    """Generate random samples from a delayed exponential."""
    u = rng.random(size=n)
    return -(1 / lmbda) * np.log(1 - u) + tau_0


def main():
    """Main function."""
    results = {}

    with NEURON_DATA_PATH.open("r") as f:
        spike_times = np.asarray([float(time.strip()) for time in f.readlines()])
    interspike_times = np.diff(spike_times)

    # Plot inter-spike durations
    fig, ax = plt.subplots(figsize=(4, 2), layout="constrained")
    plot_distribution_between_spikes(interspike_times, ax=ax)
    fig.savefig(
        Path("results/figures/neuron_interspiking_distribution.svg"),
        bbox_inches="tight",
        transparent=True,
    )

    # Calculate refractory period
    tau_0 = determine_refractory(interspike_times)
    (lmbda, p_tau_0), (shifted_taus, neg_log_prob_shifted_taus), (rvalue, pvalue, _) = (
        fit_exponential_function(interspike_times)
    )
    fig, ax = plt.subplots(figsize=(4, 2), layout="constrained")
    plot_linear_fit(
        shifted_taus,
        neg_log_prob_shifted_taus,
        theoretical_slope=lmbda,
        theoretical_int=-np.log(p_tau_0),
        ax=ax,
    )
    fig.savefig(
        Path("results/figures/neuron_exponential_fit.svg"),
        bbox_inches="tight",
        transparent=True,
    )

    # Generate new datapoints
    rng = np.random.default_rng(seed=RANDOM_SEED)
    interspike_samples = sample_delayed_exponential(lmbda, tau_0, n=1000, rng=rng)
    fig, ax = plt.subplots(figsize=(4, 2), constrained_layout=True)
    plot_distribution_between_spikes(
        interspike_times, alpha=0.5, label="Observed", ax=ax
    )
    plot_distribution_between_spikes(
        interspike_samples, alpha=0.3, label="Samples", show_tau_0=False, ax=ax
    )

    fig.savefig(
        "results/figures/neuron_compare_empirical_theoretical.svg",
        bbox_inches="tight",
        transparent=True,
    )

    results["n"] = len(interspike_times) - 1
    results["tau0"] = tau_0
    results["mean_tau"] = (interspike_times).mean()
    results["mean_spike_rate"] = 1 / interspike_times.mean()
    results["exponential_fit"] = {
        "lambda": lmbda,
        "p_tau_0": p_tau_0,
        "r_value": rvalue,
        "p_value": pvalue,
    }

    # Save results
    with Path("results/neuron_modelling.json").open("w") as f:
        json.dump(results, f)


if __name__ == "__main__":
    main()
