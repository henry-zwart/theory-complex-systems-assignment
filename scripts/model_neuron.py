from pathlib import Path
import json

from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import numpy.typing as npt

from tcs_assignment import configure_mpl

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


def plot_distribution_between_spikes(
    interspike_times: npt.NDArray[np.float64],
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
    sns.histplot(interspike_times, ax=ax)
    ax.axvline(
        x=tau_0,
        color="red",
        linestyle="dashed",
        linewidth=1,
        label=f"$\\tau_0 = {tau_0:.2f}ms$",
    )
    ax.set_xlabel(r"Inter-spike duration $\tau$ ($ms$)")
    ax.set_ylabel(r"$P(\tau)$")
    ax.set_title("Distribution of durations between neuron spikes")

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.legend(frameon=False)
    return ax


def main():
    """Main function."""
    results = {}

    with NEURON_DATA_PATH.open("r") as f:
        spike_times = np.asarray([float(time.strip()) for time in f.readlines()])
    interspike_times = np.diff(spike_times)

    # Plot inter-spike durations
    fig, ax = plt.subplots(figsize=(4, 2), layout="constrained")
    plot_distribution_between_spikes(interspike_times, ax=ax)
    fig.savefig(Path("results/figures/neuron_interspiking_distribution.pdf"))

    # Calculate refractory period
    tau_0 = determine_refractory(interspike_times)
    results["tau0"] = tau_0

    # Save results
    with Path("results/neuron_modelling.json").open("w") as f:
        json.dump(results, f)


if __name__ == "__main__":
    main()
