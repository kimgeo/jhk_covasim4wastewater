# plotting/plot_variant.py

"""
Variant prevalence plot:
- population-wide variant prevalence
- sequenced-only variant prevalence (optional overlay)
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_variant_prevalence(sim, variant_metrics, use_seq=False, savepath=None, show=True):
    """
    Plot variant prevalence over time.

    variant_metrics:
        - variant_names
        - prevalence[v]: list over time (population-based)
        - seq_prevalence[v]: list over time (sequenced-only, if use_seq=True)
    """

    T = sim.npts
    days = np.arange(T)

    fig, ax = plt.subplots(figsize=(10, 6))

    for v in variant_metrics.variant_names:
        prev = np.array(variant_metrics.prevalence[v])
        ax.plot(days, prev, label=f"{v} (pop)", linewidth=2)

        if use_seq:
            seq_prev = np.array(variant_metrics.seq_prevalence[v])
            ax.plot(
                days,
                seq_prev,
                linestyle="--",
                linewidth=1.5,
                label=f"{v} (seq)",
            )

    ax.set_xlabel("Day")
    ax.set_ylabel("Prevalence")
    ax.set_title("Variant prevalence over time")

    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=200)

    if show:
        plt.show()
    else:
        plt.close(fig)
