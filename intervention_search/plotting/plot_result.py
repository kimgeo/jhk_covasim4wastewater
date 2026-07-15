# plotting/plot_result.py

import matplotlib.pyplot as plt
import numpy as np

def plot_single_run(sim, results, savepath=None, show=True):
    """
    Plot:
      - new infections (blue)
      - daily diagnoses (orange)
      - daily sequences (green)
      - GT triggers (red dots)
      - D triggers (blue dots)
      - intervention days (purple dashed lines)
      - detection_day (first D trigger)
    """

    T = sim.npts
    days = np.arange(T)

    new_inf = sim.results["new_infections"]
    diag = results.diagnoses
    seqs = results.sequences

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(days, new_inf, label="New infections", color="tab:blue")
    ax.plot(days, diag, label="Diagnoses", color="tab:orange")
    ax.plot(days, seqs, label="Sequences", color="tab:green")

    # GT triggers (red dots)
    if len(results.gt_triggers) > 0:
        ax.scatter(results.gt_triggers,
                   new_inf[results.gt_triggers],
                   color="red", s=40, label="GT triggers")

    # D triggers (blue dots)
    if len(results.diag_triggers) > 0:
        ax.scatter(results.diag_triggers,
                   diag[results.diag_triggers],
                   color="blue", s=40, label="D triggers")

    # detection_day (first D trigger)
    if results.detection_day is not None:
        ax.axvline(results.detection_day,
                   color="red", linestyle="--",
                   label=f"Detection day = {results.detection_day}")

    # intervention days (multiple)
    for d in results.intervention_days:
        ax.axvline(d, color="purple", linestyle="--", alpha=0.7)

    ax.set_xlabel("Day")
    ax.set_ylabel("Count")
    ax.set_title("Single run: infections, diagnoses, sequences")

    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=200)

    if show:
        plt.show()
    else:
        plt.close(fig)
