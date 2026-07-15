# plotting/plot_policy_space.py

"""
Policy space scatter plot:
- x-axis: practicality (sequencing burden)
- y-axis: effectiveness (detection success / sensitivity)
- Color : detection_delay (optional)
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_policy_space(points, color_key=None, savepath=None, show=True):
    """
    points: list of dicts, each with at least:
        - "practicality"
        - "effectiveness"
    optionally:
        - color_key (e.g. "detection_delay")
        - "label" for annotation
    """

    x = np.array([p["practicality"] for p in points])
    y = np.array([p["effectiveness"] for p in points])

    if color_key is not None:
        c = np.array([p.get(color_key, np.nan) for p in points])
    else:
        c = "tab:blue"

    fig, ax = plt.subplots(figsize=(8, 6))

    sc = ax.scatter(x, y, c=c, cmap="viridis", s=60, edgecolor="k", alpha=0.8)

    ax.set_xlabel("Practicality (total sequences)")
    ax.set_ylabel("Effectiveness (detection success)")
    ax.set_title("Policy space: practicality vs effectiveness")

    if color_key is not None:
        cb = plt.colorbar(sc, ax=ax)
        cb.set_label(color_key)

    for p in points:
        if "label" in p:
            ax.annotate(
                p["label"],
                (p["practicality"], p["effectiveness"]),
                textcoords="offset points",
                xytext=(3, 3),
                fontsize=8,
                alpha=0.7,
            )

    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=200)

    if show:
        plt.show()
    else:
        plt.close(fig)
