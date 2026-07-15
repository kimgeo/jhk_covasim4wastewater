# plotting/plot_ensemble.py

import matplotlib.pyplot as plt
import numpy as np

def boxplot_metric(runs, key, ax=None, title=None):
    values = [r[key] for r in runs if r.get(key) is not None]

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = ax.figure

    if len(values) == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        return fig, ax

    ax.boxplot(values, vert=True)
    ax.set_ylabel(key)
    ax.set_title(title or f"Boxplot of {key}")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig, ax


def plot_ensemble_boxplots(runs, keys, savepath=None, show=True):
    n = len(keys)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4))

    if n == 1:
        axes = [axes]

    for ax, key in zip(axes, keys):
        boxplot_metric(runs, key, ax=ax, title=key)

    if savepath is not None:
        plt.savefig(savepath, dpi=200)

    if show:
        plt.show()
    else:
        plt.close(fig)
