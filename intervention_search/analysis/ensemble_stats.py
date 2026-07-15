# analysis/ensemble_stats.py

"""
Ensemble statistics for multiple runs.
Useful for boxplots, confidence intervals, etc.
"""

import numpy as np

def summarize_ensemble(runs, key):
    """
    runs: list of dicts
    key: metric name, e.g. "avg_resolution_time"
    """
    values = [r[key] for r in runs if r[key] is not None]

    if len(values) == 0:
        return None

    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "std": float(np.std(values)),
        "q25": float(np.percentile(values, 25)),
        "q75": float(np.percentile(values, 75)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "n": len(values),
    }
