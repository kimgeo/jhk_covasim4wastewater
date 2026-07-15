# analysis/outbreak_impact.py

"""
Outbreak impact metrics:
- peak infections
- total infections
- outbreak duration
- infections after intervention
"""

import numpy as np

def peak_infections(sim):
    return int(np.max(sim.results["new_infections"]))


def total_infections(sim):
    return int(sim.results["cum_infections"][-1])


def outbreak_duration(sim):
    """
    Duration until new infections drop to zero.
    """
    inf = sim.results["new_infections"]
    for t in range(len(inf) - 1, -1, -1):
        if inf[t] > 0:
            return t
    return 0


def infections_after_intervention(sim, interv):
    if interv.intervention_day is None:
        return None
    return int(np.sum(sim.results["new_infections"][interv.intervention_day:]))


def summarize_outbreak(sim, interv):
    return {
        "peak_infections": peak_infections(sim),
        "total_infections": total_infections(sim),
        "outbreak_duration": outbreak_duration(sim),
        "infections_after_intervention": infections_after_intervention(sim, interv),
    }
