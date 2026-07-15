# analysis/effectiveness.py

"""
Effectiveness vs practicality metrics
"""

import numpy as np

def practicality(results):
    """
    Total number of sequences generated.
    """
    return int(np.sum(results.sequences))


def effectiveness(results):
    """
    Effectiveness = inverse of avg_resolution_time
    Smaller resolution_time → faster detection → higher effectiveness
    """
    if results.avg_resolution_time is None:
        return 0.0

    return 1.0 / (1.0 + results.avg_resolution_time)


def policy_point(results):
    """
    Return a dict for scatter plot:
      x = practicality
      y = effectiveness
    """
    return {
        "practicality": practicality(results),
        "effectiveness": effectiveness(results),
    }
