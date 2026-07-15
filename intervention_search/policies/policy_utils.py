# policies/policy_utils.py

import numpy as np

def compute_detection_delay(sim, obs):
    """
    Returns detection delay:
      detection_day - true_emergence_day
    """
    true_day = sim.true_emergence_day if hasattr(sim, "true_emergence_day") else None
    if true_day is None or obs.detection_day is None:
        return None
    return obs.detection_day - true_day


def compute_sensitivity(sim, obs, window=14):
    """
    Sensitivity = probability of detecting variant within 'window' days
    after true emergence.
    """
    true_day = getattr(sim, "true_emergence_day", None)
    if true_day is None:
        return None

    if obs.detection_day is None:
        return 0.0

    return 1.0 if obs.detection_day <= true_day + window else 0.0


def compute_sequencing_burden(obs):
    """
    Total number of sequences generated.
    """
    return int(np.sum(obs.daily_sequences))


def summarize_policy(sim, obs):
    """
    Return a dict summarizing policy performance.
    """
    return {
        "detection_day": obs.detection_day,
        "sequencing_burden": compute_sequencing_burden(obs),
        "sensitivity_14d": compute_sensitivity(sim, obs, window=14),
    }
