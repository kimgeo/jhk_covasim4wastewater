"""
evaluation.py

Evaluate the outcome of a policy simulation.
Takes the output of runner.run_simulation() and computes summary metrics.
"""

import numpy as np

def evaluate_results(results):
    """
    Compute evaluation metrics from simulation results.

    Parameters
    ----------
    results : dict
        Output from run_simulation():
        {
            "diagnoses": array,
            "sequences": array,
            "trigger_events": list,
            "sim": sim
        }

    Returns
    -------
    metrics : dict
        Dictionary of evaluation metrics.
    """

    sim = results["sim"]
    diagnoses = results["diagnoses"]
    sequences = results["sequences"]
    triggers = results["trigger_events"]

    # Basic epidemic metrics
    inf = sim.results["new_infections"].values

    total_infections = inf.sum()
    peak_infections = inf.max()
    peak_day = int(np.argmax(inf))
    
    # Observation metrics
    total_diagnoses = diagnoses.sum()
    total_sequences = sequences.sum()
    peak_diagnoses = diagnoses.max()

    # Trigger metrics
    num_triggers = len(triggers)
    trigger_days = [day for day, name in triggers]

    metrics = {
        "total_infections": total_infections,
        "peak_infections": peak_infections,
        "peak_day": peak_day,
        "total_diagnoses": total_diagnoses,
        "total_sequences": total_sequences,
        "peak_diagnoses": peak_diagnoses,
        "num_triggers": num_triggers,
        "trigger_days": trigger_days,
    }

    return metrics
