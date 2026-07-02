"""
evaluation.py

Expanded evaluation module:
- Effectiveness metrics
- Practicality metrics
- Surveillance cost metrics
- Variant metrics
"""

import numpy as np

def evaluate_results(results):
    """
    Compute evaluation metrics from simulation + observation + intervention tracking.

    Parameters
    ----------
    results : dict
        {
            "diagnoses": array,
            "sequences": array,
            "trigger_events": list,
            "sim": sim,
            "intervention_log": list of (start, end, intensity),
            "variant_info": dict {
                "variant_infections": array,
                "variant_detected_day": int or None,
                "variant_dominance_day": int or None
            },
            "surveillance_pars": dict {
                "p_test": float,
                "p_seq": float,
                "delay_pmf": array
            }
        }

    Returns
    -------
    metrics : dict
    """

    sim = results["sim"]
    diagnoses = results["diagnoses"]
    sequences = results["sequences"]
    triggers = results["trigger_events"]

    intervention_log = results.get("intervention_log", [])
    variant_info = results.get("variant_info", {})
    surv = results.get("surveillance_pars", {})

    # 1) effectiveness metrics 

    inf = sim.results["new_infections"].values
    total_infections = inf.sum()
    peak_infections = inf.max()
    peak_day = int(np.argmax(inf))

    # outbreak duration = last day of infection happening
    outbreak_duration = int(np.where(inf > 0)[0][-1]) if inf.sum() > 0 else 0

    # Variant metrics
    variant_infections = variant_info.get("variant_infections", None)
    if variant_infections is not None:
        total_variant_infections = variant_infections.sum()
        variant_dominance_day = variant_info.get("variant_dominance_day", None) # The day the dominant variant got detected for the first time 
        variant_detected_day = variant_info.get("variant_detected_day", None)
    else:
        total_variant_infections = None
        variant_dominance_day = None
        variant_detected_day = None

    # 2) practicality metrics 
    # intervention_log: list of (start_day, end_day, intensity)
    num_interventions = len(intervention_log) # how many intervention episodes

    if num_interventions > 0:
        durations = [end - start for start, end, _ in intervention_log] # length of each intervention
        total_intervention_duration = sum(durations)
        mean_intervention_duration = np.mean(durations)

        intensities = [intensity for _, _, intensity in intervention_log] # how intense the intervention was (ex. beta reduction factor)
        mean_intervention_intensity = np.mean(intensities)

        # intervention interval
        starts = [start for start, _, _ in intervention_log]
        starts_sorted = sorted(starts)
        if len(starts_sorted) > 1:
            intervals = np.diff(starts_sorted)
            mean_intervention_interval = np.mean(intervals)
        else:
            mean_intervention_interval = None
    else:
        total_intervention_duration = 0
        mean_intervention_duration = 0
        mean_intervention_intensity = 0
        mean_intervention_interval = None

    # feasibility score
    # if an intervention is applied for too long, it would be feasible
    T = sim.npts
    if total_intervention_duration > 0.4 * T: # low feasibility if interventions were held during 40% of the simulation duration
        feasibility_score = "low"
    elif num_interventions > 20: # low feasibility if interventions were held more than 20 times
        feasibility_score = "low"
    elif total_intervention_duration > 0.2 * T:
        feasibility_score = "medium"
    else:
        feasibility_score = "high"

    # 3) surveillance cost metrics

    p_test = surv.get("p_test", None)
    p_seq = surv.get("p_seq", None)
    delay_pmf = surv.get("delay_pmf", None)

    if delay_pmf is not None: # calculate the mean detection delay
        mean_detection_delay = np.sum(np.arange(len(delay_pmf)) * delay_pmf)
    else:
        mean_detection_delay = None

    # burden category
    if p_test is not None and p_seq is not None:
        if p_test >= 0.3 or p_seq >= 0.2: # high p_test and p_seq -> high burden
            surveillance_burden = "high"
        elif p_test >= 0.15 or p_seq >= 0.1:
            surveillance_burden = "medium"
        else:
            surveillance_burden = "low"
    else:
        surveillance_burden = None

    # 4) trigger metrics 
    
    num_triggers = len(triggers)
    trigger_days = [day for day, name in triggers]

    if len(trigger_days) > 1:
        trigger_intervals = np.diff(sorted(trigger_days))
        mean_trigger_interval = np.mean(trigger_intervals)
    else:
        mean_trigger_interval = None

    # return all metrics

    metrics = {
        # Effectiveness
        "total_infections": total_infections,
        "peak_infections": peak_infections,
        "peak_day": peak_day,
        "outbreak_duration": outbreak_duration,

        # Variant
        "total_variant_infections": total_variant_infections,
        "variant_detected_day": variant_detected_day,
        "variant_dominance_day": variant_dominance_day,

        # Practicality
        "num_interventions": num_interventions,
        "total_intervention_duration": total_intervention_duration,
        "mean_intervention_duration": mean_intervention_duration,
        "mean_intervention_intensity": mean_intervention_intensity,
        "mean_intervention_interval": mean_intervention_interval,
        "feasibility_score": feasibility_score,

        # Surveillance
        "diagnosis_rate": p_test,
        "sequencing_rate": p_seq,
        "mean_detection_delay": mean_detection_delay,
        "surveillance_burden": surveillance_burden,

        # Trigger
        "num_triggers": num_triggers,
        "trigger_days": trigger_days,
        "mean_trigger_interval": mean_trigger_interval,
    }

    return metrics
