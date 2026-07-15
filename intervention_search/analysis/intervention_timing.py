# analysis/intervention_timing.py

"""
Intervention timing analysis (GT→D based)
Uses:
  - obs.diag_triggers  (first detection)
  - interv.intervention_days  (intervention start days)
"""

def get_detection_day(obs):
    """
    First D trigger day.
    """
    return obs.diag_triggers[0] if len(obs.diag_triggers) else None


def get_first_intervention_day(interv):
    """
    First intervention day.
    """
    return interv.intervention_days[0] if len(interv.intervention_days) else None


def compute_intervention_delay(obs, interv):
    """
    Delay = first intervention day - first detection day.
    """
    det = get_detection_day(obs)
    intv = get_first_intervention_day(interv)

    if det is None or intv is None:
        return None

    return intv - det


def intervention_success(obs, interv, max_delay=7):
    """
    Success if intervention starts within max_delay days after detection.
    """
    det = get_detection_day(obs)
    intv = get_first_intervention_day(interv)

    if det is None or intv is None:
        return None

    return 1.0 if intv <= det + max_delay else 0.0


def summarize_intervention(obs, interv, max_delay=7):
    det = get_detection_day(obs)
    first_intv = get_first_intervention_day(interv)
    delay = compute_intervention_delay(obs, interv)
    success = intervention_success(obs, interv, max_delay=max_delay)

    return {
        "detection_day": det,
        "first_intervention_day": first_intv,
        "intervention_delay": delay,
        "intervention_success": success,
    }
