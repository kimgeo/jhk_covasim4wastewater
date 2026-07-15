"""
GT -> D based detection delay analysis (ObservationModel-based)
"""

def get_first_detection_day(obs):
    """
    First day where D trigger occurred.
    """
    return obs.diag_triggers[0] if len(obs.diag_triggers) else None


def compute_detection_delay(obs):
    """
    Detection delay = first D trigger day.
    (No emergence baseline.)
    """
    return get_first_detection_day(obs)


def detection_success(obs, window=14):
    """
    Detection success = first D trigger occurs within `window` days.
    """
    det = get_first_detection_day(obs)
    if det is None:
        return 0.0
    return 1.0 if det <= window else 0.0


def summarize_detection(sim, obs, window=14):
    """
    Summary dict for detection metrics.
    """
    det = get_first_detection_day(obs)
    delay = compute_detection_delay(obs)
    success = detection_success(obs, window=window)

    return {
        "first_detection_day": det,
        "detection_delay": delay,
        "detection_success_window": success,
    }
