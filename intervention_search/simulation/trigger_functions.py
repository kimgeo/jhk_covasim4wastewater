# intervention_search/simulation/trigger_functions.py

import numpy as np
from typing import Sequence

def _to_array(x: Sequence[float]) -> np.ndarray:
    return np.asarray(x, dtype=float)


# Relative threshold

def relative_threshold(
    diagnoses: Sequence[float],
    t: int,
    proportion: float,
    pop_size: int,
    window: int = 1,
) -> bool:
    arr = _to_array(diagnoses)

    if t < window - 1:
        return False

    avg = arr[t-window+1:t+1].mean()
    rel = avg / pop_size

    return rel >= proportion


# Weekly growth

def weekly_growth(
    diagnoses: Sequence[float],
    t: int,
    ratio: float,
    window: int = 7,
) -> bool:
    arr = _to_array(diagnoses)

    if t < 2 * window - 1:
        return False

    curr = arr[t-window+1:t+1].mean()
    prev = arr[t-2*window+1:t-window+1].mean()

    if prev <= 0:
        return False

    return (curr / prev) >= ratio


# Slope trigger

def slope_trigger(
    diagnoses: Sequence[float],
    t: int,
    slope_threshold: float,
) -> bool:
    arr = _to_array(diagnoses)

    if t < 1:
        return False

    slope = arr[t] - arr[t-1]
    return slope >= slope_threshold


# Sustained increase

def sustained_increase(
    diagnoses: Sequence[float],
    t: int,
    days: int = 5,
) -> bool:
    arr = _to_array(diagnoses)

    if t < days:
        return False

    for i in range(days):
        if arr[t-i] <= arr[t-i-1]:
            return False

    return True


# Unified dispatcher

def check_triggers(
    t: int,
    diagnoses: Sequence[float],
    sequences: Sequence[float],
    variant_prev: Sequence[float],
    trigger_list,
    sim=None,
):
    """
    trigger_list = [
        ("weekly_growth", {"ratio": 2.0}),
        ("slope_trigger", {"slope_threshold": 5}),
        ("sustained_increase", {"days": 5}),
        ("relative_threshold", {"proportion": 0.0001, "window": 3}),
    ]
    """

    for trigger_name, pars in trigger_list:

        if trigger_name == "weekly_growth":
            if weekly_growth(diagnoses, t, pars["ratio"], pars.get("window", 7)):
                return True

        elif trigger_name == "slope_trigger":
            if slope_trigger(diagnoses, t, pars["slope_threshold"]):
                return True

        elif trigger_name == "sustained_increase":
            if sustained_increase(diagnoses, t, pars.get("days", 5)):
                return True

        elif trigger_name == "relative_threshold":
            if sim is None:
                raise ValueError("relative_threshold requires sim")
            pop = sim.pars["pop_size"]
            if relative_threshold(
                diagnoses,
                t,
                pars["proportion"],
                pop,
                pars.get("window", 1)
            ):
                return True

        else:
            raise ValueError(f"Unknown trigger: {trigger_name}")

    return False
