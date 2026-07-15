# Trigger functions for diagnosis-based intervention policies.
# Each trigger function takes:
# - diagnoses: 1D array-like of daily observed diagnoses (output of observe())
# - t:         integer index for the current day (0-based)
# - other parameters specific to the trigger
# and returns:
# - True  if the trigger condition is satisfied on day t
# - False otherwise

from typing import Sequence
import numpy as np

def _to_array(x: Sequence[float]) -> np.ndarray: 
    # Convert input to a NumPy float array.
    return np.asarray(x, dtype=float)

def relative_threshold(
    diagnoses,
    t,
    proportion,
    sim,
    window=1,
):
    arr = np.asarray(diagnoses, dtype=float)

    # Not enough data for window
    if t < window - 1:
        return False

    # Average diagnoses over the last 'window' days
    avg = arr[t-window+1:t+1].mean()

    # Population size
    pop = sim.pars["pop_size"]

    # Relative incidence
    rel = avg / pop

    return rel >= proportion


def weekly_growth(
    diagnoses: Sequence[float],
    t: int,
    ratio: float,
    window: int = 7,
) -> bool:
    #Trigger when this week's avg diagnoses / last week's avg >= ratio.

    arr = _to_array(diagnoses)

    if t < 2 * window - 1:
        return False

    curr = arr[t-window+1:t+1].mean()
    prev = arr[t-2*window+1:t-window+1].mean()

    if prev <= 0:
        return False

    return (curr / prev) >= ratio


def slope_trigger(
    diagnoses: Sequence[float],
    t: int,
    slope_threshold: float,
) -> bool:
    # Trigger when diagnoses[t] - diagnoses[t-1] >= slope_threshold.

    arr = _to_array(diagnoses)

    if t < 1:
        return False

    slope = arr[t] - arr[t-1]
    return slope >= slope_threshold


def sustained_increase(
    diagnoses: Sequence[float],
    t: int,
    days: int = 5,
) -> bool:
    # Trigger when diagnoses increase for days consecutive days.
    arr = _to_array(diagnoses)

    if t < days:
        return False

    for i in range(days):
        if arr[t-i] <= arr[t-i-1]:
            return False

    return True

def variant_prevalence_threshold(
    variant_prev: Sequence[float],
    t: int,
    threshold: float,
) -> bool:
    """
    Trigger when variant prevalence exceeds threshold.
    variant_prev[t] is between 0 and 1.
    """
    arr = _to_array(variant_prev)
    if t < 0:
        return False
    return arr[t] >= threshold


def variant_growth_rate_trigger(
    variant_prev: Sequence[float],
    t: int,
    growth_threshold: float,
    window: int = 7,
) -> bool:
    """
    Trigger when variant prevalence increases rapidly:
    (prev[t] - prev[t-window]) >= growth_threshold
    """
    arr = _to_array(variant_prev)

    if t < window:
        return False

    growth = arr[t] - arr[t-window]
    return growth >= growth_threshold

def check_triggers(
    t,
    diagnoses,
    sequences,
    variant_prev,
    trigger_list,
):
    """
    Unified trigger dispatcher.

    runner.py expects:
        check_triggers(t, diagnoses, sequences, variant_prev, trigger_list)

    trigger_list is a list of:
        (trigger_name, trigger_pars)
    """

    for trigger_name, pars in trigger_list:

        # 1) sustained increase
        if trigger_name == "sustained_increase":
            if sustained_increase(diagnoses, t, pars.get("days", 5)):
                return True

        # 2) slope trigger
        elif trigger_name == "slope_trigger":
            if slope_trigger(diagnoses, t, pars["slope_threshold"]):
                return True

        # 3) weekly growth
        elif trigger_name == "weekly_growth":
            if weekly_growth(diagnoses, t, pars["ratio"], pars.get("window", 7)):
                return True

        # 4) relative threshold
        elif trigger_name == "relative_threshold":
            # sim is not passed here → runner.py must not use this trigger
            # or you can modify runner.py to pass sim
            raise NotImplementedError("relative_threshold requires sim, runner.py does not pass sim")

        # 5) variant prevalence threshold
        elif trigger_name == "variant_prevalence_threshold":
            if variant_prevalence_threshold(variant_prev, t, pars["threshold"]):
                return True

        # 6) variant growth rate trigger
        elif trigger_name == "variant_growth_rate_trigger":
            if variant_growth_rate_trigger(
                variant_prev,
                t,
                pars["growth_threshold"],
                pars.get("window", 7)
            ):
                return True

        # 7) sequencing-based trigger (optional)
        elif trigger_name == "sequence_threshold":
            if sequences[t] >= pars["threshold"]:
                return True

        else:
            raise ValueError(f"Unknown trigger: {trigger_name}")

    return False
