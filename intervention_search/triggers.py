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
