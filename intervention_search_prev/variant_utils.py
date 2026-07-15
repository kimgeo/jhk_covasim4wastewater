import numpy as np

# 1) Haplotype reconstruction (placeholder)

def reconstruct_haplotype(genome):
    """
    Placeholder haplotype reconstruction.
    In real use: apply actual reconstruction logic.
    """
    return genome

# 2) Variant detection

def detect_variant(haplotype, reference):
    """
    Return True if haplotype differs from reference.
    """
    return haplotype != reference

# 3) Variant prevalence

def compute_variant_prevalence(sim, t, variant_index=1):
    """
    Compute variant prevalence at day t.
    variant_index: index of variant in Covasim results (default = 1)
    """

    # If variant module not enabled
    if "variant" not in sim.results:
        return 0.0

    total_inf = sim.results["new_infections"][t]

    # Avoid division by zero
    if total_inf <= 0:
        return 0.0

    variant_inf = sim.results["variant"]["new_infections_by_variant"][variant_index, t]

    return float(variant_inf / total_inf)

# 4) Variant spread speed

def compute_variant_spread_speed(variant_prev, window=7):
    """
    Compute variant spread speed as slope over last 'window' days.
    variant_prev: array of daily variant prevalence
    """

    if len(variant_prev) < window + 1:
        return 0.0

    # Simple slope: prev[t] - prev[t-window]
    return variant_prev[-1] - variant_prev[-window]


# 5) Variant detection delay

def compute_detection_delay(intro_day, detected_day):
    """
    Compute delay between variant introduction and detection.
    """

    if detected_day is None:
        return None

    return detected_day - intro_day


# 6) Variant prevalence at detection

def prevalence_at_detection(variant_prev, detected_day):
    """
    Return variant prevalence at detection day.
    """

    if detected_day is None:
        return None

    if detected_day >= len(variant_prev):
        return None

    return variant_prev[detected_day]
