import os
import numpy as np
import matplotlib.pyplot as plt

from plot_result import plot_result
from variant_utils import compute_variant_spread_speed, compute_detection_delay, prevalence_at_detection


def _compute_metrics_for_result(result, variant_intro_day=None):
    """
    result: single run dict from runner.py
    Returns a dict of metrics for that run.
    """

    sim                 = result["sim"]
    diagnoses           = np.array(result["diagnoses"])
    sequences           = np.array(result["sequences"])
    variant_prev        = np.array(result["variant_prevalence"])
    variant_infections  = np.array(result["variant_infections"])
    intervention_log    = result["intervention_log"]
    variant_detected_day = result["variant_detected_day"]

    # Effectiveness
    cum_inf   = sim.results["new_infections"].values.sum()
    peak_inf  = sim.results["new_infections"].values.max()
    time_peak = int(np.argmax(sim.results["new_infections"].values))

    # Practicality
    n_interventions = len(intervention_log)
    if n_interventions > 0:
        durations = [end - start for (start, end, _) in intervention_log]
        avg_duration = float(np.mean(durations))
        avg_intensity = float(np.mean([intensity for (_, _, intensity) in intervention_log]))
    else:
        durations    = []
        avg_duration = 0.0
        avg_intensity = 0.0

    # Surveillance burden
    burden_diagnoses = diagnoses.sum()
    burden_sequences = sequences.sum()
    total_burden     = burden_diagnoses + burden_sequences

    # Variant metrics
    spread_speed = compute_variant_spread_speed(variant_prev) if len(variant_prev) > 0 else 0.0
    detection_delay = compute_detection_delay(variant_intro_day, variant_detected_day) if variant_intro_day is not None else None
    prev_at_det = prevalence_at_detection(variant_prev, variant_detected_day)

    return dict(
        cum_inf=cum_inf,
        peak_inf=peak_inf,
        time_peak=time_peak,
        n_interventions=n_interventions,
        avg_duration=avg_duration,
        avg_intensity=avg_intensity,
        burden_diagnoses=burden_diagnoses,
        burden_sequences=burden_sequences,
        total_burden=total_burden,
        spread_speed=spread_speed,
        detection_delay=detection_delay,
        prev_at_det=prev_at_det,
        variant_detected_day=variant_detected_day,
    )


def _aggregate_strategy_metrics(strategy_name, results, variant_intro_day=None):
    """
    Aggregate metrics across multiple runs of the same strategy.
    """

    per_run = [_compute_metrics_for_result(r, variant_intro_day=variant_intro_day) for r in results]

    def mean_of(key):
        vals = [m[key] for m in per_run if m[key] is not None]
        return float(np.mean(vals)) if vals else None

    return dict(
        name=strategy_name,
        n_runs=len(results),

        # Effectiveness
        mean_cum_inf=mean_of("cum_inf"),
        mean_peak_inf=mean_of("peak_inf"),
        mean_time_peak=mean_of("time_peak"),

        # Practicality
        mean_n_interventions=mean_of("n_interventions"),
        mean_avg_duration=mean_of("avg_duration"),
        mean_avg_intensity=mean_of("avg_intensity"),

        # Surveillance burden
        mean_burden_diagnoses=mean_of("burden_diagnoses"),
        mean_burden_sequences=mean_of("burden_sequences"),
        mean_total_burden=mean_of("total_burden"),

        # Variant metrics
        mean_spread_speed=mean_of("spread_speed"),
        mean_detection_delay=mean_of("detection_delay"),
        mean_prev_at_det=mean_of("prev_at_det"),
    )


def compute_policy_tables(results_dict, variant_intro_day=None):
    """
    results_dict: {
        "baseline": [result_run1, result_run2, ...],
        "strategy_A": [...],
        "strategy_B": [...],
        ...
    }

    Returns:
        metrics_by_strategy: list of aggregated metric dicts.
    """

    metrics_by_strategy = []

    for strat_name, strat_results in results_dict.items():
        agg = _aggregate_strategy_metrics(strat_name, strat_results, variant_intro_day=variant_intro_day)
        metrics_by_strategy.append(agg)

    return metrics_by_strategy


def print_policy_tables(metrics_by_strategy):
    """
    Pretty-print tables to stdout.
    """

    print("\n=== Effectiveness metrics ===")
    print(f"{'Strategy':15s}  {'Cum inf':>10s}  {'Peak inf':>10s}  {'Time to peak':>13s}")
    for m in metrics_by_strategy:
        print(f"{m['name']:15s}  {m['mean_cum_inf']:10.1f}  {m['mean_peak_inf']:10.1f}  {m['mean_time_peak']:13.1f}")

    print("\n=== Practicality metrics ===")
    print(f"{'Strategy':15s}  {'Interventions':>13s}  {'Avg dur':>10s}  {'Avg intensity':>14s}")
    for m in metrics_by_strategy:
        print(f"{m['name']:15s}  {m['mean_n_interventions']:13.1f}  {m['mean_avg_duration']:10.1f}  {m['mean_avg_intensity']:14.2f}")

    print("\n=== Surveillance burden ===")
    print(f"{'Strategy':15s}  {'Diag':>10s}  {'Seq':>10s}  {'Total':>10s}")
    for m in metrics_by_strategy:
        print(f"{m['name']:15s}  {m['mean_burden_diagnoses']:10.1f}  {m['mean_burden_sequences']:10.1f}  {m['mean_total_burden']:10.1f}")

    print("\n=== Variant detection summary ===")
    print(f"{'Strategy':15s}  {'Spread speed':>13s}  {'Detect delay':>13s}  {'Prev@det':>10s}")
    for m in metrics_by_strategy:
        dd = m['mean_detection_delay']
        pa = m['mean_prev_at_det']
        print(f"{m['name']:15s}  {m['mean_spread_speed']:13.3f}  {dd if dd is not None else float('nan'):13.1f}  {pa if pa is not None else float('nan'):10.3f}")


def generate_policy_report(results_dict, outdir="policy_report", variant_intro_day=None):
    """
    High-level entry point.

    - Computes metrics per strategy
    - Prints tables
    - Saves plots per run
    """

    os.makedirs(outdir, exist_ok=True)

    # 1) Compute aggregated metrics
    metrics_by_strategy = compute_policy_tables(results_dict, variant_intro_day=variant_intro_day)

    # 2) Print tables to console
    print("\n================ POLICY REPORT ================")
    print_policy_tables(metrics_by_strategy)
    print("===============================================")

    # 3) Save plots for each run
    for strat_name, strat_results in results_dict.items():
        for i, result in enumerate(strat_results):
            fname = os.path.join(outdir, f"{strat_name}_run{i+1}.png")
            plot_result(result, save_path=fname)

    # 4) Optional: summary scatter of practicality vs effectiveness
    fig, ax = plt.subplots(figsize=(7, 5))
    for m in metrics_by_strategy:
        ax.scatter(m["mean_n_interventions"], m["mean_cum_inf"], label=m["name"], s=80)
    ax.set_xlabel("Mean intervention count (practicality)")
    ax.set_ylabel("Mean cumulative infections (effectiveness)")
    ax.set_title("Practicality vs effectiveness across strategies")
    ax.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "practicality_vs_effectiveness.png"), dpi=200)
    plt.close(fig)

    return metrics_by_strategy
