import numpy as np
import matplotlib.pyplot as plt

def plot_top3_and_baselines(
    baseline_sim,
    policy_results_list,
    best_total_inf,
    best_total_deaths,
    best_shortest_intervention,
    ax=None,
    title="Baseline, Mean, and Top-3 Best Runs (Cumulative Infections)"
):
    # Plot cumulative infections for:
    # - baseline (no intervention)
    # - policy mean of all runs during the search
    # - best run by total infections
    # - best run by total deaths
    # - best run by shortest intervention duration (least cost)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.figure

    # Time axis
    T = baseline_sim.npts
    days = np.arange(T)

    # Baseline cumulative infections
    baseline_cum = np.cumsum(baseline_sim.results["new_infections"])

    # Policy mean cumulative infections
    policy_inf_all = np.array([
        res["sim"].results["new_infections"] for res in policy_results_list
    ])
    policy_mean = policy_inf_all.mean(axis=0)
    policy_mean_cum = np.cumsum(policy_mean)

    # Best runs cumulative infections
    best_inf_cum = np.cumsum(best_total_inf["sim"].results["new_infections"])
    best_deaths_cum = np.cumsum(best_total_deaths["sim"].results["new_infections"])
    best_interv_cum = np.cumsum(best_shortest_intervention["sim"].results["new_infections"])

    # Plot baseline
    ax.plot(days, baseline_cum, c="black", lw=.8, label="Baseline (no intervention)")

    # Plot policy mean
    ax.plot(days, policy_mean_cum, c="blue", lw=.8, label="Policy mean")

    # Plot best runs
    ax.plot(days, best_inf_cum, c="red", lw=.8, label="Best: total infections")
    ax.plot(days, best_deaths_cum, c="orange", lw=.8, label="Best: total deaths")
    ax.plot(days, best_interv_cum, c="green", lw=.8, label="Best: shortest intervention duration")

    ax.set_xlabel("Day")
    ax.set_ylabel("Cumulative infections")
    ax.set_title(title)
    ax.legend()
    ax.grid(True)

    return fig, ax
