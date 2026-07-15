"""
run_grid.py

Run a grid search over trigger x intervention combinations.
Compute metrics, select best runs, and visualize results.
"""

import numpy as np
import covasim as cv

from runner import run_simulation
from jhk_covasim4wastewater.intervention_search_prev.evaluation import evaluate_results
from plot_result import plot_top3_and_baselines

# Trigger functions
from triggers import weekly_growth, relative_threshold

# Intervention functions
from interventions import reduce_beta, increase_testing

# 1) Define baseline simulation

def run_baseline(sim_pars):
    baseline = cv.Sim(sim_pars)
    baseline.run()
    return baseline


# 2) Define trigger grid

def build_trigger_grid(sim):
    """
    Build a list of trigger configurations.
    Each element is (name, fn, params)
    """

    trigger_grid = []

    # Example weekly growth triggers
    trigger_grid.append(("growth2", weekly_growth,
                         {"ratio": 2, "window": 7}))

    trigger_grid.append(("growth5", weekly_growth,
                         {"ratio": 5, "window": 7}))

    # Example relative threshold triggers
    trigger_grid.append(("rel0.1", relative_threshold,
                         {"proportion": 0.1, "window": 7, "sim": sim}))

    trigger_grid.append(("rel0.5", relative_threshold,
                         {"proportion": 0.5, "window": 7, "sim": sim}))

    return trigger_grid


# 3) Define intervention grid

def build_intervention_grid():
    """
    Build a list of intervention functions.
    Each element is a callable: lambda sim: ...
    """

    interventions = []

    interventions.append(lambda sim: reduce_beta(sim, 0.7))
    interventions.append(lambda sim: reduce_beta(sim, 0.5))
    interventions.append(lambda sim: increase_testing(sim, 0.4))

    return interventions


# 4) Grid search

def run_grid_search(sim_pars, obs_pars):
    """
    Run all trigger x intervention combinations.
    Return list of results dicts.
    """

    # Build baseline sim (for plotting)
    baseline_sim = run_baseline(sim_pars)

    # Build trigger grid (needs sim object for relative_threshold)
    dummy_sim = cv.Sim(sim_pars)
    trigger_grid = build_trigger_grid(dummy_sim)

    # Build intervention grid
    intervention_grid = build_intervention_grid()

    all_results = []

    # Loop over all combinations
    for (name, fn, params) in trigger_grid:
        for intervention_fn in intervention_grid:

            results = run_simulation(
                sim_pars=sim_pars,
                obs_pars=obs_pars,
                trigger_list=[(name, fn, params)],
                intervention_fn=intervention_fn,
            )

            metrics = evaluate_results(results)

            all_results.append({
                "name": name,
                "intervention_fn": intervention_fn,
                "results": results,
                "metrics": metrics,
            })

            print(f"Finished: trigger={name}, intervention={intervention_fn}")


    return baseline_sim, all_results


# 5) Select best runs

def select_best_runs(all_results):
    """
    Select best runs by:
    - total infections
    - total deaths
    - shortest intervention duration
    """

    # total infections
    best_total_inf = min(all_results,
                         key=lambda r: r["metrics"]["total_infections"])

    # total deaths
    best_total_deaths = min(all_results,
                            key=lambda r: r["metrics"]["total_sequences"])

    # shortest intervention duration
    best_shortest_intervention = min(all_results,
                                     key=lambda r: len(r["results"]["trigger_events"]))

    # ---- summary printing ----
    def summarize(run):
        return {
            "trigger_name": run["name"],
            "total_infections": run["metrics"]["total_infections"],
            "total_sequences": run["metrics"]["total_sequences"],
            "num_triggers": len(run["results"]["trigger_events"]),
            "trigger_days": run["results"]["trigger_events"],
            "intervention_fn": run["intervention_fn"],
        }

    print("\n=== Best by total infections ===")
    print(summarize(best_total_inf))

    print("\n=== Best by total deaths ===")
    print(summarize(best_total_deaths))

    print("\n=== Best by shortest intervention duration ===")
    print(summarize(best_shortest_intervention))

    return best_total_inf, best_total_deaths, best_shortest_intervention



# 6) Main entry

def main():

    # Example sim parameters
    sim_pars = {
        "pop_size": 100000,
        "pop_infected": 50,
        "beta": 0.015,
        "n_days": 120,
    }

    # Example observation model parameters
    obs_pars = {
        "p_test": 0.2,
        "p_seq": 0.1,
        "delay_pmf": [0.5, 0.3, 0.2],
    }

    # Run grid search
    baseline_sim, all_results = run_grid_search(sim_pars, obs_pars)

    # Extract only results dicts for plotting mean
    policy_results_list = [r["results"] for r in all_results]

    # Select best runs
    best_total_inf, best_total_deaths, best_shortest_intervention = \
        select_best_runs(all_results)

    # Plot
    fig, ax = plot_top3_and_baselines(
        baseline_sim=baseline_sim,
        policy_results_list=policy_results_list,
        best_total_inf=best_total_inf["results"],
        best_total_deaths=best_total_deaths["results"],
        best_shortest_intervention=best_shortest_intervention["results"],
    )

    fig.savefig("grid_search_results.png")
    print("Saved plot: grid_search_results.png")


if __name__ == "__main__":
    main()
