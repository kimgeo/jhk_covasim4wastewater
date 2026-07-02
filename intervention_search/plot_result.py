import matplotlib.pyplot as plt
import numpy as np

def plot_result(result, save_path=None):
    """
    result: dict from runner.py
        - diagnoses
        - sequences
        - variant_prevalence
        - variant_detected_day
        - variant_infections
        - trigger_days
        - intervention_log
        - sim
    """

    diagnoses = np.array(result["diagnoses"])
    sequences = np.array(result["sequences"])
    variant_prev = np.array(result["variant_prevalence"])
    variant_inf = np.array(result["variant_infections"])
    trigger_days = result["trigger_days"]
    intervention_log = result["intervention_log"]
    variant_detected_day = result["variant_detected_day"]
    sim = result["sim"]

    days = np.arange(len(diagnoses))

    fig, axes = plt.subplots(5, 1, figsize=(12, 22))
    fig.suptitle(f"Policy Evaluation — {result['name']}", fontsize=18)

    # 1) Trigger timeline
    ax = axes[0]
    ax.plot(days, diagnoses, label="Diagnoses", color="blue")
    ax.scatter(trigger_days, diagnoses[trigger_days], color="red", label="Triggers", s=60)
    ax.set_title("Trigger Timeline")
    ax.set_xlabel("Day")
    ax.set_ylabel("Diagnoses")
    ax.legend()

    # 2) Intervention timeline
    ax = axes[1]
    ax.plot(days, diagnoses, label="Diagnoses", color="blue")

    for (start, end, intensity) in intervention_log:
        ax.axvspan(start, end, color="orange", alpha=0.3)
        ax.text(start, diagnoses[start], f"{intensity:.2f}", color="black")

    ax.set_title("Intervention Timeline")
    ax.set_xlabel("Day")
    ax.set_ylabel("Diagnoses")
    ax.legend()

    # 3) Variant detection timeline
    ax = axes[2]
    ax.plot(days, variant_prev, label="Variant Prevalence", color="purple")
    ax.plot(days, variant_inf, label="Variant Infections", color="green")

    if variant_detected_day is not None:
        ax.axvline(variant_detected_day, color="red", linestyle="--", label="Variant Detected")

    ax.set_title("Variant Detection Timeline")
    ax.set_xlabel("Day")
    ax.set_ylabel("Prevalence / Infections")
    ax.legend()

    # 4) Practicality vs Effectiveness scatter

    # Effectiveness: cumulative infections
    cum_inf = sim.results["new_infections"].values.sum()

    # Practicality: number of intervention episodes
    practicality = len(intervention_log)

    ax = axes[3]
    ax.scatter([practicality], [cum_inf], color="black", s=120)
    ax.set_title("Practicality vs Effectiveness")
    ax.set_xlabel("Practicality (Intervention Count)")
    ax.set_ylabel("Effectiveness (Cumulative Infections)")
    ax.grid(True)


    # 5) Surveillance burden vs infections
    burden = diagnoses.sum() + sequences.sum()

    ax = axes[4]
    ax.scatter([burden], [cum_inf], color="brown", s=120)
    ax.set_title("Surveillance Burden vs Infections")
    ax.set_xlabel("Surveillance Burden (Diagnoses + Sequences)")
    ax.set_ylabel("Cumulative Infections")
    ax.grid(True)

    # Save or show
    if save_path:
        plt.savefig(save_path, dpi=200)
    plt.show()
