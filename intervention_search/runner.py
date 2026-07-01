# This script connects the following steps
# 1. Covasim step() execution
# 2. ObservationModel (diagnosis generation)
# 3. Trigger functions
# 4. Intervention functions

import numpy as np
import covasim as cv

from observation_model import ObservationModel
from triggers import relative_threshold, weekly_growth, slope_trigger, sustained_increase
from interventions import reduce_beta


def run_simulation(
    sim_pars, # Parameters passed to cv.Sim() in dict
    obs_pars, # Parameters for ObservationModel
    trigger_list, # List of (name, fn, params). 
                  # Ex. [
                  #     ("abs50", absolute_threshold, {"threshold": 50, "window": 7}),
                  #     ("growth1.3", weekly_growth, {"ratio": 1.3})
                  #   ]
                  # Checking absolute threshold of 50 and name it 'abs50'
                  # and weekly growth of ratio 1.3 and name it 'growth1.3'.
    intervention_fn, # Function that modifies sim when a trigger fires. Ex. increase_testing(sim, test_prob)
):

    # Returns a dict of 
    # {
    #     "diagnoses": array,
    #     "sequences": array,
    #     "trigger_events": list of (day, trigger_name),
    #     "sim": sim
    # }


    # 1) Initialize simulation
    sim = cv.Sim(sim_pars)
    sim.initialize()

    T = sim.npts

    # 2) Initialize observation model
    obs = ObservationModel(
        T=T,
        p_test=obs_pars["p_test"],
        p_seq=obs_pars["p_seq"],
        delay_pmf=obs_pars["delay_pmf"],
    )

    # 3) Prepare arrays
    diagnoses = np.zeros(T)
    sequences = np.zeros(T)
    trigger_events = []

    # 4) Step loop
    for t in range(T):
        sim.step()

        # Observation model
        new_inf_today = sim.results["new_infections"][t]
        diag_today, seq_today = obs.step(t, new_inf_today)

        diagnoses[t] = diag_today
        sequences[t] = seq_today

        # Trigger check
        for name, fn, params in trigger_list:
            if "sim" in params:
                if fn(diagnoses, t, **params):
                    trigger_events.append((t, name))
                    intervention_fn(sim)
            else:
                if fn(diagnoses, t, **params):
                    trigger_events.append((t, name))
                    intervention_fn(sim)   # apply intervention immediately

    # 5) Return results
    return {
        "diagnoses": diagnoses,
        "sequences": sequences,
        "trigger_events": trigger_events,
        "sim": sim,
    }
