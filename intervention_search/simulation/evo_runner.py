# intervention_search/simulation/evo_runner.py

import covasim as cv
from typing import Dict, Any, Tuple

from .observation_model import ObservationModel
from .intervention import InterventionController
from .variant_utils import VariantMetrics


def run_evo_simulation(sim_pars: Dict[str, Any],
                       evo_pars: Dict[str, Any],
                       policy: Dict[str, Any],
                       rand_seed: int = 0):

    # 1) build sim with evo_pars
    pars = sim_pars.copy()
    pars["evo_pars"] = evo_pars
    pars["rand_seed"] = rand_seed

    sim = cv.Sim(pars=pars)
    sim.initialize()

    # 2) attach observation + intervention
    obs = ObservationModel(policy=policy)
    obs.initialize(sim)

    interv = InterventionController(policy=policy)

    # 3) variant metrics
    variant_metrics = None
    # variant_metrics = VariantMetrics()
    # variant_metrics.initialize(sim)

    # 4) trigger logs
    gt_triggers = []
    diag_triggers = []
    intervention_days = []

    # policy options
    max_interventions = policy.get("max_interventions", 1)
    apply_every_diag_trigger = policy.get("apply_every_diag_trigger", False)
    quiet_period = policy.get("quiet_period", 0)

    # first detection (first D trigger)
    first_diag_trigger_day = None

    # 5) run day-by-day
    n_days = int(sim_pars["n_days"])

    for t in range(n_days):
        sim.step()

        # observation model
        obs.apply(sim, t)

        # GT trigger
        if obs.gt_trigger_today:
            gt_triggers.append(t)

        # D trigger
        if obs.diag_trigger_today:
            diag_triggers.append(t)

            if first_diag_trigger_day is None:
                first_diag_trigger_day = t

            # intervention option B: every diag trigger
            if apply_every_diag_trigger:
                if len(intervention_days) < max_interventions:
                    if t >= (t + quiet_period):  # quiet_period relative to this D trigger
                        interv.start(sim, t)
                        intervention_days.append(t)

            # option A: quiet_period relative to first detection
            else:
                if first_diag_trigger_day is not None and t >= first_diag_trigger_day + quiet_period:
                    if len(intervention_days) < max_interventions:
                        interv.start(sim, t)
                        intervention_days.append(t)

        # update intervention (duration only, no restore)
        interv.update(sim, t)

        # variant metrics
        # variant_metrics.update(sim, t, sequenced_agents=obs.daily_sequenced_agents[t])

    # return everything
    return sim, obs, interv, variant_metrics
