import covasim as cv
import numpy as np

from observation_model import ObservationModel
from interventions import reduce_beta, maybe_end_intervention, tracker
from triggers import check_triggers  
from strategy import Strategy      
from variant_utils import compute_variant_prevalence


def run_strategy(strategy: Strategy, sim_pars, n_runs=1):
    """
    Run a surveillance + trigger + variant strategy.
    Returns a dict of results for evaluation + plotting.
    """

    all_results = []

    for run in range(n_runs):
        # 1) Build simulation
        sim = cv.Sim(sim_pars)
        # 2) Variant scenario
        vp = strategy.variant_pars

        if vp is not None:
            variant = cv.variant(
                label='variantX',
                rel_beta=vp['advantage'],
                start_day=vp['intro_day']
            )
            sim.pars['variants'] = [variant]

        # 3) Run simulation day-by-day
        sim.initialize()

        # Observation model
        obs = ObservationModel(
            p_test=strategy.obs_pars['p_test'],
            p_seq=strategy.obs_pars['p_seq'],
            delay_pmf=strategy.obs_pars.get('delay_pmf', None),
            variant_ref=strategy.variant_pars.get('reference', "REF")
        )

        # Trigger tracking
        trigger_days = []
        last_trigger_day = -999

        # Intervention tracking (tracker imported)
        tracker.log = []
        tracker.active = False

        # Storage
        daily_diagnoses = []
        daily_sequences = []
        daily_variant_prev = []

        # Simulation loop
        for t in range(sim.npts):

            sim.step()  # advance Covasim one day

            # 1) Apply observation model for this day
            true_inf = sim.results['new_infections'][t]
            diag = np.random.binomial(true_inf, strategy.obs_pars['p_test'])
            seq = np.random.binomial(diag, strategy.obs_pars['p_seq'])

            daily_diagnoses.append(diag)
            daily_sequences.append(seq)

            # 2) Variant prevalence
            if 'variant' in sim.results:
                prev = compute_variant_prevalence(sim, t)
            else:
                prev = 0.0
            daily_variant_prev.append(prev)

            # 3) Check triggers
            trigger_fired = check_triggers(
                t=t,
                diagnoses=daily_diagnoses,
                sequences=daily_sequences,
                variant_prev=daily_variant_prev,
                trigger_list=strategy.trigger_list
            )

            if trigger_fired:
                trigger_days.append(t)
                last_trigger_day = t

                # Apply intervention
                reduce_beta(sim, factor=strategy.obs_pars['intervention_factor'], day=t)

            # 4) Check intervention end condition
            maybe_end_intervention(t, last_trigger_day, quiet_period=strategy.obs_pars['quiet_period'])

        # 4) After sim: full observation model pass
        obs_results = obs.apply(sim)

        # 5) Collect results
        result = {
            "name": strategy.name,
            "sim": sim,
            "diagnoses": daily_diagnoses,
            "sequences": daily_sequences,
            "variant_prevalence": daily_variant_prev,
            "variant_detected_day": obs_results["variant_detected_day"],
            "variant_infections": obs_results["variant_infections"],
            "trigger_days": trigger_days,
            "intervention_log": tracker.log,
        }

        all_results.append(result)

    return all_results
