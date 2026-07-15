# intervention_search/simulation/runner.py

import os
import json
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List

import covasim as cv
import numpy as np

from .observation_model import ObservationModel
from .intervention import InterventionController
from .variant_utils import VariantMetrics
from .trigger_functions import check_triggers


@dataclass
class RunConfig:
    sim_pars: Dict[str, Any]
    evo_pars: Optional[Dict[str, Any]] = None
    policy: Dict[str, Any] = None
    run_id: str = "run_000"
    out_dir: str = "results/raw"


@dataclass
class RunResults:
    run_id: str
    days: List[int] = field(default_factory=list)

    true_infections: List[int] = field(default_factory=list)
    diagnoses: List[int] = field(default_factory=list)
    sequences: List[int] = field(default_factory=list)

    variant_prevalence: Dict[str, List[float]] = field(default_factory=dict)

    detection_day: Optional[int] = None
    intervention_day: Optional[int] = None

    beta: List[float] = field(default_factory=list)
    cum_infections: List[int] = field(default_factory=list)

    # trigger logs
    gt_triggers: List[int] = field(default_factory=list)
    diag_triggers: List[int] = field(default_factory=list)
    intervention_days: List[int] = field(default_factory=list)

    # GT→D resolution metrics
    resolution_times: List[int] = field(default_factory=list)
    resolution_gt_counts: List[int] = field(default_factory=list)
    avg_resolution_time: Optional[float] = None
    avg_cumulative_gt: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "run_id": self.run_id,
            "days": self.days,
            "true_infections": self.true_infections,
            "diagnoses": self.diagnoses,
            "sequences": self.sequences,
            "variant_prevalence": self.variant_prevalence,
            "detection_day": self.detection_day,
            "intervention_day": self.intervention_day,
            "beta": self.beta,
            "cum_infections": self.cum_infections,
            "gt_triggers": self.gt_triggers,
            "diag_triggers": self.diag_triggers,
            "intervention_days": self.intervention_days,
            "resolution_times": self.resolution_times,
            "resolution_gt_counts": self.resolution_gt_counts,
            "avg_resolution_time": self.avg_resolution_time,
            "avg_cumulative_gt": self.avg_cumulative_gt,
        }


class Runner:
    """
    Full GT/D trigger + intervention + resolution logic,
    with optional evo_pars support.
    Advanced trigger is applied to both:
      - GT: true infections (sim.results["new_infections"])
      - D: diagnoses (ObservationModel.daily_diagnoses)
    """

    def __init__(self, config: RunConfig):
        self.config = config
        self.sim: cv.Sim = self._build_sim(config.sim_pars, config.evo_pars)

        self.obs = ObservationModel(policy=config.policy)
        self.obs.initialize(self.sim)

        self.intervention = InterventionController(policy=config.policy)

        # self.variant_metrics = VariantMetrics()
        # self.variant_metrics.initialize(self.sim)

        self.results = RunResults(run_id=config.run_id)

        # trigger logs
        self.gt_triggers: List[int] = []
        self.diag_triggers: List[int] = []
        self.intervention_days: List[int] = []

        # policy options
        self.max_interventions = config.policy.get("max_interventions", 1)
        self.quiet_period = config.policy.get("quiet_period", 0)

        # first detection (first D advanced trigger)
        self.first_diag_trigger_day: Optional[int] = None

    def _build_sim(self, sim_pars: Dict[str, Any], evo_pars: Optional[Dict[str, Any]]) -> cv.Sim:
        pars = sim_pars.copy()
        if evo_pars is not None:
            pars["evo_pars"] = evo_pars
        sim = cv.Sim(pars=pars)
        sim.initialize()
        return sim

    def run(self) -> RunResults:
        n_days = int(self.sim.pars["n_days"])
        trigger_list = self.config.policy.get("trigger_list", [])

        for t in range(n_days):
            self.sim.step()

            # 1) observation model (diagnoses/sequences generation)
            self.obs.apply(self.sim, t)

            diagnoses = self.obs.daily_diagnoses
            sequences = self.obs.daily_sequences

            # 2) variant prevalence series (per variant; choose one or default)
            variant_prev = None
            # variant_prev = self.variant_metrics.prevalence.get(
            #     self.config.policy.get("variant_name", "variantX"),
            #     []
            # )

            # 3) GT advanced trigger: apply trigger to true infections
            true_inf_series = self.sim.results["new_infections"]
            if check_triggers(
                t,
                true_inf_series,
                sequences,
                variant_prev,
                trigger_list,
                sim=self.sim,
            ):
                self.gt_triggers.append(t)

            # 4) D advanced trigger: apply trigger to diagnoses and start intervention
            if check_triggers(
                t,
                diagnoses,
                sequences,
                variant_prev,
                trigger_list,
                sim=self.sim,
            ):
                self.diag_triggers.append(t)

                if self.first_diag_trigger_day is None:
                    self.first_diag_trigger_day = t

                if len(self.intervention_days) < self.max_interventions:
                    if t >= self.first_diag_trigger_day + self.quiet_period:
                        self._start_intervention(t)

            # 5) update intervention (duration only, no restore)
            self.intervention.update(self.sim, t)

            # 6) variant metrics update
            # self.variant_metrics.update(
            #     self.sim,
            #     t,
            #     sequenced_agents=self.obs.daily_sequenced_agents[t],
            # )

            # 7) record daily metrics
            self._record_day(t)

        # compute GT→D resolution metrics
        self._compute_resolution_metrics()

        # expose triggers on intervention for plotting
        self.intervention.gt_triggers = self.gt_triggers
        self.intervention.diag_triggers = self.diag_triggers
        self.intervention.intervention_days = self.intervention_days

        # save raw
        self._save_raw()
        return self.results

    def _start_intervention(self, t: int):
        if self.intervention.active:
            return
        self.intervention.start(self.sim, t)
        self.intervention_days.append(t)

    def _record_day(self, t: int) -> None:
        self.results.days.append(t)

        new_inf = int(self.sim.results["new_infections"][t])
        self.results.true_infections.append(new_inf)

        cum_inf = int(self.sim.results["cum_infections"][t])
        self.results.cum_infections.append(cum_inf)

        self.results.diagnoses.append(self.obs.daily_diagnoses[t])
        self.results.sequences.append(self.obs.daily_sequences[t])

        # for vname, series in self.variant_metrics.prevalence.items():
        #     if vname not in self.results.variant_prevalence:
        #         self.results.variant_prevalence[vname] = []
        #     self.results.variant_prevalence[vname].append(series[t])

        self.results.beta.append(float(self.sim.pars["beta"]))

        # detection_day = first D advanced trigger
        if self.results.detection_day is None and self.first_diag_trigger_day is not None:
            self.results.detection_day = self.first_diag_trigger_day

        # intervention_day = first intervention
        if self.results.intervention_day is None and len(self.intervention_days) > 0:
            self.results.intervention_day = self.intervention_days[0]

    def _compute_resolution_metrics(self) -> None:
        gt = self.gt_triggers
        diag = self.diag_triggers

        resolution_times: List[int] = []
        resolution_gt_counts: List[int] = []

        i = 0
        j = 0

        while i < len(gt) and j < len(diag):
            g_first = gt[i]
            d = diag[j]

            group_count = 0
            while i < len(gt) and gt[i] <= d:
                group_count += 1
                i += 1

            resolution_times.append(d - g_first)
            resolution_gt_counts.append(group_count)

            j += 1

        self.results.gt_triggers = gt
        self.results.diag_triggers = diag
        self.results.intervention_days = self.intervention_days
        self.results.resolution_times = resolution_times
        self.results.resolution_gt_counts = resolution_gt_counts

        if resolution_times:
            self.results.avg_resolution_time = sum(resolution_times) / len(resolution_times)
        if resolution_gt_counts:
            self.results.avg_cumulative_gt = sum(resolution_gt_counts) / len(resolution_gt_counts)

    def _save_raw(self) -> None:
        out_dir = self.config.out_dir
        os.makedirs(out_dir, exist_ok=True)

        out_path = os.path.join(out_dir, f"{self.config.run_id}.json")

        def convert(o):
            if isinstance(o, np.integer):
                return int(o)
            if isinstance(o, np.floating):
                return float(o)
            if isinstance(o, np.ndarray):
                return o.tolist()
            return o

        with open(out_path, "w") as fh:
            json.dump(self.results.to_dict(), fh, indent=2, default=convert)


def run_single(sim_pars: Dict[str, Any], policy: Dict[str, Any], run_id: str) -> RunResults:
    cfg = RunConfig(sim_pars=sim_pars, policy=policy, run_id=run_id)
    runner = Runner(cfg)
    return runner.run()


def run_evo_simulation(sim_pars: Dict[str, Any],
                       evo_pars: Dict[str, Any],
                       policy: Dict[str, Any],
                       rand_seed: int = 0):
    pars = sim_pars.copy()
    pars["rand_seed"] = rand_seed

    cfg = RunConfig(sim_pars=pars, evo_pars=evo_pars, policy=policy, run_id=f"seed_{rand_seed}")
    runner = Runner(cfg)
    results = runner.run()

    return runner.sim, runner.obs, runner.intervention, None, runner.results