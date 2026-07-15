# intervention_search/simulation/intervention.py

class InterventionController:
    """
    Handles:
      - beta reduction when intervention starts
      - intervention duration
      - (NO RESTORE) beta stays reduced; Covasim dynamics handle future changes
    """

    def __init__(self, policy):
        self.beta_reduction = policy.get("beta_reduction", 0.3)
        self.beta_min = policy.get("beta_min", 0.005)
        self.intervention_duration = policy.get("intervention_duration", 14)

        # internal state
        self.active = False
        self.intervention_day = None

        self.intervention_days = []
        self.gt_triggers = []
        self.diag_triggers = []

    def start(self, sim, t):
        """
        Called by Runner when an intervention is triggered.
        """
        if self.active:
            return

        self.intervention_day = t
        self.intervention_days.append(t)   # NEW: record intervention start

        self._activate_intervention(sim)
        self.active = True

    def update(self, sim, t):
        """
        Called every day by Runner.
        Ends intervention after duration, but DOES NOT restore beta.
        """
        if not self.active:
            return

        if t >= self.intervention_day + self.intervention_duration:
            # simply end intervention; do NOT restore beta
            self.active = False

    def _activate_intervention(self, sim):
        """
        Reduce beta when intervention starts.
        """
        new_beta = max(sim.pars["beta"] * (1 - self.beta_reduction), self.beta_min)
        sim.pars["beta"] = new_beta
