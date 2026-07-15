import numpy as np

class ObservationModel:
    """
    Observation model:
      - diagnoses sampling (Binomial)
      - delay PMF (Multinomial)
      - sequencing sampling (Binomial)
      - GT trigger (true infections)
      - D trigger (diagnoses)
      - GT/D trigger logs (lists)
    """

    def __init__(self, policy, sim_npts=None):
        self.p_test = policy.get("p_test", 0.1)
        self.p_seq = policy.get("p_seq", 0.1)


        # delay PMF
        self.delay_pmf = np.array(policy.get("delay_pmf", [1.0]))
        self.max_delay = len(self.delay_pmf)

        self.sim_npts = sim_npts

        self.daily_diagnoses = None
        self.daily_sequences = None
        self.daily_sequenced_agents = None

        # daily flags
        self.gt_trigger_today = False
        self.diag_trigger_today = False

        # trigger logs (NEW)
        self.gt_triggers = []
        self.diag_triggers = []

    def initialize(self, sim):
        T = sim.npts
        self.sim_npts = T

        self.daily_diagnoses = np.zeros(T, dtype=int)
        self.daily_sequences = np.zeros(T, dtype=int)
        self.daily_sequenced_agents = [[] for _ in range(T)]

    def apply(self, sim, t):
        """
        Called once per day by runner.py.
        """

        # reset daily flags
        self.gt_trigger_today = False
        self.diag_trigger_today = False

        # 1) true infections
        true_inf = int(sim.results["new_infections"][t])

        # 2) diagnoses sampling
        diagnoses_raw = np.random.binomial(true_inf, self.p_test)
        if diagnoses_raw == 0:
            return

        # 3) distribute diagnoses across future days
        diag_by_delay = np.random.multinomial(diagnoses_raw, self.delay_pmf)

        for k, count in enumerate(diag_by_delay):
            if count == 0:
                continue

            diag_day = t + k
            if diag_day >= self.sim_npts:
                break

            # diagnoses
            self.daily_diagnoses[diag_day] += count


            # sequences
            seq_count = np.random.binomial(count, self.p_seq)
            self.daily_sequences[diag_day] += seq_count

            if seq_count > 0:
                infectious_agents = np.where(sim.people.infectious)[0]
                if len(infectious_agents) > 0:
                    chosen = np.random.choice(
                        infectious_agents,
                        size=min(seq_count, len(infectious_agents)),
                        replace=False,
                    )
                    self.daily_sequenced_agents[diag_day].extend(chosen.tolist())
