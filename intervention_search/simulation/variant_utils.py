# intervention_search/simulation/variant_utils.py

import numpy as np

class VariantMetrics:
    """
    Tracks variant prevalence over time.
    Two views:
      - population-wide variant prevalence (label-based, from sim.people.variant)
      - sequencing-based variant counts (restricted to sequenced agents)
    """

    def __init__(self):
        # variant_name -> list of prevalence over time (population-based)
        self.prevalence = {}

        # variant_name -> list of counts over time (population-based)
        self.counts = {}

        # variant_name -> list of counts over time among sequenced agents
        self.seq_counts = {}

        # variant_name -> list of prevalence over time among sequenced agents
        self.seq_prevalence = {}

        self.variant_names = None
        self.N = None

    def initialize(self, sim):
        """
        Called once at the start of the run.
        Extract variant names from sim.
        """
        variants = sim.people.variant
        self.variant_names = sorted(list(set(variants)))
        self.N = sim.pars["pop_size"]

        for v in self.variant_names:
            self.prevalence[v] = []
            self.counts[v] = []
            self.seq_counts[v] = []
            self.seq_prevalence[v] = []

    def update(self, sim, t, sequenced_agents=None):
        """
        Called once per day by runner.py.
        Compute prevalence for each variant.
        Optionally, also compute variant distribution among sequenced agents.
        """

        # today's variant labels
        variants_today = sim.people.variant

        # population-wide counts and prevalence
        for v in self.variant_names:
            count_v = np.sum(variants_today == v)
            self.counts[v].append(int(count_v))

            prev_v = count_v / self.N
            self.prevalence[v].append(float(prev_v))

        # sequencing-based view (if sequenced_agents provided)
        if sequenced_agents is not None and len(sequenced_agents) > 0:
            # total sequenced today
            S = len(sequenced_agents)

            # variant labels for sequenced agents
            seq_variants = variants_today[sequenced_agents]

            for v in self.variant_names:
                seq_count_v = np.sum(seq_variants == v)
                self.seq_counts[v].append(int(seq_count_v))
                self.seq_prevalence[v].append(seq_count_v / S)
        else:
            # if no sequenced agents today, append zeros
            for v in self.variant_names:
                self.seq_counts[v].append(0)
                self.seq_prevalence[v].append(0.0)
