import numpy as np

class ObservationModel:
    # Observation model for step-based Covasim simulation.

    # This class maintains internal queues for delayed diagnoses and sequencing.
    # Each day, you call step(t, new_inf_today) to process today's infections
    # and retrieve today's observed diagnoses and sequences.

    def __init__(self, T, p_test, p_seq, delay_pmf):
        self.T = T # Total simulation length (sim.npts).
        self.p_test = p_test
        self.p_seq = p_seq
        self.delay_pmf = np.asarray(delay_pmf, dtype=float)

        # Queues for scheduled future diagnoses and sequences
        self.diagnosis_queue = np.zeros(T)
        self.sequence_queue = np.zeros(T)

    def step(self, t, new_inf_today):
        # Process today's infections and return today's observed diagnoses & sequences.  

        # 1) Determine how many will be diagnosed from today's infections
        if new_inf_today > 0:
            n_diag = np.random.binomial(int(new_inf_today), self.p_test)

            # 2) Assign diagnoses delay according to delay_pmf
            diag_by_delay = np.random.multinomial(n_diag, self.delay_pmf)

            for k, count in enumerate(diag_by_delay):
                day = t + k
                if day < self.T:
                    self.diagnosis_queue[day] += count

        # 3) Diagnoses that occur today
        diag_today = self.diagnosis_queue[t] # Observed diagnoses on day t

        # 4) Sequencing among today's diagnoses
        if diag_today > 0:
            seq_today = np.random.binomial(int(diag_today), self.p_seq) # Observed sequences on day t
        else:
            seq_today = 0

        return diag_today, seq_today
