import numpy as np

class ObservationModel:
    """
    Observation model:
    - diagnoses (testing)
    - sequencing (genomic surveillance)
    - variant detection
    """

    def __init__(self, p_test=0.1, p_seq=0.05, delay_pmf=None, variant_ref=None):
        """
        p_test: probability of diagnosis
        p_seq: probability of sequencing among diagnosed
        delay_pmf: diagnosis delay distribution
        variant_ref: reference haplotype (string or array)
        """
        self.p_test = p_test
        self.p_seq = p_seq
        self.delay_pmf = delay_pmf
        self.variant_ref = variant_ref

        # Variant detection tracking
        self.variant_detected_day = None
        self.variant_infections = None  # filled after sim run

    def reconstruct_haplotype(self, genome):
        """
        Placeholder haplotype reconstruction.
        In real use: apply actual reconstruction logic.
        """
        return genome  # assume genome is already haplotype-like

    def is_variant(self, haplotype):
        """
        Compare reconstructed haplotype with reference.
        If different → variant detected.
        """
        return haplotype != self.variant_ref

    def apply(self, sim):
        """
        Apply observation model to simulation results.
        Called AFTER sim.run().
        """

        n_days = sim.npts

        # 1) Diagnoses
        true_infections = sim.results['new_infections'].values
        diagnoses = np.random.binomial(true_infections, self.p_test)

        # 2) Sequencing among diagnosed
        sequences = np.random.binomial(diagnoses, self.p_seq)

        # 3) Variant infections per day (true infections)
        if 'variant' in sim.results:
            # Covasim variant tracking
            self.variant_infections = sim.results['variant']['new_infections_by_variant'][:, 1]
        else:
            # If variant tracking not enabled
            self.variant_infections = np.zeros(n_days)

        # 4) Variant detection from sequencing
        variant_detected = False

        for day in range(n_days):
            if sequences[day] > 0:
                # For each sequenced sample, reconstruct haplotype
                for _ in range(sequences[day]):
                    # Fake genome: assume sim stores variant genomes
                    genome = sim.people.genomes[sim.t] if hasattr(sim.people, 'genomes') else "REF"

                    hap = self.reconstruct_haplotype(genome)

                    if self.is_variant(hap):
                        variant_detected = True
                        if self.variant_detected_day is None:
                            self.variant_detected_day = day
                        break

            if variant_detected:
                break

        return {
            "diagnoses": diagnoses,
            "sequences": sequences,
            "variant_detected_day": self.variant_detected_day,
            "variant_infections": self.variant_infections,
        }
