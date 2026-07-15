# policies/policy_classes.py

class Policy:
    """
    Defines a single surveillance policy:
      - testing probability
      - sequencing probability
      - delay PMF
      - detection threshold
      - quiet period
      - intervention parameters (optional)
    """

    def __init__(
        self,
        p_test=0.1,
        p_seq=0.1,
        delay_pmf=None,
        detection_threshold=5,
        quiet_period=0,
        beta_reduction=0.3,
        beta_min=0.005,
        beta_restore=1.0,
        intervention_duration=14,
        label=None,
    ):
        self.p_test = p_test
        self.p_seq = p_seq

        self.delay_pmf = delay_pmf if delay_pmf is not None else [1.0]
        self.detection_threshold = detection_threshold
        self.quiet_period = quiet_period

        # intervention parameters
        self.beta_reduction = beta_reduction
        self.beta_min = beta_min
        self.beta_restore = beta_restore
        self.intervention_duration = intervention_duration

        # optional label for logging
        self.label = label or f"policy_pt{p_test}_ps{p_seq}_thr{detection_threshold}"

    def to_dict(self):
        """
        Convert to dict for ObservationModel and InterventionController.
        """
        return {
            "p_test": self.p_test,
            "p_seq": self.p_seq,
            "delay_pmf": self.delay_pmf,
            "detection_threshold": self.detection_threshold,
            "quiet_period": self.quiet_period,
            "beta_reduction": self.beta_reduction,
            "beta_min": self.beta_min,
            "beta_restore": self.beta_restore,
            "intervention_duration": self.intervention_duration,
        }
