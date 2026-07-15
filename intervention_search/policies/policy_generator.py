# policies/policy_generator.py

import itertools
from .policy_classes import Policy

class PolicyGenerator:
    """
    Generates multiple policy combinations for grid search or random search.
    """

    def __init__(self):
        pass

    def grid_search(
        self,
        p_test_values,
        p_seq_values,
        threshold_values,
        delay_pmfs=None,
        quiet_periods=None,
    ):
        """
        Generate all combinations of policies.
        """
        delay_pmfs = delay_pmfs or [[1.0]]
        quiet_periods = quiet_periods or [0]

        policies = []

        for pt, ps, thr, pmf, qp in itertools.product(
            p_test_values,
            p_seq_values,
            threshold_values,
            delay_pmfs,
            quiet_periods,
        ):
            policies.append(
                Policy(
                    p_test=pt,
                    p_seq=ps,
                    delay_pmf=pmf,
                    detection_threshold=thr,
                    quiet_period=qp,
                    label=f"pt{pt}_ps{ps}_thr{thr}_qp{qp}"
                )
            )

        return policies

    def random_search(self, n=20):
        """
        Randomly sample policies.
        """
        import numpy as np

        policies = []
        for i in range(n):
            pt = np.random.uniform(0.01, 0.5)
            ps = np.random.uniform(0.01, 0.5)
            thr = np.random.randint(1, 20)
            qp = np.random.randint(0, 7)

            policies.append(
                Policy(
                    p_test=pt,
                    p_seq=ps,
                    detection_threshold=thr,
                    quiet_period=qp,
                    label=f"rand_{i}"
                )
            )
        return policies
