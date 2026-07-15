import covasim as cv

class InterventionTracker:
    """
    Tracks intervention episodes:
    - start day
    - end day
    - intensity (beta reduction factor)
    """

    def __init__(self):
        self.active = False
        self.current_start = None
        self.current_intensity = None
        self.log = []  # list of (start_day, end_day, intensity)

    def start(self, day, intensity): 
        # Start a new intervention episode
        self.active = True
        self.current_start = day
        self.current_intensity = intensity

    def end(self, day):
        # End the current intervention episode
        if self.active:
            self.log.append((self.current_start, day, self.current_intensity))
        self.active = False
        self.current_start = None
        self.current_intensity = None

    def is_active(self):
        return self.active


# Global tracker instance (runner.py will import this)
tracker = InterventionTracker()


def reduce_beta(sim, factor, day=None):
    """
    Reduce transmission rate by multiplying beta with a factor.
    Example: factor = 0.7 → 30% reduction.
    """

    # Apply intervention
    sim.pars['beta'] *= factor

    # Record intensity
    if day is not None:
        if not tracker.is_active():
            tracker.start(day, factor)
        else:
            # If already active, update intensity if needed
            tracker.current_intensity = factor


def increase_testing(sim, test_prob, day=None):
    """
    Increase daily testing probability.
    """

    sim.pars['test_prob'] = test_prob

    # Testing intervention is also recorded as an episode
    if day is not None:
        if not tracker.is_active():
            tracker.start(day, test_prob)
        else:
            tracker.current_intensity = test_prob


def maybe_end_intervention(day, last_trigger_day, quiet_period=5):
    """
    End intervention if no triggers have occurred for 'quiet_period' days.
    """

    if tracker.is_active():
        if day - last_trigger_day >= quiet_period:
            tracker.end(day)
