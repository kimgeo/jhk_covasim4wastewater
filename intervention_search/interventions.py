import covasim as cv

# These functions are designed to be called inside a step() loop when a trigger fires. 

def reduce_beta(sim, factor): # direct parameter modification for beta
    # Multiply the transmission rate beta by a given factor.
    # Ex. factor = 0.7 means 30% reduction.
    sim.pars['beta'] *= factor


def increase_testing(sim, test_prob): # direct parameter modification for test_prop
    # Increase daily testing probability.
    sim.pars['test_prob'] = test_prob