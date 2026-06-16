import covasim as cv
import sciris as sc
import numpy as np
import matplotlib.pyplot as plt

def make_korea_sim(
    start_day='2020-02-15',
    end_day='2020-12-31',
    pop_size=200_000,      # Number of simulated agents (scaled population)
    pop_scale=250,         # Real population = pop_size * pop_scale
    pop_infected=50,       # Initial number of infected individuals
    rand_seed=42,
    datafile=None,         # Path to Korea case/death CSV (optional, for calibration)
):
    """
    Covasim simulation template configured for Korea.
    Real Korean data and policy timelines can be added later.
    """

    pars = dict(
        # Population and timeline settings
        pop_type     = 'hybrid',      # Includes household, school, workplace, community layers
        pop_size     = pop_size,
        pop_scale    = pop_scale,
        pop_infected = pop_infected,
        rand_seed    = rand_seed,

        # Dates
        start_day    = start_day,
        end_day      = end_day,

        # Epidemiological parameters
        beta         = 0.015,         # Baseline transmission probability (placeholder)
        n_days       = None,          # Automatically computed from end_day

        # Output and performance
        verbose      = 0,
    )

    # Time‑dependent interventions (example placeholders)
    interventions = []

    # Example: Social distancing strengthened on March 1, 2020 → reduce transmission
    interventions.append(
        cv.change_beta(
            days=['2020-03-01'],
            changes=[0.6],
        )
    )

    # Example: Increased mobility and clusters around Aug 15, 2020 → increase transmission
    interventions.append(
        cv.change_beta(
            days=['2020-08-15'],
            changes=[1.2],
        )
    )

    # Attach real data if provided
    if datafile is not None:
        pars['datafile'] = datafile

    sim = cv.Sim(pars=pars, interventions=interventions, label='Korea baseline')

    return sim


if __name__ == '__main__':
    sim = make_korea_sim(
        start_day='2020-02-15',
        end_day='2020-12-31',
        pop_size=200_000,
        pop_scale=250,
        pop_infected=50,
        rand_seed=42,
        datafile=None,  # Example: "data/korea_timeseries.csv"
    )

    # Run simulation
    sim.run()

    # Plot results
    fig = sim.plot()

    # Save plot
    fig.savefig("korea_simulation_wo_cal.png", dpi=300, bbox_inches='tight')

    print("Plot saved as korea_simulation_wo_cal.png")
