import covasim as cv
import numpy as np
import matplotlib.pyplot as plt

def make_boarding_school_sim(
    start_day='2020-01-01',
    end_day='2021-12-31',
    n_students=200,
    dorm_size=20,
    class_size=25,
    n_teachers=20,
    pop_infected=1,
    rand_seed=42,
):
    """
    Covasim simulation for a 200-student boarding school (2020–2021).
    Demonstrates multiple custom contact layers:
    dorm, class, teacher, teach, common, random.
    """

    total_pop = n_students + n_teachers

    pars = dict(
        pop_size=total_pop,
        pop_scale=1,
        pop_type='random',      # we override with custom layers
        start_day=start_day,
        end_day=end_day,
        pop_infected=pop_infected,
        rand_seed=rand_seed,
        beta=0.012,
        verbose=0,
    )

    # Create simulation
    sim = cv.Sim(pars=pars)

    # IMPORTANT: initialize BEFORE modifying contacts
    sim.initialize()

    # ---------------------------
    # 1. Build custom layers
    # ---------------------------

    # Student indices: 0 ~ n_students-1
    # Teacher indices: n_students ~ total_pop-1

    # --- Dorm layer ---
    dorm_layer = []
    n_dorms = n_students // dorm_size
    for d in range(n_dorms):
        members = list(range(d*dorm_size, (d+1)*dorm_size))
        for i in members:
            for j in members:
                if i < j:
                    dorm_layer.append([i, j])

    # --- Class layer ---
    class_layer = []
    n_classes = n_students // class_size
    for c in range(n_classes):
        members = list(range(c*class_size, (c+1)*class_size))
        for i in members:
            for j in members:
                if i < j:
                    class_layer.append([i, j])

    # --- Teacher–teacher layer ---
    teacher_layer = []
    teacher_ids = list(range(n_students, total_pop))
    for i in teacher_ids:
        for j in teacher_ids:
            if i < j:
                teacher_layer.append([i, j])

    # --- Teacher–student (teaching) layer ---
    teach_layer = []
    for c in range(n_classes):
        class_members = list(range(c*class_size, (c+1)*class_size))
        teacher = n_students + (c % n_teachers)
        for stu in class_members:
            teach_layer.append([stu, teacher])

    # --- Common area layer ---
    common_layer = []
    for _ in range(500):
        i = np.random.randint(0, total_pop)
        j = np.random.randint(0, total_pop)
        if i != j:
            common_layer.append([i, j])

    # --- Random background layer ---
    random_layer = []
    for _ in range(300):
        i = np.random.randint(0, total_pop)
        j = np.random.randint(0, total_pop)
        if i != j:
            random_layer.append([i, j])

    # Assign layers AFTER initialize()
    sim.people.contacts['dorm'] = np.array(dorm_layer)
    sim.people.contacts['class'] = np.array(class_layer)
    sim.people.contacts['teacher'] = np.array(teacher_layer)
    sim.people.contacts['teach'] = np.array(teach_layer)
    sim.people.contacts['common'] = np.array(common_layer)
    sim.people.contacts['random'] = np.array(random_layer)

    # ---------------------------
    # 2. Interventions
    # ---------------------------

    interventions = []

    # Mask mandate in fall 2020
    interventions.append(cv.change_beta(days='2020-09-01', changes=0.7))

    # Winter break reduces class/common contacts
    interventions.append(cv.clip_edges(
        days=['2020-12-20', '2021-12-20'],
        changes=[0.2, 0.2],
        layers=['class', 'common']
    ))

    # No outside infections
    sim.pars['importation_rate'] = 0.00001

    sim.update_pars(interventions=interventions)

    return sim


if __name__ == '__main__':
    sim = make_boarding_school_sim()
    sim.run()

    # Plot and save
    fig = sim.plot()
    fig.savefig("boarding_school_simulation.png", dpi=300, bbox_inches='tight')

    print("Saved: boarding_school_simulation.png")
