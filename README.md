# jhk_covasim4wastewater — intervention_search

A Covasim-based epidemic surveillance simulation project covering:
1. Diagnosis/sequencing observation model (how much, and how late, do we observe true infections)
2. Diagnosis-driven intervention triggers (when to fire an intervention)
3. Reconstructing transmission/mutation trees from the infection log and comparing them phylogenetically (how well can sequence data alone recover the true transmission structure)

## Directory structure

```
intervention_search/
├── configs/
│   ├── sim_pars.yaml              # Base covasim Sim parameters (pop_size, n_days, beta, etc.)
│   ├── evo_pars.yaml              # Genomic evolution extension parameters (reference sequence, molecular clock, substitution model)
│   └── intervention_policies.yaml # Definitions for baseline / policy_A / policy_B / policy_C
│
├── simulation/
│   ├── observation_model.py       # ObservationModel — diagnosis/sequencing sampling
│   ├── intervention.py            # InterventionController — beta-reduction intervention
│   ├── trigger_functions.py       # Trigger rules + check_triggers dispatcher
│   └── runner.py                  # RunConfig / RunResults / Runner — simulation orchestrator
│
├── analysis/
│   ├── mutation_tree.py           # infection log -> transmission/mutation tree, FASTA/Newick export, IQ-TREE
│   └── metrics.py                 # detection/intervention timing, outbreak impact, ensemble summary metrics
│
├── plotting/
│   ├── plot_mutation_tree.py      # tree visualizations (event tree / pruned tree / phylogenetic tree)
│   └── run_plots.py               # single-run timeline, ensemble boxplots
│
└── cophylogeny/
    └── polytomy_tanglegram.py     # polytomy-aware tanglegram — side-by-side comparison of two trees
```

## File-by-file description

### configs/

- **sim_pars.yaml** — Base parameters passed straight into covasim's `Sim`. `pop_size`, `n_days`, `beta`, `pop_infected`, `rand_seed`, etc.
- **evo_pars.yaml** — Settings for covasim's genomic evolution extension. `enable: true`, reference sequence path (`reference`), molecular clock rate (`mol_clock_rate`), substitution model (`sub_model`), `store_infection_log: true` (keep the full infection log).
- **intervention_policies.yaml** — Definitions for the 4 policies. Each policy is a **plain dict** carrying test/sequencing rates (`p_test`, `p_seq`), intervention strength (`beta_reduction`) and duration (`duration`), and trigger conditions (`trigger_list`) — there's no separate Policy class; `ObservationModel` / `InterventionController` / `Runner` consume this dict directly.

### simulation/

- **observation_model.py** — `ObservationModel` class. Each day, draws the number of diagnosed cases from a Binomial(`p_test`) on true infections, spreads diagnosis dates forward via a Multinomial(`delay_pmf`), then draws how many of those get sequenced via a Binomial(`p_seq`). Also records which specific infection events were sequenced (`daily_sequenced_agents`).
- **intervention.py** — `InterventionController` class. When triggered, multiplies `beta` by `(1 - beta_reduction)`, holds it there for `duration` days, then leaves it as-is (no restoration).
- **trigger_functions.py** — Four trigger rules (`weekly_growth`, `slope_trigger`, `sustained_increase`, `relative_threshold`) plus `check_triggers()`, which walks `trigger_list` and dispatches by name.
- **runner.py** — The core engine of the project.
  - `RunConfig` — bundles simulation parameters, policy, and a run ID.
  - `RunResults` — daily records (infections/diagnoses/sequences/beta), the GT trigger (based on true infections) and D trigger (based on observed diagnoses), detection/intervention day, and GT→D delay statistics.
  - `Runner` — steps the covasim `Sim` day by day, applies `ObservationModel`, checks the GT and D triggers separately, and starts the intervention on a D trigger (subject to `quiet_period` and `max_interventions`).
  - `run_single()` / `run_evo_simulation()` — convenience wrappers around `Runner`. `run_evo_simulation()` takes `evo_pars` to run the genomic-evolution simulation, and is what the notebooks actually call.

### analysis/

- **mutation_tree.py**
  - `build_mutation_tree(sim, daily_sequenced_agents)` — a tree where each node is one infection event.
  - `build_infection_seq_tree(sim, daily_sequenced_agents)` — a tree collapsed down to unique mutation states (haplotypes).
  - `extract_sequenced_subtree(G)` — prunes branches that don't lead to a sequenced event.
  - `mutation_tree_to_newick(G, root=...)` — converts to Newick format.
  - `add_internal_leaves(G)` — duplicates internal nodes as pendant leaves (so ancestral states can be compared in a tanglegram).
  - `export_fasta_from_G()` / `infer_tree_jc69()` — FASTA export + running IQ-TREE (JC69).
  - `set_all_branch_lengths_to_one()` — sets all branch lengths to 1 so only topology is compared.
- **metrics.py**
  - `summarize_outbreak(results, interv)` — total infections / peak infections / outbreak duration / infections after intervention (all computed from `RunResults.true_infections`).
  - `summarize_detection(results)` / `summarize_intervention(results)` — detection day, intervention day, GT→D delay, etc., read straight off `RunResults` (which already computes them).
  - `summarize_ensemble(runs, key)` — mean/median/std/quartiles for a given metric across multiple runs.

### plotting/

- **plot_mutation_tree.py** — `plot_mutation_tree` (event tree), `plot_erase_mutation_tree` (fades/removes branches that don't lead anywhere observable), `plot_phylo_tree` (IQ-TREE Newick result).
- **run_plots.py** — `plot_single_run` (one run's infection/diagnosis/sequence timeline plus GT/D triggers and intervention timing), `plot_ensemble_boxplots` / `grouped_boxplot` (comparing metrics across runs/policies).

### cophylogeny/

- **polytomy_tanglegram.py** — parses two Newick trees, supports multi-child (polytomy) nodes, minimizes crossings, and draws them side by side. A single call to `run_polytomy_tanglegram(nwk1, nwk2)` handles parsing, optimization, and rendering.

## Policy dict field reference

| Key | Used by | Description |
|---|---|---|
| `p_test`, `p_seq`, `delay_pmf` | `ObservationModel` | test rate, sequencing rate, diagnosis delay distribution |
| `beta_reduction`, `beta_min`, `duration` | `InterventionController` | beta reduction fraction on intervention, floor beta, how many days it holds |
| `trigger_list`, `quiet_period`, `max_interventions` | `Runner` | list of trigger rules, cooldown before re-triggering, max times an intervention can fire per policy (default 1) |

## Weekly notebooks

| Notebook | What it demonstrates |
|---|---|
| **W4 — Observation Model** | How much the `ObservationModel` "hides" the true infection curve. Varies `p_test` and `delay_pmf` to see how far the observed diagnosis curve drifts from the true one. |
| **W5 — Intervention Threshold Pipeline** | A grid search over trigger types/values x intervention strengths/durations to find the best policy. Along the way, discovers that repeated interventions unrealistically drive beta to zero and always "win" -- this motivated later adding the `max_interventions` cap. |
| **W6 — Evolutionary Simulation Analysis** | Runs an ensemble (multiple seeds) of the baseline/A/B/C policies with `evo_pars` enabled, comparing detection delay, total infections, diagnoses, and sequences across policies via boxplots. Also visualizes the true transmission tree for one example run. |
| **W7 — Mutation Tree Inference** | Reconstructs the full genome of each sequenced infection event, exports FASTA, and infers a phylogenetic tree with IQ-TREE (JC69). Checks how well pure sequence data recovers the true transmission structure. |
| **W8 — Mutation Tree Topology** | Introduces a tree collapsed by mutation state (haplotype) instead of by event, and turns the FASTA-export/IQ-TREE pipeline into reusable functions. |
| **W9 — Polytomy Tanglegram** | Compares the "true transmission tree" against the "IQ-TREE-inferred tree" side by side using a polytomy-aware tanglegram, with crossing minimization to keep the comparison legible. |

## Design notes

- Policies are intentionally kept as plain dicts (no Policy class, no automatic grid-search generator).
- Variant-prevalence tracking and a practicality-vs-effectiveness Pareto analysis are out of scope for now.
- `Runner` always computes both the ground-truth trigger (GT, on true infections) and the diagnosis-based trigger (D, on observed diagnoses) side by side -- "how late, and how distorted, is our observation of the truth" is the question running through the whole project.