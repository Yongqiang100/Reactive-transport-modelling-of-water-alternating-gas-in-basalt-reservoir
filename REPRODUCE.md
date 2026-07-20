# Reproducing the results

This guide walks through reproducing every simulation, figure, and derived
number in the paper from this repository. There are two routes:

- **Route A — full re-run:** regenerate the input decks, run all PFLOTRAN
  simulations on an HPC cluster, then post-process and plot. This reproduces
  everything from scratch (~104 simulations).
- **Route B — figures only:** download the archived simulation outputs from
  Zenodo and run only the analysis/figure scripts. Use this if you want the
  figures and tables without the (substantial) cost of re-running PFLOTRAN.

The scientific model is identical in both routes; only whether you regenerate
the raw outputs differs.

---

## 1. What you need

**PFLOTRAN** — the reactive-transport solver. This study used the **v6 build**
run as a compiled MPI binary on Setonix (Pawsey). It is open source
(https://www.pflotran.org, source at https://bitbucket.org/pflotran/pflotran).
Only the compiled binary is needed to run the simulations; no Python is
involved at solve time.

**A thermodynamic database** — every deck reads `hanford.dat` via a plain
`DATABASE hanford.dat` line. Each run directory gets a symlink to a single copy,
pointed to by `HANFORD_DB` (see `setonix_env.sh`).

**Python (conda)** — used only to (a) generate the input decks and (b) read the
HDF5/observation output and make figures. The deck generators are pure standard
library; the analysis/figure scripts need `numpy`, `scipy`, `matplotlib`,
`h5py`, and `pandas`. Create the environment with:

```
conda env create -f environment.yml     # creates the env, then:
conda activate <env-name>                # e.g. co2conv on Setonix
```

**A Slurm cluster** — the 104 simulations are submitted as Slurm job arrays
(one array per study, one task per deck). Most studies finish within ~1.5 h on
16 MPI ranks; the high-rate WAG runs in `08_rate_sweep` are the most computationally demanding and use
more ranks and a longer walltime by default.

---

## 2. Route A — full re-run

### Step 0 — prepare and configure

Put the repository on scratch and check the paths in **`setonix_env.sh`**. This
file only *defines* variables (no module loads, no conda activation), so it is
safe to source anywhere. On Setonix the PFLOTRAN paths are already set; the one
value to confirm is where you keep the database:

```
PFLOTRAN_SOFT = /software/.../pflotran-v6          # your PFLOTRAN build
PFLOTRAN_ENV  = $PFLOTRAN_SOFT/env.sh              # runtime libraries
PFLOTRAN_EXE  = $PFLOTRAN_SOFT/pflotran/src/pflotran/pflotran
HANFORD_DB    = $MYSCRATCH/WAG/hanford.dat         # the database copy
```

Every value is overridable from the environment, so on a different machine you
can export new paths instead of editing the file.

### Step 1 & 2 — generate decks and submit all jobs

A single driver regenerates every deck and submits every study:

```
cd $MYSCRATCH/WAG
chmod +x run_all_setonix.sh run_study_setonix.sh
./run_all_setonix.sh
```

This activates the Python env, runs each set's `generate_*.py` (writing
`<set>/decks/*.in`), then submits one Slurm **job array** per study via
`run_study_setonix.sh` (one array task per deck). Useful options (all via
environment variables):

```
./run_all_setonix.sh --gen-only          # only write the decks; print the sbatch commands
NTASKS_PER_RUN=32 ./run_all_setonix.sh   # MPI ranks per deck (default 16)
WALLTIME=00:45:00 ./run_all_setonix.sh   # per-task walltime
ARRAY_THROTTLE=8  ./run_all_setonix.sh   # cap concurrent array tasks per study
```

Each array task runs one deck in its own run directory (`<set>/runs/<label>/`),
symlinks the database in, and executes a single
`srun -N 1 -n <ranks> $PFLOTRAN_EXE -input_prefix run`.

### Step 3 — monitor and resume

```
squeue -u $USER
```

Re-running the driver is safe: each task **skips a deck that already completed**
(it looks for `Wall Clock Time` in the run log), so if some jobs time out you can
resubmit, and only the unfinished runs execute.

### Step 4 — generate all figures

Once the runs are complete, produce the full figure set with the single entry
point (run it from the folder that holds the analysis modules and the study
directories):

```
cd $MYSCRATCH/WAG
python3 generate_all_figures.py          # -> ./figures/  (PDF + PNG)
```

`generate_all_figures.py` reads the PFLOTRAN outputs, performs the analysis
(reusing `make_manuscript_figures.py` and `analyse_transport_limitation.py`),
and writes every figure to `./figures/`. Set `WAG_ROOT` if your runs live in a
different tree.

For the set-specific tables and diagnostics, each study also has its own
post-processor, e.g.:

```
( cd 03_dape           && python3 analyse_dape_sweep.py )
( cd 05_kinetic_crossover && python3 analyse_kappa_sweep.py )
( cd 08_rate_sweep     && python3 analyse_rate_sweep.py )
```

---

## 3. Route B — figures only (from archived outputs)

If you do not want to re-run PFLOTRAN:

1. Download the simulation-output archive from **Zenodo**
   (https://doi.org/10.5281/zenodo.21390676; a single ZIP, ~33.9 GB, containing
   the full study) and unpack it so each `<set>/runs/<label>/` directory
   contains its PFLOTRAN output (`.h5`, observation, and mass-balance files).
2. Create and activate the Python environment (Section 1).
3. Run `python3 generate_all_figures.py` and the per-set `analyse_*.py`
   scripts exactly as in Step 4 above.

The analysis scripts only read the outputs, so this reproduces every figure and
number without a cluster.

---

## 4. Which command makes which figure

`generate_all_figures.py` writes all of these to `./figures/`:

| Figure (output name)        | Content                                   | Source routine |
|-----------------------------|-------------------------------------------|----------------|
| `fig_domain_schematic`      | model domain / setup                      | `make_manuscript_figures.fig_domain` |
| `fig_study_design`          | study-design schematic (no data)          | inline |
| `fig_baseline_comparison`   | baseline S1–S6 comparison                 | `make_manuscript_figures.fig_baseline` |
| `fig_carbonate_breakdown`   | carbonate mineral breakdown               | `make_manuscript_figures.fig_carb_bar` |
| `fig_spatial_profiles`      | near-well spatial profiles                | `make_manuscript_figures.fig_spatial` |
| `fig_gas_saturation_2d`     | 2-D gas-saturation fields                 | `make_manuscript_figures.fig_gas2d` |
| `fig_kappa_crossover`       | kinetic-sensitivity (κ) crossover         | `analyse_transport_limitation` |
| `fig_damkohler_sweep`       | injection-rate / Damköhler sweep          | `make_manuscript_figures.fig_da_sweep` |
| `fig_da_sigma_regime`       | Damköhler–Σ regime map                     | inline |

---

## 5. Study / directory map

Each `NN_*` folder is one set: a `generate_*.py` (writes `decks/`), the
runs it produces (`runs/`), and an `analyse_*.py`.

```
01_baseline            6 sims  — S1–S6 injection strategies at the baseline rate
03_dape               18 sims  — kinetic-rate and diffusion sweeps (Da–Pe separation)
04_mechanisms          6 sims  — phase-partitioning and buoyancy-by-position cases
05_kinetic_crossover  28 sims  — global kinetic-rate (κ) sweep, 7 orders of magnitude
06_grid_resolution     5 sims  — grid-convergence check
07_da_consistency     11 sims  — single-Damköhler collapse test
08_rate_sweep         30 sims  — injection-rate (Damköhler) sweep, 6 strategies × 5 rates
                     -------
                     104 simulations
```

Supporting sets (`09_scco2_kappa_controls`, `10_batch_conservation`,
`11_kappa_mu30`, `12_boundary_tests`, `supcrt_standalone`) provide the control
and consistency checks referenced in the text; each follows the same
`generate_* → run → analyse_*` pattern.

Top-level analysis modules: `make_manuscript_figures.py`,
`analyse_transport_limitation.py`, `generate_all_figures.py` (entry point),
`extract_two_timepoints.py`, `compare_to_paper.py`, `selfseal_check.py`.

---

## 6. Reproducing on a non-Setonix cluster

The physics and decks are cluster-independent; only the submission layer is
Setonix-specific. To port:

1. In `setonix_env.sh`, export `PFLOTRAN_EXE` (and `PFLOTRAN_ENV` if your build
   needs a runtime setup) and `HANFORD_DB` to your paths.
2. In `run_study_setonix.sh`, change the Slurm directives at the top
   (`--account`, `--partition`, node/rank/walltime) to your scheduler's values.
3. If you don't use Slurm, generate the decks with `--gen-only`, then run each
   `<set>/decks/*.in` with your PFLOTRAN
   (`mpirun -n <ranks> pflotran -input_prefix <deck-without-.in>`) in a run
   directory containing a `hanford.dat` symlink, and proceed to Step 4.

Runtime guide: on 16 ranks most studies complete in under ~1.5 h each;
`08_rate_sweep` is the exception (fast WAG slug transitions force small
timesteps) and benefits from 32 ranks and a multi-hour walltime.
