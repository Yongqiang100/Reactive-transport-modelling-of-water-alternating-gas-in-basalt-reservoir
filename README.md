# Transport and Gas-Water Phase Partitioning Govern CO₂ Mineralization in Basaltic Reservoirs

Reproducibility package — reactive-transport simulations (PFLOTRAN) and analysis for the manuscript *Transport and Gas-Water Phase Partitioning Govern CO₂ Mineralization in Basaltic Reservoirs* submitted to *JGR: Solid Earth*.

This repository contains everything needed to regenerate every input deck, run
every simulation, and reproduce every figure in the article, as well as the full
set of supporting **diagnostic and verification** results. Two external
dependencies are required and not bundled (they are not our code): the
**PFLOTRAN v6** solver and the **`hanford.dat`** thermodynamic database — see
[Prerequisites](#2-prerequisites).

---

## 1. What this package produces

```
generators (*.py)  ──►  input decks (decks/*.in)
                              │  PFLOTRAN v6 (Slurm job arrays)
                              ▼
                        raw outputs (runs/<label>/*.h5, *-mas.dat, *-obs*.tec)
                              │  analysis / figure scripts (*.py)
                              ▼
              manuscript figures (figures/*.pdf,*.png)  +  diagnostic tables
```

- **153 input decks** across **11 sets**; **104 simulations** define the
  seven core studies of the paper, with a further set of supporting/diagnostic
  runs.
- **9 manuscript figures**, produced by a single entry point
  (`generate_all_figures.py`).
- **Diagnostic results**: grid convergence, single-Damköhler collapse,
  boundary-condition sensitivity, kinetic controls at high rate, and
  single-cell carbon-conservation checks.

> **Reproducibility note — carbonate seeding.** All decks use the corrected
> **zero-seed** initial condition (`seed_carb_vf = 0.0d0`): the secondary
> carbonates start at zero volume fraction, so reported carbonate reflects only
> injection-driven precipitation. (An earlier configuration seeded them at
> `1.0d-4`, which added injection-independent apparent carbonate; that is
> superseded.) You can confirm the whole tree is on the corrected seed with
> `bash scan_seed.sh .` (see [§8](#8-verification--integrity)).

---

## 2. Prerequisites

Prepared once; see `BUILD_AND_RUN.md` for step-by-step build instructions.

- **PFLOTRAN v6** — compiled MPI solver. Use your cluster's module
  (`module load pflotran/6.0`) or build PETSc (3.18–3.21, with HDF5) + PFLOTRAN
  from source (https://www.pflotran.org, https://petsc.org). Build the **same
  version** used here (v6); a different build can shift the numbers.
- **`hanford.dat`** — thermodynamic database read by every deck
  (`DATABASE hanford.dat`). Use the copy distributed with PFLOTRAN
  (`.../pflotran/database/hanford.dat`) if unmodified, or the copy archived on
  Zenodo if it was customised. Point `HANFORD_DB` at it.
- **Python (conda)** — for writing decks and making figures only, not for
  solving. `conda env create -f environment.yml` then activate it
  (`numpy`, `scipy`, `matplotlib`, `h5py`, `pandas`). Deck generators are pure
  standard library.
- **A Slurm cluster** — simulations run as job arrays. Most studies finish in
  ≲1.5 h on 16 MPI ranks; `08_rate_sweep` is the exception (small timesteps
  during fast WAG slugs) and uses more ranks / longer walltime.

---

## 3. Quick start

```bash
# 0. configure paths for your machine
#    edit setonix_env.sh -> PFLOTRAN_EXE, HANFORD_DB
#    (non-Setonix) edit the SBATCH directives in run_study_setonix.sh

# 1. generate all decks + submit all simulations as Slurm job arrays
./run_all_setonix.sh

# 2. once the runs finish, produce every manuscript figure
python3 generate_all_figures.py          # -> ./figures/  (PDF + PNG)
```

To skip re-running the solver, download the archived outputs from Zenodo (see
[§9](#9-data-availability)) into `runs/` and run only step 2 and the analysis
scripts below.

---

## 4. Repository layout

```
setonix_env.sh              path/exports (PFLOTRAN_EXE, HANFORD_DB, ...)
run_all_setonix.sh          driver: regenerate all decks + submit every study
run_study_setonix.sh        per-deck Slurm array task (one deck per array index)
environment.yml             Python environment for decks + figures
generate_all_figures.py     ONE entry point -> all 9 manuscript figures
make_manuscript_figures.py  figure/reader library (imported by the above)
analyse_transport_limitation.py  kinetic-crossover analysis + readers

NN_*/                       one set each:
   generate_*.py  (or make_*.py)   builds this set's decks/*.in
   decks/*.in                       the committed input decks
   analyse_*.py                     set-specific post-processing
   runs/<label>/                    created at run time (excluded from the repo)

# diagnostic / analysis scripts (top level)
cation_balance.py  gas_buoyancy.py  nearwell_precip.py  selfseal_check.py
extract_two_timepoints.py  compare_to_paper.py

# verification tools
compare_code.py  verify_decks.sh  verify_package.sh  scan_seed.sh

# documentation
README.md  BUILD_AND_RUN.md  REPRODUCE.md  PUBLICATION_CHECKLIST.md
```

The shared deck builder lives in `03_dape/generate_dape_decks.py`
(`build_deck`); the other set generators import it, so decks must be
generated within the full tree layout (the drivers and `verify_decks.sh` handle
this automatically).

---

## 5. Reproducing the manuscript results

### 5.1 The study map (104 core simulations)

| Set | Sims | What it varies | Role in the paper |
|---|---:|---|---|
| `01_baseline` | 6 | Six fluid-delivery schemes (S1 dissolved, S2 scCO₂, S3 WAG-6mo, S4 WAG-3mo, S5 SWAG, S6 adaptive) at the baseline rate | Baseline comparison; phase-partitioning result |
| `03_dape` | 18 | Kinetic-rate κ (×5) and molecular diffusion D (×4), each × {dissolved, scCO₂} | Separates Damköhler from Péclet |
| `04_mechanisms` | 6 | CO₂ mole fraction (caseC ×4) and injection elevation (caseD ×2) | Phase-partitioning / buoyancy mechanism |
| `05_kinetic_crossover` | 28 | Global κ over 7 orders (10⁻⁵–10²), × {dissolved, scCO₂} | Kinetic insensitivity / transport limitation |
| `06_grid_resolution` | 5 | Mesh resolution | Grid-convergence check |
| `07_da_consistency` | 11 | Rate q and κ chosen to hold Damköhler constant | Single-Damköhler collapse test |
| `08_rate_sweep` | 30 | Injection rate μ (×5) × six schemes | Rate/Damköhler sweep; self-sealing |

The production grid is **250 × 1 × 50** (refined from the manuscript's
descriptive 140 × 1 × 25; `06_grid_resolution` demonstrates convergence).
Study: 30-yr injection + 70-yr post-injection monitoring = 100 yr.

### 5.2 Run the simulations

```bash
./run_all_setonix.sh
```

This regenerates every deck and submits one Slurm **job array per study** (one
array task per deck). Useful options (environment variables):

```bash
./run_all_setonix.sh --gen-only          # only write decks + print the sbatch lines
NTASKS_PER_RUN=32 ./run_all_setonix.sh   # MPI ranks per deck (default 16)
WALLTIME=00:45:00 ./run_all_setonix.sh   # per-task walltime
ARRAY_THROTTLE=8  ./run_all_setonix.sh   # cap concurrent tasks per study
```

Monitor with `squeue -u $USER`. **Re-running is safe**: a task whose log already
contains `Wall Clock Time` is skipped, so if some jobs time out, resubmit.

### 5.3 Generate the figures

```bash
python3 generate_all_figures.py          # run from the repo root -> ./figures/
```

| Figure (output name)       | Source routine                                | Reads |
|----------------------------|-----------------------------------------------|-------|
| `fig_domain_schematic`     | `make_manuscript_figures.fig_domain`          | — (schematic) |
| `fig_study_design`         | inline in `generate_all_figures.py`           | — (schematic) |
| `fig_baseline_comparison`  | `make_manuscript_figures.fig_baseline`        | `01_baseline` |
| `fig_carbonate_breakdown`  | `make_manuscript_figures.fig_carb_bar`        | `01_baseline` |
| `fig_spatial_profiles`     | `make_manuscript_figures.fig_spatial`         | `01_baseline` |
| `fig_gas_saturation_2d`    | `make_manuscript_figures.fig_gas2d`           | `01_baseline` |
| `fig_kappa_crossover`      | `analyse_transport_limitation.analyse_05` + `fig_crossover` | `05_kinetic_crossover` |
| `fig_damkohler_sweep`      | `make_manuscript_figures.fig_da_sweep`        | `08_rate_sweep` |
| `fig_da_sigma_regime`      | inline in `generate_all_figures.py`           | `08_rate_sweep` |

Set `WAG_ROOT` if the runs live in a different tree. Set-specific tables and
diagnostics come from each study's own post-processor, e.g.
`( cd 08_rate_sweep && python3 analyse_rate_sweep.py )`.

---

## 6. Reproducing the diagnostic & supporting results

These sets and scripts support the paper's claims and are the internal
checks referenced in the text / supporting information. Run any one after its
simulations complete.

### 6.1 Supporting simulation sets

| Set | Sims | Purpose | Post-processor |
|---|---:|---|---|
| `09_scco2_kappa_controls` | 10 | scCO₂ kinetic controls (`sk_ctrl_*`, `sk_inj_*`) — confirms the κ result under control conditions | `analyse_transport_limitation.analyse_09` |
| `11_kappa_mu30` | 28 | The κ sweep repeated at 30× injection rate — kinetic insensitivity persists off-baseline | `11_kappa_mu30/kappa_mu30_analysis.py` |
| `12_boundary_tests` | 6 | Open-top and tall-caprock boundary variants (`open_top_*`, `tallcap_*`) × {dissolved, scCO₂, WAG-6mo} — boundary-condition insensitivity | `12_boundary_tests/compare_boundary.py` |

Run a single supporting set as its own array (six decks → `--array=0-5`):

```bash
sbatch --job-name=wag_12_boundary_tests --array=0-5 --ntasks=16 --time=01:30:00 \
  --export=ALL,DECKS_DIR=$PWD/12_boundary_tests/decks,RUN_ROOT=$PWD/12_boundary_tests/runs \
  run_study_setonix.sh
```

### 6.2 Carbon-conservation batch checks (`10_batch_conservation`)

A **1 m³ closed reactor** with the identical chemistry / database / brine /
minerals as the full runs — no wells, no boundary flux, no injection — used to
verify carbon accounting in the geochemical step in isolation from flow and
transport. Five decks:

| Deck | Configuration |
|---|---|
| `batch_cell` | baseline GENERAL + GIRT closed cell |
| `batch_cell_nocarb` | carbonate precipitation disabled (rate constants → ~0): conservation control |
| `v1_noco2gas` | `CO2(g)` removed from `PASSIVE_GAS_SPECIES` |
| `v2_co2aq_primary` | `CO2(aq)` used as the primary carbon species (`HCO3⁻` secondary) |
| `v3_th_nogas` | single-phase `TH`, no gas component |

Build, run, and read:

```bash
( cd 10_batch_conservation && python3 make_batch_cell.py && \
  python3 make_batch_cell_nocarb.py && python3 make_batch_variants.py )   # writes decks/
sbatch --job-name=wag_batch --array=0-4 --ntasks=1 --time=00:20:00 \
  --export=ALL,DECKS_DIR=$PWD/10_batch_conservation/decks,RUN_ROOT=$PWD/10_batch_conservation/runs \
  run_study_setonix.sh
# read each run's four-phase carbon budget:
for r in batch_cell batch_cell_nocarb v1_noco2gas v2_co2aq_primary v3_th_nogas; do
  python3 10_batch_conservation/batch_carbon_budget.py 10_batch_conservation/runs/$r
done
```

`batch_carbon_budget.py` sums **all four carbon reservoirs** — aqueous DIC,
dissolved CO₂ (liquid "Air" component), free CO₂ gas, and mineral carbonate —
and prints a closed-system consistency check (conservative Cl⁻ tracer, cation
movement, silicate dissolution). The control variants (`batch_cell_nocarb`,
`v1`, `v2`, `v3`) isolate how carbon is partitioned among the gas, aqueous, and
reactive representations in a GENERAL + GIRT setup; the total carbon should be
compared across all four phases, not aqueous + mineral alone. These are
diagnostic checks (not manuscript figures).

### 6.3 Diagnostic analysis scripts

Run against completed `runs/` (no re-simulation needed):

- `cation_balance.py` — divalent-cation (Ca, Mg, Fe) release vs carbonate uptake.
- `gas_buoyancy.py` — vertical gas-saturation distribution / buoyancy diagnostics.
- `nearwell_precip.py` — near-well carbonate localisation (fraction within the first cells).
- `selfseal_check.py` — porosity/permeability reduction (self-sealing) at high rate.
- `extract_two_timepoints.py` — pulls fields at chosen times for spatial comparisons.
- `compare_to_paper.py` — cross-checks computed quantities against the reported values.

Each takes a run directory or set as its argument; run with `-h` or read the
header docstring for usage.

---

## 7. Running one set or one deck

Generate decks only, then run a single deck manually (no Slurm):

```bash
./run_all_setonix.sh --gen-only
cd 01_baseline/runs/base_dissolved        # create the dir; put a hanford.dat symlink in it
mpirun -n 8 $PFLOTRAN_EXE -input_prefix run
```

Each array task (via `run_study_setonix.sh`) runs one deck in
`<set>/runs/<label>/`, symlinks the database, and executes
`srun -N 1 -n <ranks> $PFLOTRAN_EXE -input_prefix run`.

---

## 8. Verification & integrity

The repository is self-checking. Keep `compare_code.py` beside the scripts.

```bash
# every committed deck is exactly what the generators produce, and the
# generators/source are unchanged (read-only):
bash verify_decks.sh .

# no deck or completed run uses the superseded 1.0d-4 carbonate seed:
bash scan_seed.sh .

# a built package reproduces this tree (source + regenerated decks):
bash verify_package.sh <package>.zip .
```

`verify_decks.sh` prints two verdicts — generator/source unchanged, and decks
reproduce — and passes only if both hold. `scan_seed.sh` should report
*"nothing on the old seed."*

---

## 9. Data availability

- **Code**: this repository
  (https://github.com/Yongqiang100/Reactive-transport-modelling-of-water-alternating-gas-in-basalt-reservoir),
  also archived on Zenodo (https://doi.org/10.5281/zenodo.21390676).
- **Simulation outputs** (HDF5 snapshots, observation and mass-balance files for
  every run) are archived on Zenodo at
  **https://doi.org/10.5281/zenodo.21390676** — a single ZIP archive (~33.9 GB
  compressed, ~65.8 GB extracted) containing the complete study: input decks,
  generators, run and analysis scripts, the full set of simulation outputs, and
  the figures. Unpack the outputs into the corresponding `runs/<label>/`
  directories to reproduce the figures without re-running PFLOTRAN.
- **PFLOTRAN**: v6, open source (https://www.pflotran.org).
- **Thermodynamic database**: `hanford.dat` (stock PFLOTRAN database; if a
  customised copy was used it is included in the Zenodo archive).

---

## 10. Runtime & resources

104 core simulations plus the supporting sets. On 16 MPI ranks most studies
complete in ≲1.5 h each; `08_rate_sweep` benefits from 32 ranks and a multi-hour
walltime. Budget several GB of disk per scenario for HDF5 snapshots. Deck
generation and figure production are seconds-to-minutes on a login node.

---

## 11. Citation

If you use this code or data, please cite both the paper and the Zenodo archive:

> Chen, Y., Xie, Q., Cao, X., Kang, Q., & Regenauer-Lieb, K. (in review).
> *Transport and Gas-Water Phase Partitioning Govern CO₂ Mineralization in
> Basaltic Reservoirs.* Journal of Geophysical Research: Solid Earth.
> [DOI to be added on acceptance]
>
> Chen, Y., Xie, Q., Cao, X., Kang, Q., & Regenauer-Lieb, K. (2026).
> *Reproducibility archive for "Transport and Gas-Water Phase Partitioning
> Govern CO₂ Mineralization in Basaltic Reservoirs"* [Data set]. Zenodo.
> https://doi.org/10.5281/zenodo.21390676

See `CITATION.cff` for machine-readable metadata.

---

## 12. License

Released under the MIT License (see `LICENSE`). PFLOTRAN and its thermodynamic
database are distributed under their own licenses.

---

## 13. Contact

Yongqiang Chen — yongqiang.chen@curtin.edu.au (corresponding author).
