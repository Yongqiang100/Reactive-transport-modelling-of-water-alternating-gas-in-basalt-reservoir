# Running on Setonix (Pawsey)

This is the Setonix variant of the WAG CO2-mineralization reproduction package.
The **model is identical** to the main package (2-D, 250×1×50, 30-yr injection +
70-yr monitoring, six injection strategies + rate sweep); only the job-submission
layer is changed to Setonix.

Submission uses **Slurm job arrays** — one array per study, one task per deck —
the same pattern as your `run_chem_sweeps.sh`. Each task is its own allocation
running a single `srun -N 1 -n <ranks> -c 1 $PFLOTRAN_EXE -input_prefix run`, so
there are **no concurrent job steps inside one allocation** (that is what caused
the earlier "step creation still disabled / Job step already running" failures).
Slurm packs the small tasks onto nodes itself. Common settings:
`--account=pawsey1284`, `--partition=work`, per-run dir + `hanford.dat` symlink,
and a "Wall Clock Time" resume guard.

Two environments are kept separate, on purpose:

- **PFLOTRAN runtime** — each array task sources `$PFLOTRAN_ENV`
  (`…/pflotran-v6/env.sh`, the same runtime your H2 / chem-sweep jobs use) so the
  binary's Cray libraries load. The task runs only the compiled PFLOTRAN, so it
  does **not** touch conda.
- **Python (co2conv)** — deck generation (and analysis) on the login node uses the
  `co2conv` conda env. The deck generators are pure-stdlib, so this is non-fatal:
  if co2conv is absent, generation falls back to the system `python3`.

## 1. One-time setup

The PFLOTRAN paths are already set to your v6 build, so for the
`pawsey1284/ychen6` account there is usually **nothing to edit**:

```
PFLOTRAN_SOFT = /software/projects/pawsey1284/ychen6/pflotran-v6
PFLOTRAN_ENV  = $PFLOTRAN_SOFT/env.sh
PFLOTRAN_EXE  = $PFLOTRAN_SOFT/pflotran/src/pflotran/pflotran
HANFORD_DB    = $MYSCRATCH/WAG/hanford.dat
```

Stage the package on scratch and check the one path that depends on where you
keep the database:

```bash
cd $MYSCRATCH/WAG                      # stage the package here (on /scratch)
# only if your hanford.dat is elsewhere:
export HANFORD_DB=/path/to/hanford.dat
```

Every value above is overridable from the environment (e.g. `export
PFLOTRAN_SOFT=…` if you move the build, or `export CO2CONV=…` for a different
conda env). `setonix_env.sh` itself only *defines* these variables — it activates
no conda and sources no other env, so it is safe to source anywhere; the job adds
`source $PFLOTRAN_ENV` and the launcher adds the co2conv activation.

## 2. Run everything

```bash
cd $MYSCRATCH/WAG
chmod +x run_all_setonix.sh run_study_setonix.sh
./run_all_setonix.sh
```

This generates all 104 decks (co2conv), then submits **one job array per study**
(seven arrays; each task is one deck on `NTASKS_PER_RUN` cores). Monitor with
`squeue -u $USER`. To stage decks and see the submit plan without submitting:
`./run_all_setonix.sh --gen-only`.

**Check that everything finished.** A run is complete when its
`run.log` contains `Wall Clock Time` (PFLOTRAN's normal-completion marker).
After the jobs end:

```bash
chmod +x check_runs.sh && ./check_runs.sh        # per-study summary + anything not done (-v lists all)
```

It reports `TOTAL: N / 104 completed` and, for each run that did not finish,
whether it errored, timed out / is still running, or never started — plus the
resubmit command to finish them (the resume guard skips the completed ones). A
one-off manual equivalent:

```bash
grep -rl "Wall Clock Time" */runs/*/run.log | wc -l     # how many finished (want 104)
grep -rL "Wall Clock Time" */runs/*/run.log             # which run.logs are not finished
```

## 3. Cores per run

Each deck is one array task on `NTASKS_PER_RUN` MPI ranks (default **16**). These
are small 2-D grids (~12,500 cells), so 16 ranks (~780 cells/rank) is well sized;
Slurm packs several tasks onto each 128-core node. Tasks are **not** exclusive, so
the node fills with other tasks automatically.

```bash
NTASKS_PER_RUN=32 ./run_all_setonix.sh     # 32 ranks/deck if you want each run faster
```

Do **not** go to 128 ranks for these decks — that over-decomposes the solve
(~98 cells/rank) and wastes the node. (An earlier version packed 8 concurrent
`srun` runs inside one 128-core allocation; that triggered Slurm "step creation
disabled / step already running" errors on Setonix and is why this now uses job
arrays instead — one `srun` per allocation.)

If you want to limit how many tasks of a study run at once (e.g. to be gentle on
the allocation), throttle the array:

```bash
ARRAY_THROTTLE=8 ./run_all_setonix.sh      # at most 8 tasks per study run concurrently (adds %8)
```

## 4. Walltime and resuming

Each task (one deck) requests **1.5 h** (`--time=01:30:00`, the default in
`run_all_setonix.sh` and `run_study_setonix.sh`; the `work` partition allows up
to 24 h). Override with the `WALLTIME` knob, e.g. `WALLTIME=00:45:00
./run_all_setonix.sh`. Because each deck is its own task now, 1.5 h is per single
run (a 100-yr 2-D run on 16 ranks finishes well inside that).

To finish anything that failed or didn't run, **just re-run the launcher** — the
resume guard skips every deck that already completed, so only the unfinished ones
actually run:

```bash
./run_all_setonix.sh
```

Or resubmit a single study's array (rate sweep, 30 decks, shown):

```bash
sbatch --job-name=wag_08_rate_sweep --array=0-29 --ntasks=16 --time=01:30:00 \
       --export=ALL,DECKS_DIR=$PWD/08_rate_sweep/decks,RUN_ROOT=$PWD/08_rate_sweep/runs \
       run_study_setonix.sh
```

(The heaviest decks are the high-rate WAG/adaptive runs in `08_rate_sweep`: at
μ=10–30 the short WAG slugs force tiny timesteps and a single run can take several
hours. The launcher therefore gives `08_rate_sweep` **32 ranks and a 12 h
walltime** by default, while the other six studies stay at 16 ranks / 1.5 h. If an
08 run is killed by walltime, its `run.log` ends in `DUE TO TIME LIMIT` followed
by `MPI_Abort(…, 59)` — that is a walltime kill, **not** a solver crash; just
resubmit, and raise `WALLTIME` further if needed.)

## 5. Outputs and analysis

Each deck produces, in `<study>/runs/<deck_name>/`:
`run.h5`, `run-mas.dat`, `run.log` (and any `.vtk`/snapshot files the deck emits).
The output prefix is `run` (the deck is copied to `run.in` and run with
`-input_prefix run`).

**Compare to the paper.** Once the runs are complete, `compare_to_paper.py`
computes the headline quantities (baseline carbonate volume fractions, CO2 budget
+ mineral-storage efficiency, the rate-sweep domain-mean vs injection-driven
scaling, and grid convergence) and prints them next to the manuscript's claimed
values:

```bash
source $MYSCRATCH/conda/bin/activate $MYSCRATCH/conda/envs/co2conv   # provides h5py + numpy
python3 compare_to_paper.py --inventory   # FIRST: list run.h5 datasets + -mas.dat header
python3 compare_to_paper.py               # the paper-vs-computed table
```

(Note: `source setonix_env.sh` only sets the PFLOTRAN/database paths for the
*jobs* — it does **not** activate conda, so for analysis you activate co2conv as
above. If `import h5py` fails, add it with co2conv active: `pip install h5py`, or
`conda install -c conda-forge h5py`.)

Run `--inventory` first to confirm the dataset names on your PFLOTRAN build match
the reader (carbonate `*_VF`, `Porosity`, `pH`, and the `Air Mass` columns); then
the table. Note the manuscript numbers it prints are the **original**
30-yr/coarse-mesh values, so differences from the revised 100-yr/refined runs are
expected — they are what the manuscript update will incorporate.

The per-study `analyse_*.py` scripts give finer detail (Da–Pe 03, mechanisms 04,
kappa 05, Da-consistency 07); point them at the Setonix layout with
`BASE_DIR=<study>/runs` (each looks for its `<deck_name>/` run dirs and the `run`
prefix). I can wire those into `compare_to_paper.py` too if you want everything in
one place.

## What changed from the main (hpc01) package

- Removed `run_all.sh` and the per-study `submit_*.sh` (hpc01 / `module load
  pflotran/6.0` / `srun --mpi=pmix`).
- Added `setonix_env.sh`, `run_study_setonix.sh`, `run_all_setonix.sh`,
  `check_runs.sh`.
- Submission is a **Slurm job array per study** (one task per deck, a single
  `srun` each, `NTASKS_PER_RUN` ranks), with the resume guard. No concurrent job
  steps inside an allocation.
- Everything else (decks, generators, geochemistry, manuscript, docs) is
  unchanged.
