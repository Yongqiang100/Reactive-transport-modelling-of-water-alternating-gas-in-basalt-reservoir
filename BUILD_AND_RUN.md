# Build and run — step by step

This repository contains the input decks, the run drivers, and the
analysis/figure code. Two things are prepared **once** and live outside the
repository — the **PFLOTRAN** solver and the **`hanford.dat`** database (see
Prerequisites). After that, the steps below regenerate every simulation and
figure in the paper.

> **Version note.** The results were produced with **PFLOTRAN v6**. Use the same
> version — a different PFLOTRAN/PETSc build can shift the numbers slightly.

---

## Prerequisites (prepare once)

- **Simulator — PFLOTRAN (v6).** A compiled Fortran/MPI code built on PETSc. If
  your cluster provides it, `module load pflotran/6.0`; otherwise build PETSc
  (3.18–3.21, with HDF5) then PFLOTRAN from source following the official docs
  (https://www.pflotran.org, https://petsc.org). Note the path to the `pflotran`
  executable — you set it in Step 2.
- **Python environment.** `conda env create -f environment.yml`, then activate it
  (e.g. `co2conv`). Used only to write the decks and make the figures
  (`numpy`, `scipy`, `matplotlib`, `h5py`, `pandas`) — not for solving.
- **Thermodynamic database — `hanford.dat`.** Every deck reads
  `DATABASE hanford.dat`. Use the copy distributed with PFLOTRAN
  (`.../pflotran/database/hanford.dat`) if unmodified, or the copy archived with
  the data if it was customised. You point `HANFORD_DB` at it in Step 2.

---

## Step 1 — Get the code

```bash
git clone https://github.com/Yongqiang100/Reactive-transport-modelling-of-water-alternating-gas-in-basalt-reservoir.git
cd Reactive-transport-modelling-of-water-alternating-gas-in-basalt-reservoir   # repository root
```

On an HPC system, put this on scratch (fast parallel filesystem), not `$HOME`.

---

## Step 2 — Point the run scripts at your setup

Edit **`setonix_env.sh`** so these resolve on your system (each is also
overridable by exporting it in your shell):

```bash
PFLOTRAN_EXE = <path to the pflotran executable>
PFLOTRAN_ENV = <optional: a script that sets runtime libraries for your build>
HANFORD_DB   = <path to hanford.dat>
```

If you are **not** on Setonix, also edit the Slurm directives at the top of
**`run_study_setonix.sh`** (`--account`, `--partition`, nodes/ranks/walltime).

---

## Step 3 — Run the 104 simulations

One driver regenerates every deck and submits every study as a Slurm job array
(one array task per deck):

```bash
chmod +x run_all_setonix.sh run_study_setonix.sh
./run_all_setonix.sh
```

Useful options (all via environment variables):

```bash
./run_all_setonix.sh --gen-only          # only write decks + print the sbatch lines
NTASKS_PER_RUN=32 ./run_all_setonix.sh   # MPI ranks per deck (default 16)
WALLTIME=00:45:00 ./run_all_setonix.sh   # per-task walltime
ARRAY_THROTTLE=8  ./run_all_setonix.sh   # cap concurrent tasks per study
```

Monitor with `squeue -u $USER`. **Re-running is safe**: a task that already
finished (its log contains `Wall Clock Time`) is skipped, so if some jobs time
out, resubmit; only the unfinished runs execute. Each task runs one deck
in `<set>/runs/<label>/`, symlinks the database, and executes
`srun -N 1 -n <ranks> $PFLOTRAN_EXE -input_prefix run`.

**No Slurm / a single run on a workstation:** generate the decks with
`./run_all_setonix.sh --gen-only`, then run any deck directly:
```bash
cd 01_baseline/runs/base_dissolved      # create the dir; put a hanford.dat symlink in it
mpirun -n 8 $PFLOTRAN_EXE -input_prefix run
```

---

## Step 4 — Post-process and make the figures

Once the runs are done, produce the full figure set with one command (run it
from the repository root):

```bash
python3 generate_all_figures.py          # -> ./figures/  (PDF + PNG)
```

It reads the PFLOTRAN outputs, runs the analysis, and writes every manuscript
figure to `./figures/`. Set `WAG_ROOT` if the runs live elsewhere. For
set-specific tables/diagnostics, each study also has its own
post-processor, e.g. `( cd 08_rate_sweep && python3 analyse_rate_sweep.py )`.

---

## Runtime and resources

104 simulations across seven core studies (`01_baseline`, `03_dape`,
`04_mechanisms`, `05_kinetic_crossover`, `06_grid_resolution`,
`07_da_consistency`, `08_rate_sweep`), plus supporting control sets
(`09`–`12`, `supcrt_standalone`). On 16 MPI ranks most studies finish within
~1.5 h each; `08_rate_sweep` is the exception (fast WAG slug transitions force
small timesteps) and benefits from 32 ranks and a multi-hour walltime. Budget a
few GB of disk per scenario for the HDF5 snapshots.

---

## Troubleshooting

- **`Cannot open file hanford.dat`** — set `HANFORD_DB` (Prerequisites / Step 2);
  confirm a `hanford.dat` symlink appears in the run directory.
- **PFLOTRAN executable not found** — set `PFLOTRAN_EXE` in `setonix_env.sh`.
- **Jobs rejected by the scheduler** — the `--account`/`--partition` in
  `run_study_setonix.sh` are Setonix-specific; change them for your cluster.
- **Figure script errors on missing files** — a run hasn't finished; check
  `squeue` and resubmit (completed runs are skipped).
- **Numbers differ slightly from the paper** — confirm the same PFLOTRAN version
  and the same database.
