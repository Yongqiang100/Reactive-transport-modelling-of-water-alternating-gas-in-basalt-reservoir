# SUPCRT-HPT parallel case

A complete, ready-to-run parallel case mirroring every core simulation in the WAG study but
using the **supcrt-hpt.dat** thermodynamic database instead of hanford.dat, for a direct
Hanford-vs-SUPCRT-HPT comparison. Only the DATABASE line changes -- geochemistry inputs, grid,
injection, kinetics and mineral surface areas are identical -- so any difference is
attributable to the thermodynamic dataset alone.

**Motivation.** Under Hanford, secondary diopside precipitates from the CO2-supersaturated
fluid and intercepts most of the released Ca (and ~1/3 of the Mg), suppressing calcite.
SUPCRT-HPT is a Helgeson-lineage high-P/T dataset with better-constrained silicate equilibria;
running the same reversible-kinetics model on it tests whether diopside stays undersaturated
and dissolves on its own. This is a thermodynamic (data) correction -- nature deciding on
better data -- not an imposed dissolution-only assignment.

## Contents
- `make_supcrt_case.py`     generate the SUPCRT decks from your existing Hanford decks
- `get_supcrt_database.sh`  obtain supcrt-hpt.dat (PFLOTRAN install, or Bitbucket download)
- `submit_supcrt_all.sh`    submit every study as a Slurm array (smoke-test gated)
- `cation_balance.py`       per-mineral change + Mg/Ca budget for ANY case
- `compare_cases.py`        Hanford vs SUPCRT-HPT side-by-side (carbonate, diopside, calcite)

## Prerequisites
Your existing repo at `$MYSCRATCH/WAG`, containing the study folders (01_baseline, 03_dape,
..., 11_kappa_mu30), `run_study_setonix.sh`, and `make_manuscript_figures.py`.

## Setup
    cd $MYSCRATCH/WAG
    unzip supcrt_hpt_case.zip
    cp supcrt_hpt_case/*.py supcrt_hpt_case/*.sh .
    chmod +x *.sh

## Workflow
**1. Generate the parallel decks** (auto-discovers all core studies, incl. 09 and 11):

    python3 make_supcrt_case.py
    # -> supcrt_case/<study>/decks/*.in   (DATABASE supcrt-hpt.dat, chemistry identical)

**2. Get the database:**

    ./get_supcrt_database.sh
    # checks your PFLOTRAN install first; else downloads from Bitbucket (login node)

**3. SMOKE-TEST one deck** before committing ~147 runs. This checks BOTH species-name
compatibility and the diopside behaviour:

    sbatch --job-name=supcrt_test --array=0-0 --ntasks=16 --time=01:30:00 \
      --export=ALL,DECKS_DIR=$PWD/supcrt_case/01_baseline/decks,\
    RUN_ROOT=$PWD/supcrt_case/01_baseline/runs,HANFORD_DB=$PWD/supcrt-hpt.dat \
      run_study_setonix.sh

If it stops with **"... not found in database"**, a mineral/species is named differently in
SUPCRT-HPT. Add the mapping to the `RENAME` dict at the top of `make_supcrt_case.py`, re-run
`python3 make_supcrt_case.py`, and the whole case regenerates with corrected names.

**4. Run the full case** (only after the smoke test clears):

    ./submit_supcrt_all.sh $PWD/supcrt-hpt.dat

**5. Analyse and compare to Hanford:**

    python3 cation_balance.py                        # Hanford baseline (default)
    python3 cation_balance.py supcrt_case/01_baseline  # SUPCRT baseline
    python3 compare_cases.py                         # Hanford vs SUPCRT side-by-side

## What to look for
- Does diopside now **dissolve** (dn < 0) under SUPCRT-HPT instead of precipitate? If so, the
  supersaturation was a Hanford artifact and the Ca closes into calcite + dolomite.
- Do carbonate volumes **rise** (freed Ca/Mg)? Expected if diopside stops competing.
- Does the **S1 > WAG > S2 ordering hold**? It should -- it is a transport / phase-partitioning
  result, not a database one.

## Notes
- ~147 runs is a real allocation cost. Consider starting with 01_baseline + 08_rate_sweep to
  confirm the database behaves sensibly, then submit the rest.
- `supcrt-hpt.dat` is not bundled (obtain via step 2); it is too large to ship inline.
- All scripts assume they sit in `$MYSCRATCH/WAG` alongside the study folders.
