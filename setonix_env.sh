#!/bin/bash
# setonix_env.sh — paths for the WAG CO2-mineralization package on Setonix
# (Pawsey). This file only DEFINES variables (no conda activation, no module
# loads, no sourcing of other envs), so it is safe to source both on the login
# node and inside a Slurm job:
#   - run_all_setonix.sh sources this, then activates co2conv for deck generation.
#   - run_study_setonix.sh (the job) sources this, then sources $PFLOTRAN_ENV to
#     load the PFLOTRAN runtime (Cray libraries) before running the binary.
# Every value is overridable from the environment.

: "${MYSCRATCH:?MYSCRATCH is not set — are you on Setonix? (expected /scratch/<project>/<user>)}"

# --- PFLOTRAN v6 build (under /software, same as the H2 runs) -------------
# PFLOTRAN_ENV is the runtime setup the job sources so the binary's libraries
# resolve (this is the env your dfnworks_env.sh / run_chem_sweeps.sh use).
export PFLOTRAN_SOFT="${PFLOTRAN_SOFT:-/software/projects/pawsey1284/ychen6/pflotran-v6}"
export PFLOTRAN_ENV="${PFLOTRAN_ENV:-$PFLOTRAN_SOFT/env.sh}"
export PFLOTRAN_EXE="${PFLOTRAN_EXE:-$PFLOTRAN_SOFT/pflotran/src/pflotran/pflotran}"

# --- co2conv (Python for deck generation + analysis) ----------------------
export CONDA_BASE_SETONIX="${CONDA_BASE_SETONIX:-$MYSCRATCH/conda}"
export CO2CONV="${CO2CONV:-$MYSCRATCH/conda/envs/co2conv}"

# --- geochemical database referenced by every deck (DATABASE hanford.dat) -
# Each run directory gets a symlink to this, so the decks' bare
# `DATABASE hanford.dat` resolves.
export HANFORD_DB="${HANFORD_DB:-$MYSCRATCH/WAG/hanford.dat}"

echo "[setonix_env] PFLOTRAN_EXE = $PFLOTRAN_EXE"
echo "[setonix_env] PFLOTRAN_ENV = $PFLOTRAN_ENV"
echo "[setonix_env] HANFORD_DB   = $HANFORD_DB"
