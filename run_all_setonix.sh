#!/bin/bash
# run_all_setonix.sh — full WAG CO2 reproduction on Setonix (Pawsey).
# Generates every deck under co2conv, then submits ONE Slurm job ARRAY per study
# (one array task per deck) via run_study_setonix.sh. Each task is its own
# allocation with a single srun — no concurrent steps in one allocation.
# 104 simulations across 7 studies.
#
#   cd $MYSCRATCH/WAG
#   chmod +x run_all_setonix.sh run_study_setonix.sh
#   ./run_all_setonix.sh
#
# Options (environment):
#   NTASKS_PER_RUN=32 ./run_all_setonix.sh   # MPI ranks per deck (default 16; 2-D grid ~12.5k cells)
#   WALLTIME=00:45:00 ./run_all_setonix.sh   # per-task walltime (default 01:30:00)
#   ARRAY_THROTTLE=8  ./run_all_setonix.sh   # cap concurrent array tasks per study (adds %N)
#   ./run_all_setonix.sh --gen-only          # generate decks, print sbatch cmds, do not submit

set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
source "$ROOT/setonix_env.sh"

# co2conv for the Python deck generators (stdlib-only, so non-fatal).
if [ -f "$CONDA_BASE_SETONIX/bin/activate" ]; then
    # shellcheck disable=SC1091
    source "$CONDA_BASE_SETONIX/bin/activate" "$CO2CONV" 2>/dev/null \
        || echo "WARN: could not activate co2conv ($CO2CONV); using $(command -v python3)"
else
    echo "WARN: conda not found at $CONDA_BASE_SETONIX; using $(command -v python3) for deck generation"
fi

GEN_ONLY=0
[ "${1:-}" = "--gen-only" ] && GEN_ONLY=1

THROTTLE=""
[ -n "${ARRAY_THROTTLE:-}" ] && THROTTLE="%${ARRAY_THROTTLE}"

# Per-study resources. The high-rate WAG/adaptive runs in 08_rate_sweep are by
# far the stiffest (tiny timesteps through the near-well cells during fast WAG
# slug transitions at mu=10-30), so they get more ranks and a longer walltime by
# default. The other studies all complete inside 1.5 h on 16 ranks. Setting
# WALLTIME / NTASKS_PER_RUN in the environment overrides these for ALL studies.
study_walltime() { case "$1" in 08_rate_sweep) echo "${WALLTIME:-12:00:00}";; *) echo "${WALLTIME:-01:30:00}";; esac; }
study_ranks()    { case "$1" in 08_rate_sweep) echo "${NTASKS_PER_RUN:-32}";;    *) echo "${NTASKS_PER_RUN:-16}";;    esac; }

# study : generator script   (decks land in <study>/decks/)
STUDIES=(
  "01_baseline:generate_baseline.py"
  "03_dape:generate_dape_decks.py"
  "04_mechanisms:generate_mechanism_decks.py"
  "05_kinetic_crossover:generate_kappa_sweep.py"
  "06_grid_resolution:generate_grid_decks.py"
  "07_da_consistency:generate_da_consistency_decks.py"
  "08_rate_sweep:generate_rate_sweep.py"
)

echo ">>> [1/2] Generating decks ..."
total=0
for entry in "${STUDIES[@]}"; do
    IFS=":" read -r dir gen <<< "$entry"
    ( cd "$ROOT/$dir" && python3 "$gen" >/dev/null )
    n=$(ls "$ROOT/$dir"/decks/*.in 2>/dev/null | wc -l)
    printf "    %-22s %s decks\n" "$dir" "$n"
    total=$((total + n))
done
echo "    TOTAL decks: $total"

echo ""
echo ">>> [2/2] Submitting one job ARRAY per study (per-study ranks/walltime; 08_rate_sweep heavier) ..."
for entry in "${STUDIES[@]}"; do
    IFS=":" read -r dir gen <<< "$entry"
    n=$(ls "$ROOT/$dir"/decks/*.in 2>/dev/null | wc -l)
    [ "$n" -gt 0 ] || { echo "    $dir: no decks, skipping"; continue; }
    JOB="wag_${dir}"
    ARRAY="0-$((n-1))${THROTTLE}"
    WALL="$(study_walltime "$dir")"
    NT="$(study_ranks "$dir")"
    EXPORTS="ALL,DECKS_DIR=$ROOT/$dir/decks,RUN_ROOT=$ROOT/$dir/runs"
    if [ "$GEN_ONLY" -eq 1 ]; then
        echo "    sbatch --job-name=$JOB --array=$ARRAY --ntasks=$NT --time=$WALL --export=$EXPORTS run_study_setonix.sh"
    else
        sbatch --job-name="$JOB" --array="$ARRAY" --ntasks="$NT" \
               --time="$WALL" --export="$EXPORTS" "$ROOT/run_study_setonix.sh"
    fi
done

echo ""
if [ "$GEN_ONLY" -eq 1 ]; then
    echo "Gen-only: decks built, submission commands printed above (nothing submitted)."
else
    echo "All studies submitted as job arrays. Monitor:  squeue -u \$USER"
    echo "Check completion:  ./check_runs.sh"
fi
