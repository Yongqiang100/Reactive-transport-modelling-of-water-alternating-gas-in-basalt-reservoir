#!/bin/bash
# run_study_setonix.sh — run ONE PFLOTRAN deck per Slurm array task.
# Submitted as a job ARRAY by run_all_setonix.sh (one array per study, one task
# per deck): each task is its own allocation running a single `srun`, so there
# are NO concurrent job steps inside one allocation (which is what made the
# earlier "step creation still disabled / Job step already running" failures).
# This mirrors the working run_chem_sweeps.sh array pattern.
#
# Each array task:
#   - picks decks/<DECKS_DIR>/*.in[SLURM_ARRAY_TASK_ID]  (glob is sorted)
#   - skips it if already completed (run.log has "Wall Clock Time")
#   - else runs it on $SLURM_NTASKS ranks in its own run dir
#
# run_all_setonix.sh provides --array, --ntasks, --time, --job-name, --export.
# Manual single-study submit (30-deck study shown):
#   sbatch --job-name=wag_08_rate_sweep --array=0-29 --ntasks=16 --time=01:30:00 \
#     --export=ALL,DECKS_DIR=$PWD/08_rate_sweep/decks,RUN_ROOT=$PWD/08_rate_sweep/runs \
#     run_study_setonix.sh
#
#SBATCH --account=pawsey1284
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -u
set -o pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$ROOT"
source "$ROOT/setonix_env.sh"

# PFLOTRAN runtime (Cray libraries); the task runs only the compiled binary.
if [ -f "$PFLOTRAN_ENV" ]; then
    source "$PFLOTRAN_ENV"
else
    echo "WARNING: PFLOTRAN_ENV not found: $PFLOTRAN_ENV" >&2
fi
export PFLOTRAN_EXE="${PFLOTRAN_EXE:-$PFLOTRAN_SOFT/pflotran/src/pflotran/pflotran}"

NTASKS="${SLURM_NTASKS:-16}"
DECKS_DIR="${DECKS_DIR:?set DECKS_DIR to the decks directory of the study to run}"
RUN_ROOT="${RUN_ROOT:-$(cd "$DECKS_DIR/.." && pwd)/runs}"
DB_SRC="$HANFORD_DB"
DB_NAME="$(basename "$DB_SRC")"

[ -x "$PFLOTRAN_EXE" ] || { echo "ERROR: PFLOTRAN_EXE not executable: $PFLOTRAN_EXE" >&2; exit 1; }
[ -f "$DB_SRC" ]       || { echo "ERROR: database not found: $DB_SRC"            >&2; exit 1; }

# This task's deck (bash sorts the glob lexically; same order run_all counted).
shopt -s nullglob
DECKS=( "$DECKS_DIR"/*.in )
[ "${#DECKS[@]}" -gt 0 ] || { echo "ERROR: no .in decks in $DECKS_DIR" >&2; exit 1; }
idx="${SLURM_ARRAY_TASK_ID:-0}"
DECK="${DECKS[$idx]:-}"
[ -n "$DECK" ] || { echo "ERROR: no deck at array index $idx (have ${#DECKS[@]})" >&2; exit 1; }
name="$(basename "${DECK%.in}")"
RUN="$RUN_ROOT/$name"

echo "=========================================================="
echo "  study deck   : $name   (array index $idx)"
echo "  ranks        : $NTASKS   node: $(hostname)"
echo "  run dir      : $RUN"
echo "  PFLOTRAN     : $PFLOTRAN_EXE"
echo "=========================================================="

# resume guard: skip if this deck already finished
if [ -f "$RUN/run.log" ] && grep -q "Wall Clock Time" "$RUN/run.log" 2>/dev/null; then
    echo "[$name] already completed; skipping"
    exit 0
fi

mkdir -p "$RUN"
ln -sf "$DB_SRC" "$RUN/$DB_NAME"
cp "$DECK" "$RUN/run.in"

cd "$RUN"
echo "[$name] launching PFLOTRAN on $NTASKS ranks  $(date +%T)"
srun -N 1 -n "$NTASKS" -c 1 "$PFLOTRAN_EXE" -input_prefix run > run.log 2>&1
rc=$?

if grep -q "Wall Clock Time" run.log 2>/dev/null; then
    echo "[$name] COMPLETED"
else
    echo "[$name] DID NOT COMPLETE (srun rc=$rc) — see $RUN/run.log" >&2
    exit 1
fi
