#!/bin/bash
# check_runs.sh — report which PFLOTRAN simulations completed.
# A run is COMPLETED when its run.log contains "Wall Clock Time" (PFLOTRAN's
# normal-completion marker). Run from the package root after jobs finish:
#     ./check_runs.sh        summary per study + list of anything not done
#     ./check_runs.sh -v     also list every completed run
set -u
ROOT="$(cd "$(dirname "$0")" && pwd)"
VERBOSE=0; [ "${1:-}" = "-v" ] && VERBOSE=1

STUDIES=(01_baseline 03_dape 04_mechanisms 05_kinetic_crossover \
         06_grid_resolution 07_da_consistency 08_rate_sweep 09_scco2_kappa_controls)

tot_exp=0; tot_ok=0; tot_bad=0; tot_miss=0
problems=()

echo "=========================================================="
echo "  WAG CO2 run status  ($(date))"
echo "  COMPLETED = run.log contains 'Wall Clock Time'"
echo "=========================================================="
for dir in "${STUDIES[@]}"; do
    decks_dir="$ROOT/$dir/decks"; runs_dir="$ROOT/$dir/runs"
    if [ ! -d "$decks_dir" ] || [ -z "$(ls -A "$decks_dir"/*.in 2>/dev/null)" ]; then
        printf "  %-22s no decks (generate first)\n" "$dir"; continue
    fi
    exp=0; ok=0; bad=0; miss=0
    for deck in "$decks_dir"/*.in; do
        [ -e "$deck" ] || continue
        exp=$((exp+1)); name="$(basename "${deck%.in}")"; log="$runs_dir/$name/run.log"
        if [ ! -f "$log" ]; then
            miss=$((miss+1)); problems+=("$dir/$name : NOT STARTED (no run.log)")
        elif grep -q "Wall Clock Time" "$log" 2>/dev/null; then
            ok=$((ok+1)); [ "$VERBOSE" -eq 1 ] && echo "  OK    $dir/$name"
        else
            bad=$((bad+1))
            hint="incomplete (timed out / crashed / still running)"
            grep -qiE "error|abort|fault|killed|out of memory|oom" "$log" 2>/dev/null && hint="ERROR in log"
            last="$(tail -n 1 "$log" 2>/dev/null | cut -c1-70)"
            problems+=("$dir/$name : $hint  [last: $last]")
        fi
    done
    printf "  %-22s %2d/%-2d done" "$dir" "$ok" "$exp"
    [ "$bad"  -gt 0 ] && printf "   %d incomplete" "$bad"
    [ "$miss" -gt 0 ] && printf "   %d not-started" "$miss"
    echo
    tot_exp=$((tot_exp+exp)); tot_ok=$((tot_ok+ok)); tot_bad=$((tot_bad+bad)); tot_miss=$((tot_miss+miss))
done

echo "----------------------------------------------------------"
echo "  TOTAL: $tot_ok / $tot_exp completed   ($tot_bad incomplete, $tot_miss not-started)"
if [ "${#problems[@]}" -eq 0 ] && [ "$tot_exp" -gt 0 ]; then
    echo "  All simulations completed successfully."
    exit 0
fi
if [ "${#problems[@]}" -gt 0 ]; then
    echo ""
    echo "  Needs attention (${#problems[@]}):"
    printf "    - %s\n" "${problems[@]}"
    echo ""
    echo "  To finish them, just re-run the launcher — the resume guard skips"
    echo "  everything already completed, so only the unfinished decks run:"
    echo "      ./run_all_setonix.sh"
    echo "  Or resubmit a single study's array (rate sweep, 30 decks, shown):"
    echo "    sbatch --job-name=wag_08_rate_sweep --array=0-29 --ntasks=16 --time=01:30:00 \\"
    echo "      --export=ALL,DECKS_DIR=\$PWD/08_rate_sweep/decks,RUN_ROOT=\$PWD/08_rate_sweep/runs \\"
    echo "      run_study_setonix.sh"
    echo ""
    echo "  Diagnose a specific failure:  tail -50 <study>/runs/<name>/run.log"
    echo "  Slurm-side (timeouts/OOM):     sacct -X --name=wag_<study> --format=JobID,State,Elapsed,MaxRSS"
fi
exit 1
