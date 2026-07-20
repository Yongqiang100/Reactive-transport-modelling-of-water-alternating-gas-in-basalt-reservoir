#!/usr/bin/env bash
# =====================================================================
# check_dape_jobs.sh — Status report for the Da-Pe sweep
#
# Reads submitted_jobs.txt, queries Slurm for each job's state,
# and reports which ones completed/failed/are still running.
# =====================================================================

BASE_DIR="${BASE_DIR:-$HOME/WAG/DaPe-disentangling}"
SUMMARY="$BASE_DIR/submitted_jobs.txt"

if [ ! -f "$SUMMARY" ]; then
    echo "ERROR: $SUMMARY not found. Did you run submit_dape_sweep.sh yet?"
    exit 1
fi

echo "============================================================"
echo "  Da-Pe SWEEP — Job Status Report"
echo "  $(date)"
echo "============================================================"
printf "%-12s  %-40s  %s\n" "JobID" "Name" "State"
echo "------------------------------------------------------------"

total=0; done=0; failed=0; running=0; pending=0; other=0

while read -r line; do
    [[ "$line" =~ ^#.*$ ]] && continue
    [[ -z "$line" ]] && continue
    JID=$(echo "$line" | awk '{print $1}')
    NAME=$(echo "$line" | awk '{print $2}')
    [ -z "$JID" ] && continue
    total=$((total+1))

    # Query Slurm (sacct is more reliable for finished jobs)
    STATE=$(sacct -j "$JID" --format=State --noheader 2>/dev/null \
              | head -1 | awk '{print $1}' | sed 's/[+]$//')
    if [ -z "$STATE" ]; then
        STATE="UNKNOWN"
    fi

    case "$STATE" in
        COMPLETED)  done=$((done+1)) ;;
        FAILED|TIMEOUT|CANCELLED|NODE_FAIL) failed=$((failed+1)) ;;
        RUNNING)    running=$((running+1)) ;;
        PENDING)    pending=$((pending+1)) ;;
        *)          other=$((other+1)) ;;
    esac

    printf "%-12s  %-40s  %s\n" "$JID" "$NAME" "$STATE"
done < "$SUMMARY"

echo "------------------------------------------------------------"
echo "  Total: $total  |  Done: $done  |  Running: $running"
echo "         Pending: $pending  |  Failed: $failed  |  Other: $other"
echo "============================================================"

if [ "$failed" -gt 0 ]; then
    echo ""
    echo "FAILED jobs — check logs in their working directories:"
    while read -r line; do
        [[ "$line" =~ ^#.*$ ]] && continue
        JID=$(echo "$line" | awk '{print $1}')
        NAME=$(echo "$line" | awk '{print $2}')
        [ -z "$JID" ] && continue
        STATE=$(sacct -j "$JID" --format=State --noheader 2>/dev/null \
                  | head -1 | awk '{print $1}' | sed 's/[+]$//')
        case "$STATE" in
            FAILED|TIMEOUT|CANCELLED|NODE_FAIL)
                # Find the working dir
                if [[ "$NAME" == suiteA_* ]]; then
                    SUITE="suiteA"
                else
                    SUITE="suiteB"
                fi
                SUB="${NAME/suiteA_/}"
                SUB="${SUB/suiteB_/}"
                echo "  $NAME  →  $BASE_DIR/$SUITE/$SUB/${NAME}-${JID}.out"
                ;;
        esac
    done < "$SUMMARY"
fi

if [ "$done" -eq "$total" ] && [ "$total" -gt 0 ]; then
    echo ""
    echo "✓ All jobs complete. Run:"
    echo "    python3 analyse_dape_sweep.py"
fi
