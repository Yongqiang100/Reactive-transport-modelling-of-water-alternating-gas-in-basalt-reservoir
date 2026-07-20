#!/usr/bin/env python3
"""
generate_scco2_kappa_controls.py  —  control-subtracted scCO2 kinetic sweep

Closes the one gap in the transport-limitation evidence: study 05 shows the
scCO2 endmember has a larger apparent kinetic response than dissolved (a ~2.3x
low->high-kappa step), but study 05 uses (final - initial) carbonate WITHOUT a
no-injection control, so part of that step is background re-equilibration of the
initial formation water (which is itself kappa-dependent). Study 07 proved, for
DISSOLVED, that subtracting a no-injection control flattens the curve almost
completely (control-subtracted injection-driven carbonate is flat to ~0.1%).

This study provides the matching no-injection controls for scCO2 so the
injection-driven scCO2 mineralization can be isolated as a function of kappa:

    injection_driven(kappa) = (VF_inj(final) - VF_ctrl(final)) integrated.

Design (mirrors study 07, scenario = scco2, q = base):
  kappa ladder spanning the study-05 transition (low plateau, soft crossover,
  high plateau):  1e-2, 1e-1, 3e-1, 1, 10
  - injection runs  : rate_mult = 1.0   (sk_inj_<ktag>)
  - no-injection ctrl: rate_mult = 0.0   (sk_ctrl_<ktag>)  one per kappa
Total: 5 + 5 = 10 decks. scCO2 is stiff -> run like 08 (32 ranks, long walltime).

Usage:
    python3 generate_scco2_kappa_controls.py
Output:
    ./decks/sk_inj_<ktag>.in , ./decks/sk_ctrl_<ktag>.in
"""
import sys, math
from pathlib import Path

HERE = Path(__file__).resolve().parent
DAPE_DIR = HERE.parent / "03_dape"
sys.path.insert(0, str(DAPE_DIR))
from generate_dape_decks import build_deck  # noqa: E402  (verified deck builder)

KAPPAS = [1e-2, 1e-1, 3e-1, 1.0, 10.0]


def ktag(k: float) -> str:
    exp = math.floor(math.log10(k))
    mant = round(k / 10 ** exp)
    if mant == 10:
        mant, exp = 1, exp + 1
    sign = "m" if exp < 0 else ""
    return f"k{mant}e{sign}{abs(exp)}"


def main():
    out = HERE / "decks"
    out.mkdir(exist_ok=True)
    print("=" * 70)
    print("  Control-subtracted scCO2 kinetic sweep")
    print(f"  kappa ladder: {KAPPAS}")
    print("  injection (rate_mult=1) + no-injection control (rate_mult=0) per kappa")
    print("=" * 70)
    n = 0
    for k in KAPPAS:
        for tag, rm in [("inj", 1.0), ("ctrl", 0.0)]:
            deck = build_deck("scco2", kappa=k, diffusion=1.0e-9, rate_mult=rm)
            (out / f"sk_{tag}_{ktag(k)}.in").write_text(deck)
            n += 1
            print(f"  OK sk_{tag}_{ktag(k):6s}  scco2  kappa={k:<6g} rate_mult={rm}")
    print(f"\n{n} decks written to {out}/")
    print("Run (scCO2 is stiff -> 32 ranks, long walltime):")
    print("  sbatch --job-name=wag_09_scco2k --array=0-9 --ntasks=32 --time=12:00:00 \\")
    print("    --export=ALL,DECKS_DIR=$PWD/09_scco2_kappa_controls/decks,RUN_ROOT=$PWD/09_scco2_kappa_controls/runs \\")
    print("    run_study_setonix.sh")
    print("Then: python3 analyse_transport_limitation.py --only 09")


if __name__ == "__main__":
    main()
