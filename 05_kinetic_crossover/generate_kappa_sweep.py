#!/usr/bin/env python3
"""
generate_kappa_sweep.py  —  global kinetic-rate crossover sweep

Maps the FULL transport-to-reaction transition along the kinetic axis by
scaling every mineral rate constant by a global multiplier kappa over a
wide range, at the base injection rate, for the two endmember
configurations (S1 dissolved, S2 supercritical).

This extends the Da-Pe Suite A kinetic sweep (kappa in [0.1, 10], which
lies entirely on the transport-limited plateau) downward into the
reaction-limited regime so that the complete carbonation(kappa) curve is
obtained:

    high kappa  ->  transport-limited plateau (carbonation ~ constant)
    crossover   ->  Da ~ O(1)
    low  kappa  ->  reaction-limited tail (carbonation ~ kappa)

Because only the rate constants change, the injected mass, mesh, chemistry,
and Peclet number are identical across the ladder — so the crossover is a
pure Damkohler effect with none of the total-mass confound carried by the
injection-rate sweep. The kappa in [0.1, 10] points reproduce Suite A and
serve as a built-in consistency check.

The deck builder is imported directly from the Da-Pe generator so the only
thing that differs from Suite A is the kappa ladder.

Usage:
    python3 generate_kappa_sweep.py
Output:
    ./decks/kappa_<scenario>_<ktag>.in    (2 x 14 = 28 decks)
"""

import sys
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
DAPE_DIR = HERE.parent / "03_dape"
sys.path.insert(0, str(DAPE_DIR))

# Reuse the EXACT verified deck builder from the Da-Pe study. build_deck
# copies BASE_KINETICS per call (no compounding), so repeated calls with
# different kappa are independent.
from generate_dape_decks import build_deck, BASE_KINETICS  # noqa: E402

# Wide kappa ladder: plateau (>=0.1) through the crossover into the
# reaction-limited tail (down to 1e-5). The 0.1-10 points overlap Suite A.
KAPPAS = [
    1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2,
    1e-1, 3e-1, 1.0, 3.0, 10.0, 30.0, 100.0,
]
SCENARIOS = ["dissolved", "scco2"]


def ktag(k: float) -> str:
    """Filesystem-safe scientific tag: 1e-5 -> 'k1em5', 0.3 -> 'k3em1',
    1.0 -> 'k1e0', 30 -> 'k3e1', 100 -> 'k1e2'."""
    exp = math.floor(math.log10(k))
    mant = round(k / 10 ** exp)
    if mant == 10:            # guard rounding at decade boundaries
        mant, exp = 1, exp + 1
    sign = "m" if exp < 0 else ""
    return f"k{mant}e{sign}{abs(exp)}"


def main():
    out = HERE / "decks"
    out.mkdir(exist_ok=True)

    print("=" * 66)
    print("  Global kinetic-rate crossover sweep")
    print(f"  kappa ladder ({len(KAPPAS)} values): "
          f"{KAPPAS[0]:g} ... {KAPPAS[-1]:g}")
    print(f"  scenarios: {SCENARIOS}")
    print(f"  baseline forsterite(alkaline) log_k = "
          f"{BASE_KINETICS['forsterite_alk']['log_k']} (kappa=1)")
    print("=" * 66)

    n = 0
    for scenario in SCENARIOS:
        for k in KAPPAS:
            deck = build_deck(scenario, kappa=k, diffusion=1.0e-9)
            fname = out / f"kappa_{scenario}_{ktag(k)}.in"
            fname.write_text(deck)
            n += 1
            shifted = BASE_KINETICS['forsterite_alk']['log_k'] + math.log10(k)
            print(f"  ✓ {fname.name:34s} (kappa={k:<8g} "
                  f"forsterite_alk log_k={shifted:+.2f})")

    print(f"\n{n} decks written to {out}/")
    print("Note: kappa in [0.1, 10] reproduce Da-Pe Suite A (consistency check).")


if __name__ == "__main__":
    main()
