#!/usr/bin/env python3
"""
generate_grid_decks.py  —  grid-resolution & front-structure verification

Strict test of two claims that the baseline 2 m near-well mesh cannot settle
on its own:

  (1) the Pe-invariance of integrated carbonation is grid-converged, and
  (2) whether the Peclet number controls the WIDTH (sharpness) of the
      carbonate front, as opposed to its integrated magnitude.

The baseline near-well cell is 2 m, while the molecular diffusive front width
sqrt(2 D t) over 30 yr is ~1.4 m at baseline D (0.7 cell, unresolved) rising
to ~14 m at 100x D (7 cells, resolved). To resolve the baseline front the
near-well zone (0-100 m) is refined from 2 m down to 0.5 m, leaving the
far-field unchanged so the domain and the advective field are identical.

Endmember: S1 dissolved only. Single-phase aqueous injection gives a clean
horizontal advective front with no buoyant gas cap to confound the width
measurement.

Runs (kappa = 1 throughout; only grid and D vary):
  convergence (D = 1e-9):  near-well 2 m, 1 m, 0.5 m
  Pe sub-sweep (0.5 m):    D = 1e-9 (1x), 1e-8 (10x), 1e-7 (100x)
  -> 5 unique decks (0.5 m / 1e-9 is shared between the two sets).

The deck builder is imported from the Da-Pe generator; only the grid (nxyz,
dxyz) and diffusion differ from Suite B, so this is a controlled extension.

Usage:
    python3 generate_grid_decks.py
Output:
    ./decks/grid_<cell>_D<mult>.in
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DAPE_DIR = HERE.parent / "03_dape"
sys.path.insert(0, str(DAPE_DIR))
from generate_dape_decks import build_deck  # noqa: E402  (verified builder)

# Near-well refinement of the 0-100 m zone at the PRODUCTION far-field grading
# and z (50@3 + 50@9 + 50@26 in x; 50@2 in z). Each x-spec sums to 2000 m, so the
# domain and advection are unchanged. The "1m" level IS the production grid
# (canonical builder default); 2m brackets it coarser, 0.5m finer, demonstrating
# the production mesh is grid-converged.
GRIDS = {
    "2m":   ("200 1 50", "    50@2.d0 50@3.d0 50@9.d0 50@26.d0\n    1.d0\n    50@2.d0"),
    "1m":   ("250 1 50", "    100@1.d0 50@3.d0 50@9.d0 50@26.d0\n    1.d0\n    50@2.d0"),
    "0p5m": ("350 1 50", "    200@0.5d0 50@3.d0 50@9.d0 50@26.d0\n    1.d0\n    50@2.d0"),
}

# Diffusion multiples relative to the 1e-9 baseline (Pe ~ 1e7, 1e6, 1e5)
DMULT = {"1x": 1.0e-9, "10x": 1.0e-8, "100x": 1.0e-7}

# (grid_key, D_key): convergence at baseline D, then Pe sub-sweep at 0.5 m
RUNS = [
    ("2m",   "1x"),    # convergence anchor (reproduces Suite B dissolved baseline)
    ("1m",   "1x"),    # convergence
    ("0p5m", "1x"),    # convergence finest + Pe baseline
    ("0p5m", "10x"),   # Pe x10  at finest grid
    ("0p5m", "100x"),  # Pe x100 at finest grid
]


def main():
    out = HERE / "decks"
    out.mkdir(exist_ok=True)

    print("=" * 66)
    print("  Grid-resolution & front-structure verification (S1 dissolved)")
    print("  near-well refinement 2 m -> 1 m -> 0.5 m;  D = 1e-9, 1e-8, 1e-7")
    print("=" * 66)

    seen = set()
    n = 0
    for gkey, dkey in RUNS:
        if (gkey, dkey) in seen:
            continue
        seen.add((gkey, dkey))
        nxyz, dxyz = GRIDS[gkey]
        D = DMULT[dkey]
        deck = build_deck("dissolved", kappa=1.0, diffusion=D,
                          nxyz=nxyz, dxyz=dxyz)
        fname = out / f"grid_{gkey}_D{dkey}.in"
        fname.write_text(deck)
        n += 1
        nx = nxyz.split()[0]
        print(f"  ✓ {fname.name:24s} (NXYZ x={nx:>3}, D={D:.0e})")

    print(f"\n{n} decks written to {out}/")
    print("Convergence set : grid_2m_D1x, grid_1m_D1x, grid_0p5m_D1x")
    print("Pe sub-sweep    : grid_0p5m_D1x, grid_0p5m_D10x, grid_0p5m_D100x")


if __name__ == "__main__":
    main()
