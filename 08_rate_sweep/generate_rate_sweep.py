#!/usr/bin/env python3
"""
generate_rate_sweep.py  --  CORRECTED rate (Damkohler) sweep

Rebuilds the injection-rate sweep using the MANUSCRIPT configuration, after the
discovery that the study-02 generator (`01_baseline/generate_decks.py`) had
drifted badly from the paper:

    parameter        manuscript (Table 1/2)      drifted generate_decks.py
    domain           2000 m graded (140x1x25)    500 m uniform (50x1x25)
    injection fluid  pH 4.5, DIC 0.05 mol/L       pH 3.8, HCO3- 1.0 mol/L
    carbonate seed   1e-4 VF, 5000 m^2/m^3        0.0 VF, 100 m^2/m^3
    forsterite       Rimstidt -6.05/-4.65         drifted -10.0

The drifted decks injected a 20x-too-acidic, unseeded fluid into a 500 m box
that vented ~95% of the CO2 -- so dissolved carbonated zero and the sweep was
meaningless. This generator instead drives the manuscript-consistent dape
builder (verified against the kinetic-crossover / suite / study-07 runs),
extended to all six injection scenarios, scaling only the injection rate.

6 scenarios x 5 rate multipliers (0.3, 1, 3, 10, 30) = 30 decks.

Usage:
    python3 generate_rate_sweep.py
Output:
    ./decks/rs_<scenario>_mu<tag>.in
Requires ../03_dape/generate_dape_decks.py (the rate_mult + 6-scenario build_deck).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "03_dape"))
from generate_dape_decks import build_deck  # noqa: E402

SCENARIOS = [("dissolved", "S1"), ("scco2", "S2"), ("wag6mo", "S3"),
             ("wag3mo", "S4"), ("swag", "S5"), ("adaptive", "S6")]
RATE_MULTS = [("0p3", 0.3), ("1", 1.0), ("3", 3.0), ("10", 10.0), ("30", 30.0)]


def main():
    out = HERE / "decks"
    out.mkdir(exist_ok=True)
    print("=" * 68)
    print("  CORRECTED rate sweep  --  manuscript config (2000 m, pH 4.5/0.05,")
    print("  seeded carbonates, Rimstidt forsterite). 6 scenarios x 5 rates.")
    print("=" * 68)
    n = 0
    for sc, sid in SCENARIOS:
        for rt, rm in RATE_MULTS:
            deck = build_deck(sc, kappa=1.0, diffusion=1.0e-9, rate_mult=rm)
            (out / f"rs_{sc}_mu{rt}.in").write_text(deck)
            n += 1
        print(f"  {sid} {sc:10s}: mu = 0.3, 1, 3, 10, 30")
    print(f"\n{n} decks written to {out}/")
    print("Submit with ./submit_rate_sweep.sh ; analyse with analyse_rate_sweep.py")


if __name__ == "__main__":
    main()
