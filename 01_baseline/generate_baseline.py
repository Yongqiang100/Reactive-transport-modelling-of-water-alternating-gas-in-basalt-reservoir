#!/usr/bin/env python3
"""
generate_baseline.py  --  base-rate (Table-4) scenarios from the CANONICAL builder

Replaces the earlier standalone generate_decks.py, which had drifted from the
manuscript (500 m uniform grid, pH 3.8 / 1 M HCO3- injection, unseeded
carbonates). The six base-rate scenarios are now produced by the single
manuscript-consistent builder in ../03_dape/generate_dape_decks.py at
rate_mult = 1.0, guaranteeing identical geochemistry, grid (2000 m graded
140x1x25), injection chemistry (pH 4.5 / DIC 0.05), and seeded carbonates
across the entire study.

These six decks are identical to the mu=1 slice of study 08 (the rate
sweep); they are emitted here as a stand-alone baseline study for clarity
and to reproduce Table 4 / the baseline figures.

Usage:
    python3 generate_baseline.py        ->  ./decks/base_<scenario>.in
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "03_dape"))
from generate_dape_decks import build_deck  # noqa: E402

SCENARIOS = [("dissolved", "S1"), ("scco2", "S2"), ("wag6mo", "S3"),
             ("wag3mo", "S4"), ("swag", "S5"), ("adaptive", "S6")]


def main():
    out = HERE / "decks"
    out.mkdir(exist_ok=True)
    print("Base-rate scenarios (rate_mult=1.0) from the canonical builder:")
    for sc, sid in SCENARIOS:
        deck = build_deck(sc, kappa=1.0, diffusion=1.0e-9, rate_mult=1.0)
        (out / f"base_{sc}.in").write_text(deck)
        print(f"  {sid}  base_{sc}.in")
    print(f"\n6 baseline decks -> {out}/")
    print("Submit with ./submit_baseline.sh ; figures via analysis_baseline/.")


if __name__ == "__main__":
    main()
