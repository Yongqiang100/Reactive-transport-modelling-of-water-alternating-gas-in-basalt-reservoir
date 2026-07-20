#!/usr/bin/env python3
"""Kinetic-rate crossover sweep at 30x injection rate (mu = 30)."""
import sys, math
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "03_dape"))
from generate_dape_decks import build_deck, BASE_KINETICS
KAPPAS = [1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0, 3.0, 10.0, 30.0, 100.0]
SCENARIOS = ["dissolved", "scco2"]; RATE_MULT = 30.0
def ktag(k):
    exp = math.floor(math.log10(k)); mant = round(k / 10 ** exp)
    if mant == 10: mant, exp = 1, exp + 1
    return f"k{mant}e{'m' if exp<0 else ''}{abs(exp)}"
def main():
    out = HERE / "decks"; out.mkdir(parents=True, exist_ok=True); n = 0
    for sc in SCENARIOS:
        for k in KAPPAS:
            (out / f"kappa30_{sc}_{ktag(k)}.in").write_text(
                build_deck(sc, kappa=k, diffusion=1.0e-9, rate_mult=RATE_MULT)); n += 1
            print(f"  wrote kappa30_{sc}_{ktag(k)}.in  (kappa={k:g}, rate_mult={RATE_MULT:g})")
    print(f"\n{n} decks in {out}/  -> submit as a 0-{n-1} array")
if __name__ == "__main__": main()
