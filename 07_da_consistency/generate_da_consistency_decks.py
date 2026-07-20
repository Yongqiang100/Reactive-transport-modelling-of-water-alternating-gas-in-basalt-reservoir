#!/usr/bin/env python3
"""
generate_da_consistency_decks.py  —  Damkohler consistency (collapse) test

Responds to the reviewer concern that the Damkohler scaling may not be defined
or interpreted consistently: the manuscript reads the injection-rate sweep as a
transport-to-reaction-limited TRANSITION, yet the kinetic sweep shows almost no
transition. Those two are only consistent if a SINGLE Damkohler number governs
both. This study tests that directly.

Da is varied along two independent routes from the SAME deck builder, so the
only things that differ between runs are q and kappa:
    Da  ~  (reaction rate) / (advection rate)  ~  kappa / q
    -> rate axis  (kappa = 1): Da_rel = 1/q      [transport side, via residence]
    -> kinetic ax (q = 1):     Da_rel = kappa     [reaction side]

Matched-Da pairs (same Da_rel via the two routes) are the crux: if carbonation
efficiency is a function of Da alone, each pair must agree.
    q=3,kappa=1  <->  q=1,kappa=1/3     (Da_rel = 1/3)
    q=10,kappa=1 <->  q=1,kappa=1/10    (Da_rel = 1/10)
    q=30,kappa=1 <->  q=1,kappa=1/30    (Da_rel = 1/30)

No-injection controls (rate_mult=0), one per distinct kappa, isolate the
injection-driven carbonate exactly (final_injection - final_control per cell) --
robust even at high q where the plume travels far (a far-field baseline would
be contaminated there).

Endmember: S1 dissolved (clean single-phase advective transport).

Runs (11 total):
    rate axis    : da_q{1,3,10,30}_k1            (kappa=1)
    kinetic axis : da_q1_k{0p033,0p1,0p333}      (q=1; q1_k1 shared with rate ax)
    controls     : da_q0_k{1,0p333,0p1,0p033}    (rate_mult=0)

Usage:
    python3 generate_da_consistency_decks.py
Output:
    ./decks/da_q<rate>_k<kappa>.in
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DAPE_DIR = HERE.parent / "03_dape"
sys.path.insert(0, str(DAPE_DIR))
from generate_dape_decks import build_deck  # noqa: E402  (verified builder)

# (rate_tag, rate_mult, kappa_tag, kappa)
RATE_AXIS = [("1", 1.0), ("3", 3.0), ("10", 10.0), ("30", 30.0)]      # kappa = 1
KAPPA_AXIS = [("0p033", 1.0 / 30), ("0p1", 0.1), ("0p333", 1.0 / 3)]  # q = 1 (k1 shared)
CTRL_KAPPA = [("1", 1.0), ("0p333", 1.0 / 3), ("0p1", 0.1), ("0p033", 1.0 / 30)]


def main():
    out = HERE / "decks"
    out.mkdir(exist_ok=True)
    print("=" * 70)
    print("  Damkohler consistency (collapse) test  --  S1 dissolved")
    print("  Da_rel = kappa / q ; rate axis (k=1) and kinetic axis (q=1)")
    print("=" * 70)

    runs = []
    for rtag, rm in RATE_AXIS:                       # rate axis, kappa = 1
        runs.append((f"da_q{rtag}_k1", rm, 1.0))
    for ktag, kap in KAPPA_AXIS:                     # kinetic axis, q = 1
        runs.append((f"da_q1_k{ktag}", 1.0, kap))
    for ktag, kap in CTRL_KAPPA:                     # no-injection controls
        runs.append((f"da_q0_k{ktag}", 0.0, kap))

    seen = set()
    n = 0
    for name, rm, kap in runs:
        if name in seen:
            continue
        seen.add(name)
        deck = build_deck("dissolved", kappa=kap, diffusion=1.0e-9, rate_mult=rm)
        (out / f"{name}.in").write_text(deck)
        n += 1
        da_rel = (kap / rm) if rm > 0 else float("inf")
        da_s = "control" if rm == 0 else f"Da_rel={da_rel:.3f}"
        print(f"  OK {name:18s}  rate_mult={rm:<4g} kappa={kap:.4g}   {da_s}")

    print(f"\n{n} decks written to {out}/")
    print("Matched-Da pairs to compare:")
    print("  Da_rel=1/3 : da_q3_k1   vs da_q1_k0p333")
    print("  Da_rel=1/10: da_q10_k1  vs da_q1_k0p1")
    print("  Da_rel=1/30: da_q30_k1  vs da_q1_k0p033")


if __name__ == "__main__":
    main()
