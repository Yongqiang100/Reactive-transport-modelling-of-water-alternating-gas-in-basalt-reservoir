#!/usr/bin/env python3
"""
generate_boundary_tests.py -- NEW simulations testing the top-boundary / plume-geometry
critique (reviewers 2-4 and 3). All decks come from the CANONICAL builder at base rate, so
chemistry, near-well grid, injection and seeding are IDENTICAL to the baseline.

  open_top_<sc>.in : same 100 m domain, top boundary made OPEN (hydrostatic Dirichlet) so
                     free gas can migrate out of the top instead of ponding against a seal.

  tallcap_<sc>.in  : REDESIGNED taller-domain test. Domain extended to 200 m, but the added
                     upper 100 m is an INERT CAP -- same material (perm 1e-13, porosity 0.15)
                     so gas can rise into it, but all mineral surface areas set to ZERO so it
                     adds NO reactive volume. The outlet stays at z=0-100 m and the datum
                     stays at z=100 m, so through-flow and pressures in the reacting zone are
                     unchanged from the confined baseline. This isolates the buoyancy-headroom
                     effect (plume rising away from the reacting zone) from the volume/flow
                     confounds that contaminated the first tall test (dissolved control had
                     moved +18.6%). If dissolved is now ~unchanged, the test is clean.

Interpretation vs the confined baseline (01_baseline): dissolved (S1, no free gas) should be
~unchanged under both variants; if scCO2/WAG carbonation falls, the confining top was
flattering the two-phase configurations and the S1 > WAG > S2 ordering is conservative.

Usage:  python3 generate_boundary_tests.py    ->  ./decks/{open_top,tallcap}_<sc>.in
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "03_dape"))
from generate_dape_decks import build_deck  # canonical builder

SCEN = ["dissolved", "scco2", "wag6mo"]   # S1 control (no gas), S2 main, WAG intermediate

# grid for the 200 m domain (100 layers x 2 m in z); x-grading and everything else unchanged
TALL_NXYZ = "250 1 100"
TALL_DXYZ = "    100@1.d0 50@3.d0 50@9.d0 50@26.d0\n    1.d0\n    100@2.d0"

# ---------------- open-top variant (unchanged) ----------------
TOP_REGION = ("REGION top\n  FACE TOP\n  COORDINATES\n"
              "    0.d0 0.d0 100.d0\n    2000.d0 1.d0 100.d0\n  /\nEND\n\n")
TOP_BC = ("BOUNDARY_CONDITION top_open\n  FLOW_CONDITION outlet_bc\n"
          "  TRANSPORT_CONDITION outlet_tc\n  REGION top\nEND\n\n")
OUTLET_BC_BLOCK = ("BOUNDARY_CONDITION outlet_side\n  FLOW_CONDITION outlet_bc\n"
                   "  TRANSPORT_CONDITION outlet_tc\n  REGION outlet\nEND\n")

def make_open_top(deck):
    anchor = "REGION obs_10m\n  COORDINATE 10.d0 0.5d0 50.d0"
    assert deck.count(anchor) == 1, "obs_10m anchor"
    deck = deck.replace(anchor, TOP_REGION + anchor, 1)
    assert deck.count(OUTLET_BC_BLOCK) == 1, "outlet_side BC anchor"
    deck = deck.replace(OUTLET_BC_BLOCK, OUTLET_BC_BLOCK + "\n" + TOP_BC, 1)
    return deck

# ---------------- redesigned taller-domain variant (inert cap) ----------------
# concentrations copied verbatim from CONSTRAINT formation_water; minerals identical VF but
# ZERO surface area so the cap is chemically inert while keeping porosity = 0.15.
CAP_INERT = """

CONSTRAINT cap_inert
  CONCENTRATIONS
    H+        7.5d0   pH
    Ca++      2.0d-3  T
    Mg++      1.5d-3  T
    Fe++      5.0d-5  T
    Na+       1.0d-2  T
    K+        5.0d-4  T
    Al+++     1.0d-8  T
    SiO2(aq)  5.0d-4  T
    HCO3-     3.0d-3  T
    SO4--     1.0d-3  T
    Cl-       5.0d-3  T
    O2(aq)    1.0d-6  T
  /
  MINERALS
    Forsterite   0.05d0  0.d0 m^2/m^3
    Anorthite    0.30d0  0.d0 m^2/m^3
    Diopside     0.25d0  0.d0 m^2/m^3
    Kaolinite    0.25d0  0.d0 m^2/m^3
    Calcite      0.0d0   0.d0 m^2/m^3
    Magnesite    0.0d0   0.d0 m^2/m^3
    Siderite     0.0d0   0.d0 m^2/m^3
    Dolomite-ord 0.0d0   0.d0 m^2/m^3
    SiO2(am)     0.0d0   0.d0 m^2/m^3
  /
END"""

def make_tall_cap(deck):
    # (1) REGION all -> full 200 m; add REGION reactive (0-100) and REGION cap (100-200)
    old_all = ("REGION all\n  COORDINATES\n    0.d0 0.d0 0.d0\n"
               "    2000.d0 1.d0 100.d0\n  /\nEND")
    new_all = ("REGION all\n  COORDINATES\n    0.d0 0.d0 0.d0\n"
               "    2000.d0 1.d0 200.d0\n  /\nEND\n\n"
               "REGION reactive\n  COORDINATES\n    0.d0 0.d0 0.d0\n"
               "    2000.d0 1.d0 100.d0\n  /\nEND\n\n"
               "REGION cap\n  COORDINATES\n    0.d0 0.d0 100.d0\n"
               "    2000.d0 1.d0 200.d0\n  /\nEND")
    assert deck.count(old_all) == 1, "REGION all anchor"
    deck = deck.replace(old_all, new_all, 1)

    # (2) reacting zone keeps formation_water; cap gets the inert constraint
    old_ic = ("INITIAL_CONDITION\n  FLOW_CONDITION initial\n"
              "  TRANSPORT_CONDITION initial_tc\n  REGION all\nEND")
    new_ic = ("INITIAL_CONDITION\n  FLOW_CONDITION initial\n"
              "  TRANSPORT_CONDITION initial_tc\n  REGION reactive\nEND\n\n"
              "INITIAL_CONDITION cap_ic\n  FLOW_CONDITION initial\n"
              "  TRANSPORT_CONDITION cap_tc\n  REGION cap\nEND")
    assert deck.count(old_ic) == 1, "INITIAL_CONDITION anchor"
    deck = deck.replace(old_ic, new_ic, 1)

    # (3) cap transport condition -> inert constraint (after outlet_tc)
    anchor_tc = ("TRANSPORT_CONDITION outlet_tc\n  TYPE DIRICHLET_ZERO_GRADIENT\n"
                 "  CONSTRAINT_LIST\n    0.d0 formation_water\n  /\nEND")
    cap_tc = ("\n\nTRANSPORT_CONDITION cap_tc\n  TYPE DIRICHLET_ZERO_GRADIENT\n"
              "  CONSTRAINT_LIST\n    0.d0 cap_inert\n  /\nEND")
    assert deck.count(anchor_tc) == 1, "outlet_tc anchor"
    deck = deck.replace(anchor_tc, anchor_tc + cap_tc, 1)

    # (4) inert constraint definition (after formation_water block)
    anchor_fw = "    SiO2(am)     0.0d0   1.d2 m^2/m^3\n  /\nEND"
    assert deck.count(anchor_fw) == 1, "formation_water end anchor"
    deck = deck.replace(anchor_fw, anchor_fw + CAP_INERT, 1)
    return deck

def main():
    out = HERE / "decks"; out.mkdir(exist_ok=True)
    for stale in out.glob("tall_*.in"):   # superseded by tallcap_*.in
        stale.unlink(); print(f"  removed stale {stale.name}")
    for sc in SCEN:
        base100 = build_deck(sc, kappa=1.0, rate_mult=1.0)
        (out / f"open_top_{sc}.in").write_text(make_open_top(base100))
        base200 = build_deck(sc, kappa=1.0, rate_mult=1.0, nxyz=TALL_NXYZ, dxyz=TALL_DXYZ)
        (out / f"tallcap_{sc}.in").write_text(make_tall_cap(base200))
        print(f"  {sc}: open_top_{sc}.in  tallcap_{sc}.in")
    print(f"\ndecks -> {out}/  (confined baselines already in 01_baseline/)")

if __name__ == "__main__":
    main()
