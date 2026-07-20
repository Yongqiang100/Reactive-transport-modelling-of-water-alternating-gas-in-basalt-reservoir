#!/usr/bin/env python3
"""
cation_balance.py -- where do the Mg2+ and Ca2+ in the carbonates come from?

For each baseline scenario, integrate the change in every mineral's volume fraction between
the first and last HDF5 output (sum of VF * cell_volume, converted to moles via molar volume)
and build the divalent-cation budget:

  Mg source  : Forsterite (Mg2SiO4, 2 Mg) + Diopside (CaMgSi2O6, 1 Mg)     -- dissolving
  Mg sink    : Magnesite (MgCO3) + Dolomite (CaMg(CO3)2)                    -- precipitating
  Ca source  : Anorthite (CaAl2Si2O8, 1 Ca) + Diopside (1 Ca)              -- dissolving
  Ca sink    : Calcite (CaCO3) + Dolomite                                   -- precipitating
  Fe source  : none (no Fe silicate in the assemblage) -> siderite starved

Conservation check: cations released by silicate dissolution must equal those fixed in
carbonate PLUS the change in the aqueous Mg2+/Ca2+ inventory (the remainder left in solution).
Closure near 100% confirms the silicate reactants are the cation source (and validates the
molar-volume conversion). Molar volumes are standard literature values; small deviations from
100% reflect differences vs the database's internal molar volumes.

Run on Setonix:  python3 cation_balance.py
"""
import sys
import numpy as np
from pathlib import Path
import h5py

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_manuscript_figures as mmf

SCEN = ["dissolved", "scco2", "wag6mo", "wag3mo", "swag", "adaptive"]
STUDY = sys.argv[1] if len(sys.argv) > 1 else "01_baseline"   # e.g. supcrt_case/01_baseline

# molar volumes (m^3/mol), standard values (Robie & Hemingway 1995; SUPCRT92)
VM = {"Forsterite": 43.79e-6, "Anorthite": 100.79e-6, "Diopside": 66.09e-6,
      "Magnesite": 28.02e-6, "Calcite": 36.93e-6, "Siderite": 29.38e-6,
      "Dolomite-ord": 64.37e-6}
# (Mg, Ca) per formula unit
CAT = {"Forsterite": (2, 0), "Diopside": (1, 1), "Anorthite": (0, 1),
       "Magnesite": (1, 0), "Dolomite-ord": (1, 1), "Calcite": (0, 1), "Siderite": (0, 0)}

def dmoles(g0, g1, cv, m):
    a0, a1 = mmf.get_dset(g0, f"{m}_VF"), mmf.get_dset(g1, f"{m}_VF")
    if a0 is None or a1 is None:
        return None
    return float(((a1 - a0) * cv).sum()) / VM[m]        # + precipitated, - dissolved

def aq_moles(g, cv, ion):
    """aqueous component inventory (mol) = conc[mol/L] * water_volume[L]."""
    c = None
    for base in (f"Total_{ion}", f"Total {ion}", ion):
        c = mmf.get_dset(g, base)
        if c is not None:
            break
    poro = mmf.get_dset(g, "Porosity")
    if c is None or poro is None:
        return None
    gas = mmf.get_dset(g, "Gas_Saturation")
    sl = (1.0 - gas) if gas is not None else 1.0
    return float((c * poro * sl * cv * 1.0e3).sum())     # L water = poro*sat*vol*1000

def daq(g0, g1, cv, ion):
    a0, a1 = aq_moles(g0, cv, ion), aq_moles(g1, cv, ion)
    return (a1 - a0) if (a0 is not None and a1 is not None) else float("nan")

def main():
    print(f"# case: {STUDY}")
    for sc in SCEN:
        d = mmf.rdir(STUDY, f"base_{sc}")
        h5 = mmf.find_h5(d)
        if not h5:
            print(f"\n=== {sc} ===  no run.h5 at {d}"); continue
        with h5py.File(h5, "r") as f:
            tg = mmf._tg(f); g0, g1 = f[tg[0]], f[tg[-1]]
            cv = mmf.cell_vol(f, mmf.carb_total(g1).shape)
            dn = {m: dmoles(g0, g1, cv, m) for m in VM}
            print(f"\n=== {sc}  (t0={mmf._tyr(tg[0]):.0f} yr -> t1={mmf._tyr(tg[-1]):.0f} yr) ===")
            print("  d(mineral) mol [-=dissolved,+=precipitated]: "
                  + "  ".join(f"{m} {dn[m]:+.2e}" for m in VM if dn[m] is not None))
            # classify each mineral by the SIGN of its change: dissolving -> source, precipitating -> sink
            for cat, idx, ion in (("Mg2+", 0, "Mg++"), ("Ca2+", 1, "Ca++")):
                rel = sum(-dn[m]*CAT[m][idx] for m in VM if dn[m] and dn[m] < 0 and CAT[m][idx])
                upt = sum( dn[m]*CAT[m][idx] for m in VM if dn[m] and dn[m] > 0 and CAT[m][idx])
                aq  = daq(g0, g1, cv, ion)
                src = ", ".join(f"{m} {-dn[m]*CAT[m][idx]:.2e}" for m in VM if dn[m] and dn[m] < 0 and CAT[m][idx]) or "none"
                snk = ", ".join(f"{m} {dn[m]*CAT[m][idx]:.2e}"  for m in VM if dn[m] and dn[m] > 0 and CAT[m][idx]) or "none"
                aqs = f"{aq:+.2e}" if aq == aq else "n/a"
                clo = f"{100*(upt+aq)/rel:.0f}%" if (rel and aq == aq) else "n/a"
                print(f"  {cat}: source [{src}]")
                print(f"        sink   [{snk}]   + aqueous {aqs}")
                print(f"        released {rel:.3e}  =  uptake {upt:.3e} + aqueous {aqs}   ->  closure {clo}")
    print("\n(Fe2+: no Fe-silicate in the assemblage -> siderite limited to the initial pore-water Fe.)")

if __name__ == "__main__":
    main()
