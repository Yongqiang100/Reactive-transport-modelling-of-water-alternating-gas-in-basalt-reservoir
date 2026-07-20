#!/usr/bin/env python3
"""
compare_boundary.py -- quantify the boundary-test outcome against the confined baseline.

For each scenario it reports the TOTAL injection-driven carbonate VOLUME (m^3), integrated as
sum(VF * cell_volume) over the final 2-D field. Total volume (not domain-mean VF) is used so
the confined 100 m and the taller 200 m domains are directly comparable -- a domain mean would
be diluted by the extra inert rock in the tall case.

Three top-boundary treatments per scenario:
  confined : no-flow top at z=100 m         (01_baseline/runs/base_<sc>)
  open     : open Dirichlet top at z=100 m   (12_boundary_tests/runs/open_top_<sc>)
  tall     : inert-cap redesign, 200 m       (12_boundary_tests/runs/tallcap_<sc>)

Reading the result:
  * dissolved (S1, no free gas) is the CONTROL: it should be ~unchanged (<~5%) across
    boundaries. A large S1 shift -- especially under 'tall' -- would signal a grid/pressure
    artifact, not physics, and should be investigated before trusting the gas cases.
  * scCO2 (S2) and WAG: if carbonate DROPS under open/tall, the confined top was flattering
    the two-phase configurations => the reported S1 > WAG > S2 ordering is CONSERVATIVE.
  * Sanity check: the 'confined' column should reproduce the paper's ~16.3 m^3 (S1) and
    ~3.2 m^3 (S2). If it does not, stop -- the runs or paths are wrong.

Run on Setonix once the boundary runs finish, then paste the printed tables back:
    cd $MYSCRATCH/WAG && python3 12_boundary_tests/compare_boundary.py
"""
import sys
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import make_manuscript_figures as mmf  # ROOT resolves via its own file path

SCEN = ["dissolved", "scco2", "wag6mo"]
BND = [("confined", "01_baseline", "base_{}"),
       ("open",     "12_boundary_tests", "open_top_{}"),
       ("tall",     "12_boundary_tests", "tallcap_{}")]   # redesigned: inert cap, outlet at 0-100
MINS = ["Calcite", "Magnesite", "Siderite", "Dolomite-ord"]


def carb_volume_m3(r):
    """Total carbonate volume = sum(VF * cell_volume), from the final 2-D field.
    carb2d is [nx, nz]; cell widths from the grid edges; dy = 1 m (2-D slice)."""
    c = r["carb2d"]
    dx = np.diff(r["xe"]); dz = np.diff(r["ze"])
    return float((c * dx[:, None] * dz[None, :]).sum())


def load():
    data = {}
    for sc in SCEN:
        data[sc] = {}
        for lab, study, patt in BND:
            d = mmf.rdir(study, patt.format(sc))
            r = None
            try:
                r = mmf.read_run(d)
            except Exception as e:
                print(f"  [skip] {sc}/{lab}: {e}")
            if r is None:
                print(f"  [skip] {sc}/{lab}: no run.h5 at {d}")
                data[sc][lab] = None
                continue
            sg = r["series"].get("gas_mean")
            data[sc][lab] = dict(
                vol=carb_volume_m3(r),
                sg=float(sg[-1]) if sg is not None and len(sg) else float("nan"),
                ph=r["phases"],
            )
    return data


def main():
    d = load()

    print("\n=== TOTAL injection-driven carbonate volume (m^3), and change vs confined ===")
    print(f"{'scenario':<10}{'confined':>11}{'open':>11}{'tall':>11}   {'dopen%':>8}{'dtall%':>8}")
    for sc in SCEN:
        r = d[sc]
        v = {k: (r[k]["vol"] if r.get(k) else float("nan")) for k in ("confined", "open", "tall")}
        c0 = v["confined"]
        do = 100 * (v["open"] - c0) / c0 if c0 else float("nan")
        dt = 100 * (v["tall"] - c0) / c0 if c0 else float("nan")
        print(f"{sc:<10}{v['confined']:>11.3g}{v['open']:>11.3g}{v['tall']:>11.3g}   {do:>+8.1f}{dt:>+8.1f}")

    print("\n=== S1/S2 carbonate-volume ratio under each boundary (paper: ~5.1x) ===")
    for lab in ("confined", "open", "tall"):
        s1, s2 = d["dissolved"].get(lab), d["scco2"].get(lab)
        if s1 and s2 and s2["vol"]:
            print(f"  {lab:<9}: {s1['vol'] / s2['vol']:.2f}")

    print("\n=== domain-mean gas saturation, final (confirms plume escape/rise) ===")
    for sc in ("scco2", "wag6mo"):
        r = d[sc]
        print(f"  {sc:<9}: " + "  ".join(
            f"{lab} {(r[lab]['sg'] if r.get(lab) else float('nan')):.3f}" for lab in ("confined", "open", "tall")))

    print("\n=== carbonate mineral split (mean VF): is magnesite dominance robust? ===")
    for sc in SCEN:
        print(f"  {sc}:")
        for lab in ("confined", "open", "tall"):
            r = d[sc].get(lab)
            if not r:
                continue
            print(f"    {lab:<9}: " + "  ".join(f"{m} {r['ph'].get(m, 0.0):.2e}" for m in MINS))

    print("\n[interpret] S1 ~unchanged + S2/WAG lower under open/tall => ordering is conservative.")


if __name__ == "__main__":
    main()
