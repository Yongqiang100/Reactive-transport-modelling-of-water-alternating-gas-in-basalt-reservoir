#!/usr/bin/env python3
"""
diagnose_grid.py  —  localize the grid-dependence of carbonation

Re-reads the EXISTING study-06 runs (no new simulations) and splits the
carbonate volume fraction by region to test whether the non-convergence of
the domain mean is a near-well (injection-face) singularity:

  - injection region : x <= 20 m, 20 <= z <= 80 m  (the source cells)
  - reservoir        : 20 < x <= 500 m              (the swept reactive zone)
  - peak cell        : single hottest cell in the domain

If the injection-region mean and the peak cell grow with refinement while the
reservoir mean converges, the domain mean is singularity-dominated and the
reservoir mean (or front position) is the grid-robust metric to report.

Usage (run where the outputs live):
    BASE_DIR=$HOME/WAG/grid-resolution python3 diagnose_grid.py
    # or just: python3 diagnose_grid.py   (uses the same default as the analysis)
"""

import os
import re
import glob
from pathlib import Path
import numpy as np

try:
    import h5py
except ImportError:
    raise SystemExit("h5py required (conda activate geochem)")

BASE_DIR = Path(os.environ.get("BASE_DIR", Path.home() / "WAG" / "grid-resolution"))
CARB = ["Calcite", "Magnesite", "Siderite", "Dolomite-ord"]

GRID_WIDTHS = {
    "2m":   [2.0] * 50 + [10.0] * 40 + [30.0] * 50,
    "1m":   [1.0] * 100 + [10.0] * 40 + [30.0] * 50,
    "0p5m": [0.5] * 200 + [10.0] * 40 + [30.0] * 50,
}
CELL = {"2m": 2.0, "1m": 1.0, "0p5m": 0.5}
DMULT = {"1x": 1e-9, "10x": 1e-8, "100x": 1e-7}
Z_CENTERS = np.arange(25) * 4.0 + 2.0          # 25 cells of 4 m -> 2,6,...,98


def xcen(g):
    w = np.array(GRID_WIDTHS[g]); e = np.concatenate([[0], np.cumsum(w)])
    return 0.5 * (e[:-1] + e[1:])


def carb_field(sim_dir):
    h5s = sorted(glob.glob(str(sim_dir / "*.h5")))
    if not h5s:
        return None
    with h5py.File(h5s[-1], "r") as f:
        tg = sorted([g for g in f.keys() if g.startswith("Time")],
                    key=lambda s: float(s.replace("Time:", "").strip().split()[0]))
        if not tg:
            return None
        g = f[tg[-1]]
        fld = None
        for m in CARB:
            ks = [k for k in g.keys() if k.startswith(f"{m}_VF")]
            if ks:
                a = np.array(g[ks[0]], float)
                fld = a if fld is None else fld + a
        return fld          # (NX, NY, NZ)


def regions(fld, g):
    x = xcen(g)
    w = np.array(GRID_WIDTHS[g])
    f2 = fld[:, 0, :] if fld.ndim == 3 else fld          # (NX, NZ)
    inj = (x[:, None] <= 20.0) & (Z_CENTERS[None, :] >= 20.0) & (Z_CENTERS[None, :] <= 80.0)
    res_x = (x > 20.0) & (x <= 500.0)                    # x-only mask -> all z rows
    cellvol = w[:, None] * 1.0 * 4.0                     # cell volume (dy=1, dz=4 m)
    return {
        "domain_mean": float(fld.mean()),
        "inj_mean":    float(f2[inj].mean()) if inj.any() else float("nan"),
        "reservoir_mean": float(f2[res_x, :].mean()) if res_x.any() else float("nan"),
        # grid-robust integral: m^3 of carbonate precipitated in the reservoir
        "reservoir_vol": float((f2 * cellvol)[res_x, :].sum()) if res_x.any() else float("nan"),
        "peak_cell":   float(fld.max()),
    }


def main():
    print("=" * 74)
    print(f"  Carbonate localization diagnostic   ({BASE_DIR})")
    print("=" * 74)
    rows = {}
    for sub in sorted(BASE_DIR.iterdir()):
        m = re.match(r"grid_(\w+?)_D(\w+)$", sub.name)
        if not (sub.is_dir() and m):
            continue
        g, d = m.group(1), m.group(2)
        if g not in GRID_WIDTHS:
            continue
        fld = carb_field(sub)
        if fld is None:
            print(f"  (no HDF5 in {sub.name})"); continue
        rows[(g, d)] = regions(fld, g)

    # Convergence ladder at D=1e-9
    print("\n  Grid convergence at D=1e-9 (does the growth live near the well?)")
    print(f"  {'cell':>6} {'domain_mean':>13} {'inj_region':>13} {'reservoir':>13} "
          f"{'resvr_vol_m3':>13} {'peak_cell':>13}")
    for g in ("2m", "1m", "0p5m"):
        r = rows.get((g, "1x"))
        if not r:
            continue
        print(f"  {CELL[g]:>5}m {r['domain_mean']:>13.4e} {r['inj_mean']:>13.4e} "
              f"{r['reservoir_mean']:>13.4e} {r['reservoir_vol']:>13.4e} {r['peak_cell']:>13.4e}")
    # Quick verdict
    seq = [(g, rows[(g, "1x")]) for g in ("2m", "1m", "0p5m") if (g, "1x") in rows]
    if len(seq) >= 2:
        def growth(key):
            v = [r[key] for _, r in seq]
            return (v[-1] - v[0]) / abs(v[0]) * 100 if v[0] else float("nan")
        print("\n  Change 2m -> 0.5m:")
        for key, lbl in [("domain_mean", "domain mean"), ("inj_mean", "injection region"),
                          ("reservoir_mean", "reservoir mean"),
                          ("reservoir_vol", "reservoir volume m3"), ("peak_cell", "peak cell")]:
            print(f"    {lbl:22s}: {growth(key):+6.0f}%")
        print("\n  Interpretation: if 'injection region' and 'peak cell' grow strongly while")
        print("  'reservoir volume' stays roughly flat, the domain mean is singularity-")
        print("  dominated and the RESERVOIR VOLUME (m3 carbonate, injection cells excluded)")
        print("  is the grid-robust magnitude to report.")

    # Also confirm Pe-invariance is reservoir-wide (not just domain mean)
    print("\n  Pe sub-sweep at 0.5 m (reservoir mean should be flat across D):")
    for d in ("1x", "10x", "100x"):
        r = rows.get(("0p5m", d))
        if r:
            print(f"    D={DMULT[d]:.0e}: reservoir {r['reservoir_mean']:.4e}  "
                  f"domain {r['domain_mean']:.4e}")


if __name__ == "__main__":
    main()
