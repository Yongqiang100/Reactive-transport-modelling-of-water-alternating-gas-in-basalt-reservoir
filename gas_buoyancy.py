#!/usr/bin/env python3
# Standalone: gas-saturation magnitudes + buoyancy-override ratio per baseline scenario.
import os, glob
from pathlib import Path
import numpy as np, h5py

ROOT = Path(os.environ.get("WAG_ROOT", "."))
BASE = [("base_dissolved", "Dissolved (S1)"), ("base_scco2", "scCO2 (S2)"),
        ("base_wag6mo", "WAG-6mo (S3)"), ("base_wag3mo", "WAG-3mo (S4)"),
        ("base_swag", "SWAG (S5)"), ("base_adaptive", "Adaptive (S6)")]

def find_h5(d):
    hs = [h for h in sorted(glob.glob(str(Path(d) / "*.h5"))) if not h.endswith("-restart.h5")]
    return hs[-1] if hs else None
def tgroups(f):
    return sorted([g for g in f.keys() if g.startswith("Time")],
                  key=lambda s: float(s.replace("Time:", "").strip().split()[0]))
def dset(group, base):
    for k in group.keys():
        if k.startswith(base):
            try:
                return np.array(group[k], dtype=float)
            except Exception:
                pass
    return None
def edges(f, letter, n):
    for g in ("Coordinates", "Domain", "Grid"):
        if g in f and isinstance(f[g], h5py.Group):
            for k in f[g].keys():
                if k.strip().lower().startswith(letter):
                    a = np.array(f[g][k], dtype=float).ravel()
                    if a.size == n + 1:
                        return a
    return None

print("=" * 80)
print("  GAS SATURATION & BUOYANCY OVERRIDE (study 01, ~30 yr end of injection)")
print("  Sg %, near-well x<100 m; upper z=50-80 m, lower z=20-50 m; up/lo >> 1 = buoyant plume")
print("=" * 80)
print("\n  %-15s%6s%9s%9s%8s%8s%7s" % ("scenario","t(yr)","meanSg%","peakSg%","upper%","lower%","up/lo"))
print("  " + "-" * 72)
for stem, label in BASE:
    d = ROOT / "01_baseline" / "runs" / stem
    h5 = find_h5(d)
    if not h5:
        print("  %-15s (no run)" % label); continue
    with h5py.File(h5, "r") as f:
        tg = tgroups(f)
        if not tg:
            print("  %-15s (no tsteps)" % label); continue
        times = [float(g.replace("Time:", "").strip().split()[0]) for g in tg]
        gi = min(range(len(tg)), key=lambda i: abs(times[i] - 30.0))
        g = f[tg[gi]]; t = times[gi]
        sg = dset(g, "Gas_Saturation")
        if sg is None:
            print("  %-15s%6.0f   no free-gas field (single-phase aqueous)" % (label, t)); continue
        s2 = sg[:, 0, :] if sg.ndim == 3 else sg
        nx, nz = s2.shape
        xe = edges(f, "x", nx); ze = edges(f, "z", nz)
        xc = 0.5 * (xe[:-1] + xe[1:]) if xe is not None else np.arange(nx, dtype=float)
        zc = 0.5 * (ze[:-1] + ze[1:]) if ze is not None else np.arange(nz, dtype=float)
        dom = float(s2.mean()); peak = float(s2.max())
        xnw = xc < 100.0
        up = (zc >= 50) & (zc <= 80); lo = (zc >= 20) & (zc < 50)
        um = float(s2[np.ix_(xnw, up)].mean()) if (xnw.any() and up.any()) else float("nan")
        lm = float(s2[np.ix_(xnw, lo)].mean()) if (xnw.any() and lo.any()) else float("nan")
        r = (um / lm) if (lm and lm > 1e-9) else float("inf")
        print("  %-15s%6.0f%9.3f%9.2f%8.3f%8.3f%7.1f" % (label, t, dom*100, peak*100, um*100, lm*100, r))
print("  " + "-" * 72)
