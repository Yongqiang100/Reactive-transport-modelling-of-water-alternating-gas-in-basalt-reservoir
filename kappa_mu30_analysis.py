#!/usr/bin/env python3
# Standalone: kinetic-crossover carbonation(kappa) at BASE rate vs mu=30.
# Reads 11_kappa_mu30/runs (kappa30_*) and 05_kinetic_crossover/runs (kappa_*).
# Needs only h5py + numpy. Run from $MYSCRATCH/WAG:   python3 kappa_mu30_analysis.py
import os, glob, math, re
from pathlib import Path
import numpy as np, h5py

ROOT = Path(os.environ.get("WAG_ROOT", "."))
CARB = ["Calcite", "Magnesite", "Siderite", "Dolomite-ord"]
GRID_WIDTHS = np.array([1.0] * 100 + [3.0] * 50 + [9.0] * 50 + [26.0] * 50); DZ = 2.0

def find_h5(d):
    hs = [h for h in sorted(glob.glob(str(Path(d) / "*.h5"))) if not h.endswith("-restart.h5")]
    return hs[-1] if hs else None
def tg(f):
    return sorted([g for g in f.keys() if g.startswith("Time")],
                  key=lambda s: float(s.replace("Time:", "").strip().split()[0]))
def dset(g, base):
    for k in g.keys():
        if k.startswith(base):
            try: return np.array(g[k], dtype=float)
            except Exception: pass
    return None
def edges(f, letter, n):
    for gg in ("Coordinates", "Domain", "Grid"):
        if gg in f and isinstance(f[gg], h5py.Group):
            for k in f[gg].keys():
                if k.strip().lower().startswith(letter):
                    a = np.array(f[gg][k], dtype=float).ravel()
                    if a.size == n + 1: return a
    return None
def cell_vol(f, shape):
    nx, ny, nz = (list(shape) + [1, 1, 1])[:3]
    xe = edges(f, "x", nx); ze = edges(f, "z", nz); ye = edges(f, "y", ny)
    if xe is not None and ze is not None:
        dy = np.diff(ye) if ye is not None else np.array([1.0])
        return np.diff(xe)[:, None, None] * dy[None, :, None] * np.diff(ze)[None, None, :]
    if tuple(shape) == (len(GRID_WIDTHS), 1, 50):
        return GRID_WIDTHS[:, None, None] * 1.0 * DZ * np.ones(shape)
    return np.ones(shape)
def inj_m3(d):
    h5 = find_h5(d)
    if not h5: return None
    try: f = h5py.File(h5, "r")
    except Exception: return None
    with f:
        T = tg(f)
        if not T: return None
        def carb(gn):
            try: g = f[gn]
            except Exception: return None
            tot = None
            for m in CARB:
                a = dset(g, f"{m}_VF")
                if a is not None: tot = a if tot is None else tot + a
            return tot
        c1 = None
        for gn in reversed(T):
            c1 = carb(gn)
            if c1 is not None: break
        c0 = None
        for gn in T:
            c0 = carb(gn)
            if c0 is not None: break
        if c1 is None: return None
        if c0 is None: c0 = np.zeros_like(c1)
        cv = cell_vol(f, c1.shape)
        return float(((c1 - c0) * cv).sum())
def parse_ktag(tag):
    m = re.match(r"k(\d+)e(m?)(\d+)$", tag)
    if not m: return None
    mant = float(m.group(1)); sign = -1 if m.group(2) == "m" else 1; exp = int(m.group(3))
    return mant * 10.0 ** (sign * exp)
def collect(study, prefix):
    rd = ROOT / study / "runs"; series = {"dissolved": [], "scco2": []}
    if not rd.exists(): return series
    for sub in sorted(rd.glob(f"{prefix}_*")):
        m = re.match(rf"{prefix}_(dissolved|scco2)_(k\w+)$", sub.name)
        if not m: continue
        k = parse_ktag(m.group(2)); v = inj_m3(sub)
        if k is not None and v is not None:
            series[m.group(1)].append((k, v))
    for sc in series: series[sc].sort()
    return series
def tail_slope(ks, cs):
    ks = np.array(ks, float); cs = np.array(cs, float)
    lo = ks <= ks.min() * 30
    if lo.sum() >= 2 and (cs[lo] > 0).all():
        return float(np.polyfit(np.log10(ks[lo]), np.log10(cs[lo]), 1)[0])
    return float("nan")

base = collect("05_kinetic_crossover", "kappa")
hi = collect("11_kappa_mu30", "kappa30")
print("=" * 78)
print("  KINETIC-CROSSOVER TEST: carbonation(kappa) at BASE rate vs mu=30")
print("  inj_m3 = injection-driven carbonate volume (m^3).")
print("  Flat across kappa => transport-limited;  falling at low kappa (slope->1) => reaction-limited.")
print("=" * 78)
for sc in ("dissolved", "scco2"):
    print(f"\n  --- {sc} ---")
    print(f"    {'kappa':>10} {'inj_m3 (mu=1)':>16} {'inj_m3 (mu=30)':>16}")
    b = dict(base[sc]); h = dict(hi[sc])
    for k in sorted(set(list(b) + list(h))):
        bv = f"{b[k]:.4g}" if k in b else "-"
        hv = f"{h[k]:.4g}" if k in h else "-"
        print(f"    {k:>10.0e} {bv:>16} {hv:>16}")
    for lbl, s in (("mu=1 ", base[sc]), ("mu=30", hi[sc])):
        if len(s) >= 3:
            ks = [k for k, _ in s]; cs = [v for _, v in s]
            var = (max(cs) / min(cs)) if min(cs) > 0 else float("inf")
            sl = tail_slope(ks, cs)
            verdict = ("REACTION-LIMITED tail" if sl > 0.5 else
                       "transport-limited plateau" if sl < 0.2 else "partial/intermediate")
            print(f"    [{lbl}] variation {var:6.2f}x over {math.log10(max(ks)/min(ks)):.0f} decades; "
                  f"low-kappa slope {sl:+.2f} -> {verdict}")
print("\n" + "=" * 78)
