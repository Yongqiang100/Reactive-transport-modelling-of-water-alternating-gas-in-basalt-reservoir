#!/usr/bin/env python3
# Verify the high-kappa carbonation drop at mu=30: is it PHYSICAL self-sealing
# (near-well porosity collapse) or a NUMERICAL artifact (timestep cuts / early stall)?
# Reads 11_kappa_mu30/runs/kappa30_<sc>_<ktag>/ {run.h5, run.log}. Standalone (h5py+numpy).
# Run from $MYSCRATCH/WAG:   python3 selfseal_check.py
import os, glob, re
from pathlib import Path
import numpy as np, h5py

ROOT = Path(os.environ.get("WAG_ROOT", "."))
STUDY = "11_kappa_mu30"
MIN_ALL = ["Calcite", "Magnesite", "Siderite", "Dolomite-ord",
           "Forsterite", "Anorthite", "Diopside", "Kaolinite", "SiO2(am)"]
PHI0 = 0.15
NEARWELL_X = 50.0
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
def parse_ktag(tag):
    m = re.match(r"k(\d+)e(m?)(\d+)$", tag)
    if not m: return None
    return float(m.group(1)) * 10.0 ** ((-1 if m.group(2) == "m" else 1) * int(m.group(3)))

def log_health(d):
    logs = sorted(glob.glob(str(Path(d) / "*.log")))
    if not logs: return None
    txt = open(logs[0], errors="ignore").read()
    steps = len(re.findall(r"(?im)^\s*step\s+\d+", txt))
    cuts = len(re.findall(r"(?i)cut\s*time\s*step|reducing\s*time\s*step|->\s*cut|time step cut|timestep cut", txt))
    stalls = len(re.findall(r"(?i)not\s*converge|did not converge|\bstall|newton.*fail|maximum number of.*iterations", txt))
    return {"steps": steps, "cuts": cuts, "stalls": stalls, "wall": ("Wall Clock Time" in txt)}

def field_stats(d):
    h5 = find_h5(d)
    if not h5: return None
    try: f = h5py.File(h5, "r")
    except Exception: return None
    with f:
        T = tg(f)
        if not T: return None
        g0, gf = f[T[0]], f[T[-1]]
        tfin = float(T[-1].replace("Time:", "").strip().split()[0])
        probe = dset(gf, "Magnesite_VF")
        if probe is None: return {"tfin": tfin, "err": "no mineral VF"}
        nx, ny, nz = (list(probe.shape) + [1, 1, 1])[:3]
        xe = edges(f, "x", nx); ze = edges(f, "z", nz)
        xc = 0.5 * (xe[:-1] + xe[1:]) if xe is not None else np.cumsum(GRID_WIDTHS) - GRID_WIDTHS / 2
        if xe is not None and ze is not None:
            cv = np.diff(xe)[:, None, None] * np.diff(ze)[None, None, :]
        else:
            cv = (GRID_WIDTHS[:, None, None] * DZ) * np.ones(probe.shape)
        dvf = np.zeros(probe.shape); nan = False
        for m in MIN_ALL:
            a1 = dset(gf, f"{m}_VF"); a0 = dset(g0, f"{m}_VF")
            if a1 is None: continue
            if a0 is None: a0 = np.zeros_like(a1)
            if not (np.isfinite(a1).all() and np.isfinite(a0).all()): nan = True
            dvf = dvf + (a1 - a0)
        dphi = -dvf; phi = PHI0 + dphi
        nw = xc < NEARWELL_X
        dphi_dom = float((dphi * cv).sum() / cv.sum())
        dphi_nw = float((dphi[nw] * cv[nw]).sum() / cv[nw].sum())
        return {"tfin": tfin, "dphi_dom": dphi_dom, "dphi_nw": dphi_nw,
                "phi_nw_min": float(phi[nw].min()), "phi_min": float(phi.min()), "nan": nan}

rd = ROOT / STUDY / "runs"
print("=" * 104)
print(f"  HIGH-kappa SELF-SEALING / CONVERGENCE CHECK  ({STUDY}, mu=30)")
print("  dphi<0 = porosity loss (precipitation). phi_nw_min -> ~0 = near-well clogging.")
print("  PHYSICAL self-sealing: dphi_nw plunges & phi_nw_min small WITH tfin=100, few cuts, no nan, phi_min>=0.")
print("  NUMERICAL artifact:    cuts/stalls spike, tfin<100, nan flagged, or phi_min<0.")
print("=" * 104)
for sc in ("dissolved", "scco2"):
    print(f"\n  --- {sc} (mu=30) ---")
    print(f"    {'kappa':>9}{'tfin':>7}{'dphi_dom(e4)':>14}{'dphi_nw(e4)':>13}{'phi_nw_min':>12}{'phi_min':>10}{'steps':>8}{'cuts':>7}{'stall':>7}  flags")
    subs = sorted(rd.glob(f"kappa30_{sc}_*"), key=lambda p: parse_ktag(p.name.split('_')[-1]) or 0)
    if not subs:
        print("    (no runs found)"); continue
    for p in subs:
        k = parse_ktag(p.name.split('_')[-1])
        fs = field_stats(p); lg = log_health(p)
        if fs is None: print(f"    {k:>9.0e}  (no h5)"); continue
        if "err" in fs: print(f"    {k:>9.0e}{fs['tfin']:>7.0f}   {fs['err']}"); continue
        st = "" if lg is None else f"{lg['steps']:>8}{lg['cuts']:>7}{lg['stalls']:>7}"
        flags = []
        if fs["tfin"] < 99.0: flags.append("STALLED")
        if fs.get("nan"): flags.append("NaN")
        if fs["phi_min"] < 0: flags.append("phi<0")
        if lg and lg["cuts"] > 50: flags.append("many-cuts")
        print(f"    {k:>9.0e}{fs['tfin']:>7.0f}{fs['dphi_dom']*1e4:>14.2f}{fs['dphi_nw']*1e4:>13.2f}"
              f"{fs['phi_nw_min']:>12.4f}{fs['phi_min']:>10.4f}{st}  {' '.join(flags)}")
print("\n" + "=" * 104)
