#!/usr/bin/env python3
import os, glob
from pathlib import Path
import numpy as np, h5py
ROOT = Path(os.environ.get("WAG_ROOT", "."))
CARB = ["Calcite","Magnesite","Siderite","Dolomite-ord"]; NEARWELL_X = 50.0
GRID_WIDTHS = np.array([1.0]*100 + [3.0]*50 + [9.0]*50 + [26.0]*50); DZ = 2.0
BASE = [("base_dissolved","Dissolved (S1)"),("base_scco2","scCO2 (S2)"),
        ("base_wag6mo","WAG-6mo (S3)"),("base_wag3mo","WAG-3mo (S4)"),
        ("base_swag","SWAG (S5)"),("base_adaptive","Adaptive (S6)")]
def find_h5(d):
    hs=[h for h in sorted(glob.glob(str(Path(d)/"*.h5"))) if not h.endswith("-restart.h5")]
    return hs[-1] if hs else None
def tg(f):
    return sorted([g for g in f.keys() if g.startswith("Time")], key=lambda s: float(s.replace("Time:","").strip().split()[0]))
def dset(g, base):
    for k in g.keys():
        if k.startswith(base):
            try: return np.array(g[k], dtype=float)
            except Exception: pass
    return None
def edges(f, letter, n):
    for gg in ("Coordinates","Domain","Grid"):
        if gg in f and isinstance(f[gg], h5py.Group):
            for k in f[gg].keys():
                if k.strip().lower().startswith(letter):
                    a=np.array(f[gg][k], dtype=float).ravel()
                    if a.size==n+1: return a
    return None
print("="*92)
print("  NEAR-WELL CARBONATE CONCENTRATION (baseline) -- tests WAG 'redistribution' claim")
print("  near-well = x < %g m.  f_nw = near-well carbonate / domain-total carbonate." % NEARWELL_X)
print("="*92)
print(f"\n  {'scenario':<16}{'t(yr)':>6}{'total carb (m3)':>18}{'near-well (m3)':>16}{'f_nw (%)':>10}")
print("  "+"-"*66)
res={}
for stem,label in BASE:
    d=ROOT/"01_baseline"/"runs"/stem; h5=find_h5(d)
    if not h5: print(f"  {label:<16} (no run)"); continue
    with h5py.File(h5,"r") as f:
        T=tg(f)
        if not T: print(f"  {label:<16} (no tsteps)"); continue
        t=float(T[-1].replace("Time:","").strip().split()[0]); g0,gf=f[T[0]],f[T[-1]]
        probe=dset(gf,"Magnesite_VF")
        if probe is None: print(f"  {label:<16} (no mineral VF)"); continue
        nx,ny,nz=(list(probe.shape)+[1,1,1])[:3]
        xe=edges(f,"x",nx); ze=edges(f,"z",nz)
        xc=0.5*(xe[:-1]+xe[1:]) if xe is not None else np.cumsum(GRID_WIDTHS)-GRID_WIDTHS/2
        cv=(np.diff(xe)[:,None,None]*np.diff(ze)[None,None,:]) if (xe is not None and ze is not None) else (GRID_WIDTHS[:,None,None]*DZ)*np.ones(probe.shape)
        dcarb=np.zeros(probe.shape)
        for m in CARB:
            a1=dset(gf,f"{m}_VF"); a0=dset(g0,f"{m}_VF")
            if a1 is None: continue
            if a0 is None: a0=np.zeros_like(a1)
            dcarb=dcarb+(a1-a0)
        tot=float((dcarb*cv).sum()); nw=xc<NEARWELL_X; nwc=float((dcarb[nw]*cv[nw]).sum())
        frac=100.0*nwc/tot if tot>0 else float('nan'); res[label]=frac
        print(f"  {label:<16}{t:>6.0f}{tot:>18.3f}{nwc:>16.3f}{frac:>10.1f}")
print("  "+"-"*66)
d1=res.get("Dissolved (S1)"); wags=[res[k] for k in ("WAG-6mo (S3)","WAG-3mo (S4)","SWAG (S5)") if k in res]
if d1 is not None and wags:
    wmean=sum(wags)/len(wags)
    print(f"\n  dissolved f_nw = {d1:.1f}%   WAG mean f_nw = {wmean:.1f}%")
    if wmean < d1-1: print("  -> WAG LESS near-well concentrated: redistribution claim SUPPORTED.")
    elif wmean > d1+1: print("  -> WAG MORE near-well concentrated: claim NOT supported (keep softened).")
    else: print("  -> comparable: claim NOT clearly supported.")
print("\n"+"="*92)
