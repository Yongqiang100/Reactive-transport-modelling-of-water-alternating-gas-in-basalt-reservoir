#!/usr/bin/env python3
# v2: proper vertical structure of the free-gas plume (full-domain), to quantify buoyancy override.
import os, glob
from pathlib import Path
import numpy as np, h5py
ROOT = Path(os.environ.get("WAG_ROOT", "."))
BASE = [("base_dissolved","Dissolved (S1)"),("base_scco2","scCO2 (S2)"),
        ("base_wag6mo","WAG-6mo (S3)"),("base_wag3mo","WAG-3mo (S4)"),
        ("base_swag","SWAG (S5)"),("base_adaptive","Adaptive (S6)")]
def find_h5(d):
    hs=[h for h in sorted(glob.glob(str(Path(d)/"*.h5"))) if not h.endswith("-restart.h5")]
    return hs[-1] if hs else None
def tgroups(f):
    return sorted([g for g in f.keys() if g.startswith("Time")],
                  key=lambda s: float(s.replace("Time:","").strip().split()[0]))
def dset(group,base):
    for k in group.keys():
        if k.startswith(base):
            try: return np.array(group[k],dtype=float)
            except Exception: pass
    return None
def edges(f,letter,n):
    for g in ("Coordinates","Domain","Grid"):
        if g in f and isinstance(f[g],h5py.Group):
            for k in f[g].keys():
                if k.strip().lower().startswith(letter):
                    a=np.array(f[g][k],dtype=float).ravel()
                    if a.size==n+1: return a
    return None
print("="*84)
print("  GAS PLUME VERTICAL STRUCTURE (full domain, ~30 yr) -- buoyancy quantification")
print("  z_com = gas-saturation-weighted mean elevation; f(z>60) = fraction of gas above 60 m")
print("="*84)
for stem,label in BASE:
    d=ROOT/"01_baseline"/"runs"/stem; h5=find_h5(d)
    if not h5: print(f"\n  {label}: (no run)"); continue
    with h5py.File(h5,"r") as f:
        tg=tgroups(f)
        if not tg: print(f"\n  {label}: (no tsteps)"); continue
        times=[float(g.replace("Time:","").strip().split()[0]) for g in tg]
        gi=min(range(len(tg)),key=lambda i:abs(times[i]-30.0)); g=f[tg[gi]]; t=times[gi]
        sg=dset(g,"Gas_Saturation")
        if sg is None: print(f"\n  {label} (t={t:.0f}): no free-gas field (single-phase aqueous)"); continue
        s2=sg[:,0,:] if sg.ndim==3 else sg; nx,nz=s2.shape
        xe=edges(f,"x",nx); ze=edges(f,"z",nz)
        xc=0.5*(xe[:-1]+xe[1:]) if xe is not None else np.arange(nx,dtype=float)
        zc=0.5*(ze[:-1]+ze[1:]) if ze is not None else np.arange(nz,dtype=float)
        dx=np.diff(xe) if xe is not None else np.ones(nx); dz=np.diff(ze) if ze is not None else np.ones(nz)
        W=dx[:,None]*dz[None,:]; gw=s2*W; tot=float(gw.sum())
        print(f"\n  {label} (t={t:.0f} yr)   domain x=[0,{(xe[-1] if xe is not None else nx):.0f}] z=[0,{(ze[-1] if ze is not None else nz):.0f}] m")
        if tot<=0: print("    no free gas present"); continue
        z_com=float((gw*zc[None,:]).sum()/tot); x_com=float((gw*xc[:,None]).sum()/tot)
        f60=float(gw[:,zc>60].sum()/tot); f_inj=float(gw[:,(zc>=20)&(zc<=80)].sum()/tot)
        pk=np.unravel_index(int(np.argmax(s2)),s2.shape)
        print(f"    mean Sg {s2.mean()*100:6.2f}%   peak Sg {s2.max()*100:6.2f}% at (x={xc[pk[0]]:.0f} m, z={zc[pk[1]]:.0f} m)")
        print(f"    gas centre-of-mass: z_com={z_com:5.1f} m  x_com={x_com:6.0f} m   f(z>60 m)={f60*100:5.1f}%   f(20<z<80)={f_inj*100:5.1f}%")
        line="    vertical profile (mean Sg %% over all x): "
        for lo,hi in [(0,20),(20,40),(40,60),(60,80),(80,100)]:
            m=(zc>=lo)&((zc<hi) if hi<100 else (zc<=hi))
            line+=f"z{lo}-{hi}:{(s2[:,m].mean()*100 if m.any() else float('nan')):5.2f}  "
        print(line)
print("\n"+"="*84)
