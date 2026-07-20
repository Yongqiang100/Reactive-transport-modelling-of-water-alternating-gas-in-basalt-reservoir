#!/usr/bin/env python3
"""Molar carbon budget for the single-cell batch conservation test."""
import sys, glob, os
import numpy as np, h5py

M_CO2 = 0.04401
MOLARV = {"Calcite":3.6934e-5,"Magnesite":2.8018e-5,"Siderite":2.9378e-5,"Dolomite-ord":6.4365e-5}
STO = {"Calcite":1,"Magnesite":1,"Siderite":1,"Dolomite-ord":2}
CARB = list(MOLARV)

def find(d, pat):
    g = glob.glob(os.path.join(d, pat)); return g[0] if g else None
def tg(f):
    ts = [(float(k.split(":")[1].split()[0]), k) for k in f.keys() if k.startswith("Time:")]
    return [k for _, k in sorted(ts)]
def cell_vol(f):
    c = f["Coordinates"]; x=c["X [m]"][:]; y=c["Y [m]"][:]; z=c["Z [m]"][:]
    return float(np.diff(x).sum()*np.diff(y).sum()*np.diff(z).sum())
def mineral_C(d):
    h = find(d, "*.h5")
    if not h: return None
    with h5py.File(h,"r") as f:
        T=tg(f); g0,gf=f[T[0]],f[T[-1]]; V=cell_vol(f)
        def tot(g):
            s=0.0
            for m in CARB:
                for k in g.keys():
                    if k.startswith(m+"_VF"): s+=float(np.asarray(g[k]).sum())*V/MOLARV[m]*STO[m]
            return s
        return tot(g0),tot(gf)
def mas_pool(d):
    p = find(d, "*-mas.dat")
    if not p: return None
    L=[ln for ln in open(p).read().splitlines() if ln.strip()]
    cols=[c.strip().strip('"').strip() for c in L[0].split(",")]
    rows=[]
    for ln in L[1:]:
        try: rows.append([float(x) for x in ln.split()])
        except ValueError: pass
    def col(row,keys,nk=()):
        for i,c in enumerate(cols):
            cl=c.lower()
            if all(k in cl for k in keys) and not any(n in cl for n in nk) and i<len(row): return row[i]
        return 0.0
    r0,rf=rows[0],rows[-1]
    def has(keys,nk=()):
        for c in cols:
            cl=c.lower()
            if all(k in cl for k in keys) and not any(n in cl for n in nk): return True
        return False
    def pair(keys,nk=(),scale=1.0): return (col(r0,keys,nk)*scale, col(rf,keys,nk)*scale)
    if has(["global","hco3"]): dic=pair(["global","hco3"]); dic_name="HCO3-"
    elif has(["global","co2(aq)"]): dic=pair(["global","co2(aq)"]); dic_name="CO2(aq)"
    else: dic=(0.0,0.0); dic_name="??"
    return {"dic":dic,"dic_name":dic_name,
            "agas":pair(["global","air","gas","kg"],scale=1.0/M_CO2),
            "aliq":pair(["global","air","liquid","kg"],scale=1.0/M_CO2),
            "cl":pair(["global","cl-"]),"ca":pair(["global","ca++"]),
            "mg":pair(["global","mg++"]),"fe":pair(["global","fe++"])}
SIL={"Forsterite":4.302e-5,"Anorthite":1.019e-4,"Diopside":6.369e-5}
def silicate_delta(d):
    h=find(d,"*.h5")
    if not h: return None
    with h5py.File(h,"r") as f:
        T=tg(f); g0,gf=f[T[0]],f[T[-1]]; V=cell_vol(f); out={}
        for m in SIL:
            for k in gf.keys():
                if k.startswith(m+"_VF"):
                    out[m]=float(np.asarray(gf[k]).sum())*V/SIL[m]-float(np.asarray(g0[k]).sum())*V/SIL[m]
        return out
def main():
    d=sys.argv[1] if len(sys.argv)>1 else "runs/batch_cell"
    mc=mineral_C(d); mp=mas_pool(d)
    if mc is None or mp is None:
        print(f"could not read run in {d} (need *.h5 and *-mas.dat)"); return
    print("="*64); print("  SINGLE-CELL BATCH — molar carbon budget (closed system)")
    print(f"  run: {d}"); print("="*64)
    rows=[(f"aqueous DIC ({mp['dic_name']})",mp["dic"]),("dissolved CO2 (Air liq)",mp["aliq"]),
          ("free gas CO2 (Air gas)",mp["agas"]),("carbonate mineral C",mc)]
    print(f"\n  {'reservoir':<26}{'t=0':>14}{'t=final':>14}{'Δ mol C':>13}")
    t0=tf=0.0
    for name,(a,b) in rows:
        t0+=a; tf+=b; print(f"  {name:<26}{a:>14.4f}{b:>14.4f}{b-a:>13.4f}")
    print(f"  {'TOTAL C':<26}{t0:>14.4f}{tf:>14.4f}{tf-t0:>13.4f}")
    dmin=mc[1]-mc[0]; dmob=(mp["dic"][1]-mp["dic"][0])+(mp["aliq"][1]-mp["aliq"][0])+(mp["agas"][1]-mp["agas"][0])
    print(f"\n  closed cell: injection = 0, boundary flux = 0  =>  d(total C) must be ~ 0")
    print(f"  d mineral C = {dmin:+.4f} mol   d mobile C = {dmob:+.4f} mol   d total = {tf-t0:+.4f} mol")
    print(f"\n  VERDICT:")
    if abs(tf-t0)<0.05*max(abs(t0),1e-9):
        print("   d total ~ 0  =>  carbon CONSERVED in one closed cell with this configuration.")
    else:
        print(f"   d total = {tf-t0:+.4f} mol != 0  =>  carbon NOT conserved (duplicated by carbonate step).")
    cl0,clf=mp["cl"]
    print("\n  " + "-"*58); print("  CLOSED-SYSTEM SANITY"); print("  " + "-"*58)
    print(f"  tracer Cl-   : {cl0:>12.4f} -> {clf:>12.4f}   d {clf-cl0:+.4e} mol  ({'closed' if abs(clf-cl0)<1e-6*max(abs(cl0),1e-9) else 'FLUX!'})")
    for nm,key in [("Ca++","ca"),("Mg++","mg"),("Fe++","fe")]:
        a,b=mp[key]; print(f"  cation {nm:<5}: {a:>12.4f} -> {b:>12.4f}   d {b-a:+.4e} mol")
    sil=silicate_delta(d)
    if sil:
        for m,dv in sil.items(): print(f"  silicate {m:<11} d = {dv:+.4f} mol")
if __name__=="__main__": main()
