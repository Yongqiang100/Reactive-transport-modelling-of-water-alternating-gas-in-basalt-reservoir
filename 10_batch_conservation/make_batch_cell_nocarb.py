#!/usr/bin/env python3
"""Carbonate-OFF control: writes decks/batch_cell_nocarb.in from decks/batch_cell.in
with the four carbonate kinetic RATE_CONSTANTs driven to ~0 (-50 log10) so carbonate
cannot precipitate. Silicate kinetics (RATE_CONSTANT inside PREFACTOR blocks) untouched."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(HERE,"decks","batch_cell.in")
dst = os.path.join(HERE,"decks","batch_cell_nocarb.in")
CARB = {"Calcite","Magnesite","Siderite","Dolomite-ord"}
SIL  = {"Forsterite","Anorthite","Diopside","Kaolinite","SiO2(am)"}
lines = open(src).read().splitlines()
out=[]; in_mk=False; cur=None; patched=0
for ln in lines:
    s=ln.strip()
    if s=="MINERAL_KINETICS": in_mk=True
    if in_mk:
        if s in CARB: cur=s
        elif s in SIL: cur=s
        if cur in CARB and s.startswith("RATE_CONSTANT"):
            indent = ln[:len(ln)-len(ln.lstrip())]
            ln = f"{indent}RATE_CONSTANT -50.0000d0 mol/m^2-sec"
            patched+=1; cur=None
    out.append(ln)
txt = "\n".join(out)+"\n"
txt = txt.replace("SINGLE-CELL BATCH CARBON-CONSERVATION TEST",
                  "CARBONATE-OFF CONTROL (carbonate precipitation disabled)")
open(dst,"w").write(txt)
print(f"wrote {dst}  (disabled {patched} carbonate rate constants)")
