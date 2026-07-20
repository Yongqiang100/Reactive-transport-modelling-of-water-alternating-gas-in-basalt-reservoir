#!/usr/bin/env python3
"""Study-design schematic for the manuscript: three controls (columns), read top-to-bottom
as control -> parameter swept (+ run count) -> key result, with a numerical/transport
controls band underpinning all three. Left-margin labels name each tier. No data needed.
Run:  python3 fig_study_design.py  ->  writes figures/fig_study_design.pdf (+ .png)"""
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.dpi": 300,
                     "savefig.bbox": "tight", "savefig.pad_inches": 0.06})
OUT = Path(os.environ.get("WAG_ROOT", Path(__file__).resolve().parent)) / "figures"
OUT.mkdir(exist_ok=True)

# (saturated banner, light tint, border)
ORANGE = ("#c9781f", "#fbf0e2", "#c9781f")
BLUE   = ("#2c6fbb", "#eaf1fb", "#2c6fbb")
GREEN  = ("#2e8b57", "#e8f5ee", "#2e8b57")
GREY_F, GREY_E, INK = "#eef1f4", "#8794a2", "#333333"

fig, ax = plt.subplots(figsize=(8.6, 5.6))
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

CX = [0.29, 0.565, 0.84]; W = 0.245
cols = [ORANGE, BLUE, GREEN]

def rbox(cx, y0, y1, fc, ec, w=W, lw=1.2, ls="-"):
    ax.add_patch(FancyBboxPatch((cx - w/2, y0), w, y1 - y0,
                 boxstyle="round,pad=0.004,rounding_size=0.012",
                 linewidth=lw, facecolor=fc, edgecolor=ec, linestyle=ls, mutation_aspect=1.5))

def vconn(cx, ytop, ybot, color, lw=1.3):
    ax.add_patch(FancyArrowPatch((cx, ytop), (cx, ybot), arrowstyle="-|>",
                 mutation_scale=11, lw=lw, color=color, shrinkA=0, shrinkB=0))

# tier y-bands
B0, B1 = 0.860, 0.955     # banner
C0, C1 = 0.500, 0.830     # card
F0, F1 = 0.280, 0.420     # finding
V0, V1 = 0.050, 0.190     # validation band

# ---- left-margin tier labels ----
for y, lab in [((B0+B1)/2, "Control\nvaried"), ((C0+C1)/2, "Parameter\nswept"),
               ((F0+F1)/2, "Key\nresult"), ((V0+V1)/2, "Numerical\ncontrols")]:
    ax.text(0.075, y, lab, ha="center", va="center", fontsize=8.2, color="#555",
            weight="bold", linespacing=1.15)

# ---- banners: the control varied ----
TITLE = ["Injection rate", "Reaction kinetics", "Phase partitioning"]
for cx, t, (ban, _, _) in zip(CX, TITLE, cols):
    rbox(cx, B0, B1, ban, ban)
    ax.text(cx, (B0+B1)/2, t, ha="center", va="center", color="white",
            fontsize=10.5, weight="bold")

# ---- cards: parameter swept (+ run count) ----
CARD = [
 "six schemes $\\times$ five rates\n0.3\u201330$\\times$ base injection\n($\\Rightarrow$ varies Damk\u00f6hler number)",
 "$\\kappa = 10^{-5}\\!-\\!10^{2}$  (7 decades)\ndissolved and scCO$_2$, base rate\n(+ repeated at 30$\\times$ rate)",
 "six fluid-delivery schemes\nCO$_2$ mole fraction 0.04\u20130.99\nwell position: top/mid/bottom",
]
NRUN = ["30 simulations", "56 simulations", "12 simulations"]
for cx, c, n, (ban, _, ec) in zip(CX, CARD, NRUN, cols):
    rbox(cx, C0, C1, "white", ec)
    ax.text(cx, 0.700, c, ha="center", va="center", color=INK, fontsize=7.7, linespacing=1.4)
    ax.text(cx, 0.545, n, ha="center", va="center", color=ec, fontsize=7.8,
            weight="bold", style="italic")
    vconn(cx, B0 - 0.002, C1 + 0.002, ec)          # banner -> card

# ---- finding strips: key result ----
FIND = ["Sets the flow regime\n(supply vs. interface)",
        "No effect on carbonation\n(transport-limited)",
        "Sets total carbonation\n(buoyancy override)"]
for cx, f, (ban, tint, ec) in zip(CX, FIND, cols):
    vconn(cx, C0 - 0.002, F1 + 0.002, ec)          # card -> finding
    rbox(cx, F0, F1, tint, ec)
    ax.text(cx, (F0+F1)/2, f, ha="center", va="center", color=ec,
            fontsize=8.4, weight="bold", linespacing=1.3)

# ---- validation band spanning all columns ----
bx0, bx1 = CX[0] - W/2, CX[-1] + W/2
ax.add_patch(FancyBboxPatch((bx0, V0), bx1 - bx0, V1 - V0,
             boxstyle="round,pad=0.004,rounding_size=0.012",
             linewidth=1.1, facecolor=GREY_F, edgecolor=GREY_E, mutation_aspect=3.0))
ax.text((bx0+bx1)/2, (V0+V1)/2,
        "P\u00e9clet independence: molecular diffusivity $\\times10^{3}$ at fixed rate  "
        "$\\Rightarrow$  advection-dominated ($Pe\\gg1$)\n"
        "Grid convergence: five mesh resolutions  $\\Rightarrow$  results numerically converged",
        ha="center", va="center", fontsize=7.4, color=INK, linespacing=1.5)
for cx, (ban, tint, ec) in zip(CX, cols):
    ax.add_patch(FancyArrowPatch((cx, F0 - 0.002), (cx, V1 + 0.002), arrowstyle="-|>",
                 mutation_scale=8, lw=0.9, color=GREY_E, linestyle=(0, (3, 2))))

fig.savefig(OUT / "fig_study_design.pdf")
fig.savefig(OUT / "fig_study_design.png", dpi=200)
print("done")
