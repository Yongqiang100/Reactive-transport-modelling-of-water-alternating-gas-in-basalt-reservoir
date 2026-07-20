#!/usr/bin/env python3
"""Generate EVERY figure in a single run -- the one entry point for the figure set.

Reuses the readers and figure functions in make_manuscript_figures.py (mmf) and
analyse_transport_limitation.py (atl), and defines the three custom figures inline
(study-design schematic; the single-panel Da-Sigma regime = Figure 9 with the moved
labels). No separate standalone scripts are needed.

  figure                     source
  -------------------------  ---------------------------------------------
  fig_domain_schematic       mmf.fig_domain
  fig_kappa_crossover        atl.analyse_05 + atl.fig_crossover   (kinetic sensitivity)
  fig_gas_saturation_2d      mmf.fig_gas2d
  fig_baseline_comparison    mmf.fig_baseline
  fig_carbonate_breakdown    mmf.fig_carb_bar
  fig_spatial_profiles       mmf.fig_spatial
  fig_damkohler_sweep        mmf.fig_da_sweep
  fig_da_sigma_regime        fig_da_sigma_regime()   [single-panel, labelled]  (NOT mmf.fig_da_sigma)
  fig_study_design           fig_study_design()      [schematic, no data]

Run on Setonix from the folder holding the modules and the study dirs:
    cd $MYSCRATCH/WAG
    cd $MYSCRATCH/WAG && python3 generate_all_figures.py
Outputs -> ./figures/ (PDF + PNG). Set WAG_ROOT to point elsewhere if runs live in another tree.
"""
import os, sys, traceback
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

try:
    import make_manuscript_figures as mmf
    import analyse_transport_limitation as atl
except Exception as e:
    sys.exit("Could not import the figure modules from this directory:\n"
             f"  {e}\n"
             "Run from the folder containing make_manuscript_figures.py and "
             "analyse_transport_limitation.py, with the co2conv env active.")

OUT = mmf.OUT
OUT.mkdir(parents=True, exist_ok=True)

def _save(fig, name):
    fig.savefig(OUT / f"{name}.pdf"); fig.savefig(OUT / f"{name}.png", dpi=200); plt.close(fig)


# ======================================================================
#  Custom figure 1 -- study-design schematic (no data needed)
# ======================================================================
def fig_study_design():
    ORANGE = ("#c9781f", "#fbf0e2", "#c9781f"); BLUE = ("#2c6fbb", "#eaf1fb", "#2c6fbb")
    GREEN = ("#2e8b57", "#e8f5ee", "#2e8b57"); GREY_F, GREY_E, INK = "#eef1f4", "#8794a2", "#333333"
    fig, ax = plt.subplots(figsize=(8.6, 5.6)); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    CX = [0.29, 0.565, 0.84]; W = 0.245; cols = [ORANGE, BLUE, GREEN]
    def rbox(cx, y0, y1, fc, ec, w=W, lw=1.2, ls="-"):
        ax.add_patch(FancyBboxPatch((cx - w/2, y0), w, y1 - y0, boxstyle="round,pad=0.004,rounding_size=0.012",
                     linewidth=lw, facecolor=fc, edgecolor=ec, linestyle=ls, mutation_aspect=1.5))
    def vconn(cx, ytop, ybot, color, lw=1.3):
        ax.add_patch(FancyArrowPatch((cx, ytop), (cx, ybot), arrowstyle="-|>", mutation_scale=11,
                     lw=lw, color=color, shrinkA=0, shrinkB=0))
    B0, B1 = 0.860, 0.955; C0, C1 = 0.500, 0.830; F0, F1 = 0.280, 0.420; V0, V1 = 0.050, 0.190
    for y, lab in [((B0+B1)/2, "Control\nvaried"), ((C0+C1)/2, "Parameter\nswept"),
                   ((F0+F1)/2, "Key\nresult"), ((V0+V1)/2, "Numerical\ncontrols")]:
        ax.text(0.075, y, lab, ha="center", va="center", fontsize=8.2, color="#555", weight="bold", linespacing=1.15)
    for cx, t, (ban, _, _) in zip(CX, ["Injection rate", "Reaction kinetics", "Phase partitioning"], cols):
        rbox(cx, B0, B1, ban, ban)
        ax.text(cx, (B0+B1)/2, t, ha="center", va="center", color="white", fontsize=10.5, weight="bold")
    CARD = ["six schemes $\\times$ five rates\n0.3\u201330$\\times$ base injection\n($\\Rightarrow$ varies Damk\u00f6hler number)",
            "$\\kappa = 10^{-5}\\!-\\!10^{2}$  (7 decades)\ndissolved and scCO$_2$, base rate\n(+ repeated at 30$\\times$ rate)",
            "six fluid-delivery schemes\nCO$_2$ mole fraction 0.04\u20130.99\nwell position: top/mid/bottom"]
    NRUN = ["30 simulations", "56 simulations", "12 simulations"]
    for cx, c, n, (ban, _, ec) in zip(CX, CARD, NRUN, cols):
        rbox(cx, C0, C1, "white", ec)
        ax.text(cx, 0.700, c, ha="center", va="center", color=INK, fontsize=7.7, linespacing=1.4)
        ax.text(cx, 0.545, n, ha="center", va="center", color=ec, fontsize=7.8, weight="bold", style="italic")
        vconn(cx, B0 - 0.002, C1 + 0.002, ec)
    FIND = ["Sets the flow regime\n(supply vs. interface)",
            "No effect on carbonation\n(transport-limited)",
            "Sets total carbonation\n(buoyancy override)"]
    for cx, f, (ban, tint, ec) in zip(CX, FIND, cols):
        vconn(cx, C0 - 0.002, F1 + 0.002, ec); rbox(cx, F0, F1, tint, ec)
        ax.text(cx, (F0+F1)/2, f, ha="center", va="center", color=ec, fontsize=8.4, weight="bold", linespacing=1.3)
    bx0, bx1 = CX[0] - W/2, CX[-1] + W/2
    ax.add_patch(FancyBboxPatch((bx0, V0), bx1 - bx0, V1 - V0, boxstyle="round,pad=0.004,rounding_size=0.012",
                 linewidth=1.1, facecolor=GREY_F, edgecolor=GREY_E, mutation_aspect=3.0))
    ax.text((bx0+bx1)/2, (V0+V1)/2,
            "P\u00e9clet independence: molecular diffusivity $\\times10^{3}$ at fixed rate  $\\Rightarrow$  advection-dominated ($Pe\\gg1$)\n"
            "Grid convergence: five mesh resolutions  $\\Rightarrow$  results numerically converged",
            ha="center", va="center", fontsize=7.4, color=INK, linespacing=1.5)
    for cx, (ban, tint, ec) in zip(CX, cols):
        ax.add_patch(FancyArrowPatch((cx, F0 - 0.002), (cx, V1 + 0.002), arrowstyle="-|>", mutation_scale=8,
                     lw=0.9, color=GREY_E, linestyle=(0, (3, 2))))
    _save(fig, "fig_study_design")


# ======================================================================
#  Custom figure 2 -- single-panel Da-Sigma regime (Figure 9; moved labels)
# ======================================================================
def fig_da_sigma_regime():
    R_, T_, Vm = 8.314, 333.15, 3.69e-5
    sig3 = 0.7 * 2900 * 9.81 * 600; p0 = 1000 * 9.81 * 600; KOM = 1e5
    curves = mmf._rate_curve(lambda r: r["series"]["carb_mean"][-1])
    fig, ax = plt.subplots(figsize=(4.5, 4.2)); fig.subplots_adjust(bottom=0.22)
    n = 0
    for sc, (mus, vs, lbl, c, mk) in curves.items():
        dav, sigv = [], []
        for mu, tc in zip(mus, vs):
            om = max(1.01, 1 + tc * KOM); Pc = (R_ * T_ / Vm) * np.log(om)
            se = sig3 - p0 - 3e6 * min(mu, 5); sigv.append(Pc / max(se, 1e5)); dav.append(10.0 / mu)
        ax.plot(dav, sigv, color=c, ls='-', lw=0.6, alpha=0.4, zorder=3)
        ax.scatter(dav, sigv, c=c, s=28, marker=mk, edgecolors='white', linewidths=0.4, label=lbl, zorder=5)
        n += len(dav)
    ax.axvline(3.0, color='#e2e8f0', ls='-', lw=0.5, zorder=1); ax.axhline(800, color='#e2e8f0', ls='-', lw=0.5, zorder=1)
    bbox_kw = dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85); kw = dict(fontsize=5.5, fontstyle='italic', zorder=6)
    ax.text(0.03, 0.97, 'Transport-limited:\nnear-well self-sealing\n(conceptual; fractured)', transform=ax.transAxes, va='top', ha='left', color='#94a3b8', bbox=bbox_kw, **kw)
    ax.text(0.03, 0.06, 'Minimal alteration\n(all modes)', transform=ax.transAxes, va='bottom', ha='left', color='#64748b', bbox=bbox_kw, **kw)
    ax.text(0.97, 0.97, 'Supply-limited \u2192\nself-sealing / clogging\n(water-rich)', transform=ax.transAxes, va='top', ha='right', color='#dc2626', bbox=bbox_kw, **kw)
    ax.text(0.97, 0.06, 'Phase-partitioning-limited:\nsaturates (supercritical)', transform=ax.transAxes, ha='right', va='bottom', color='#16a34a', bbox=bbox_kw, **kw)
    ax1f, ay1f, ax2f, ay2f = 0.37, 0.90, 0.63, 0.90
    ax.annotate('', xy=(ax2f, ay2f), xytext=(ax1f, ay1f), xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='#94a3b8', lw=0.7))
    (PX1, PY1) = ax.transAxes.transform((ax1f, ay1f)); (PX2, PY2) = ax.transAxes.transform((ax2f, ay2f))
    ang = np.degrees(np.arctan2(PY2 - PY1, PX2 - PX1))
    ax.text((ax1f + ax2f) / 2, (ay1f + ay2f) / 2, 'Increasing injection rate', transform=ax.transAxes,
            rotation=ang, rotation_mode='anchor', ha='center', va='bottom', fontsize=5.5, color='#94a3b8')
    ax.set_xscale('log'); ax.invert_xaxis(); ax.set_xlabel('Damk\u00f6hler number (Da, decreasing)'); ax.set_ylabel('Normalised stress ratio (\u03a3)')
    if n:
        ax.legend(fontsize=5.5, ncol=6, loc='lower center', bbox_to_anchor=(0.5, -0.22), frameon=False,
                  columnspacing=0.8, handletextpad=0.3, markerscale=1.0)
    _save(fig, "fig_da_sigma_regime")


def _kappa_crossover():
    r5 = atl.analyse_05(); atl.fig_crossover(r5)


# ======================================================================
FIGS = [
    (mmf.fig_domain,       "fig_domain_schematic"),
    (_kappa_crossover,     "fig_kappa_crossover"),
    (mmf.fig_gas2d,        "fig_gas_saturation_2d"),
    (mmf.fig_baseline,     "fig_baseline_comparison"),
    (mmf.fig_carb_bar,     "fig_carbonate_breakdown"),
    (mmf.fig_spatial,      "fig_spatial_profiles"),
    (mmf.fig_da_sweep,     "fig_damkohler_sweep"),
    (fig_da_sigma_regime,  "fig_da_sigma_regime"),
    (fig_study_design,     "fig_study_design"),
]

def main():
    print("=" * 72); print("  Generating ALL figures ->", str(OUT)); print("=" * 72)
    failed = []
    for fn, name in FIGS:
        try:
            fn(); print(f"  ok    {name}")
        except Exception as e:
            failed.append((name, e)); print(f"  FAIL  {name}: {e}"); traceback.print_exc()
    print("\n=== figure check (figures/*.pdf) ===")
    missing = [name for _, name in FIGS if not (OUT / f"{name}.pdf").exists()]
    for _, name in FIGS:
        print(f"  [{'OK ' if (OUT / f'{name}.pdf').exists() else 'MISSING'}] figures/{name}.pdf")
    print("\n" + "=" * 72)
    print(f"  {len(FIGS) - len(missing)}/{len(FIGS)} figures present in {OUT}/")
    if failed:
        print("  failures:")
        for n, e in failed: print(f"    - {n}: {e}")
    print("=" * 72)

if __name__ == "__main__":
    main()
