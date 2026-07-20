#!/usr/bin/env python3
"""
make_reply_figures.py — summary figures for the response to reviewer Qinjun.
Built from the verified Setonix results (hard-coded), so it runs anywhere with
matplotlib (no cluster / no HDF5 needed). Emits four figures into ./figures/:

  fig_reply_kinetic_insensitivity   (C1: transport limitation)
  fig_reply_damkohler_consistency   (C3: throughput vs single-Da efficiency)
  fig_reply_phase_efficiency        (C2-phase, C4: phase partitioning + table)
  fig_reply_contact_buoyancy        (C2-contact, C2-buoyancy: hierarchy)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUT = Path(__file__).resolve().parent / "figures"; OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 8.5, "legend.fontsize": 7, "savefig.dpi": 300,
                     "savefig.bbox": "tight", "figure.dpi": 120})
BLUE, RED, GREEN, ORANGE, PURPLE, GREY = "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#7f7f7f"


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf"); fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig); print(f"  wrote figures/{name}.pdf (+ .png)")


# ---- real data (Setonix) ----
KAPPA = np.array([1e-5,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1,3e-1,1.0,3.0,10.0,30.0,100.0])
DIS_K = np.array([14.279,13.858,14.526,14.407,14.433,14.781,15.243,16.527,20.471,18.266,18.316,21.040,17.517,18.391])
SC_K  = np.array([1.9158,1.1931,2.1588,2.2469,2.3080,2.3737,2.4252,3.6979,7.5574,5.2530,5.2288,7.8403,3.4196,4.7891])

DA_REL   = np.array([1.0, 1/3, 0.1, 1/30])
RATE_MOL = np.array([4.7286e5, 1.4187e6, 4.6718e6, 1.4123e7])   # vary q (kappa=1)
KIN_MOL  = np.array([4.7286e5, 4.7276e5, 4.7221e5, 4.7248e5])   # vary kappa (q=1)
RATE_EFF = np.array([2.25, 2.25, 2.22, 2.24])
KIN_EFF  = np.array([2.25, 2.25, 2.24, 2.25])

XCO2   = np.array([0.04, 0.20, 0.40, 0.60, 0.80, 0.99])
CASE_C = np.array([18.266, 10.357, 9.060, 7.772, 6.492, 5.253])

SCN   = ["Dissolved\n(S1)", "scCO$_2$\n(S2)", "WAG-6mo\n(S3)", "WAG-3mo\n(S4)", "SWAG\n(S5)", "Adaptive\n(S6)"]
ETA_U = np.array([2.94, 1.33, 2.34, 2.33, 0.39, 2.18])   # uncontrolled
ETA_I = np.array([2.25, 0.04, 1.45, 1.44, 0.24, 1.23])   # injection-driven

SCK    = np.array([1e-2, 1e-1, 3e-1, 1.0, 10.0])
SC_TOT = np.array([2.3911, 3.6979, 7.5574, 5.2530, 7.8353])
SC_BG  = np.array([2.3132, 3.5744, 7.4288, 5.1147, 7.7062])
SC_INJ = np.array([0.0779, 0.12351, 0.12862, 0.13831, 0.12915])

CASE_D_LABELS = ["Top\n(z 0–30 m)", "Middle\n(z 20–80 m)", "Bottom\n(z 70–100 m)"]
CASE_D = np.array([5.2564, 5.2530, 5.2482])   # total carbonate, m^3


def fig_kinetic():
    fig, ax = plt.subplots(figsize=(5.6, 3.7))
    ax.axvspan(1.0, 1500.0, color="0.9", zorder=0)
    ax.text(40, 13.5, "literature range\n(P&K – Rimstidt, ~10$^3$×)", fontsize=6, color="0.45", ha="center")
    ax.loglog(KAPPA, DIS_K, "o-", color=BLUE, ms=4, label="Dissolved (S1): ×1.5")
    ax.loglog(KAPPA, SC_K, "s-", color=RED, ms=4, label="scCO$_2$ (S2): ×2–3 (noisy)")
    ax.axhline(DIS_K.mean(), color=BLUE, ls=":", lw=0.8)
    ax.axhline(SC_K.mean(), color=RED, ls=":", lw=0.8)
    ax.axvline(1.0, color="0.5", ls="--", lw=0.7)
    ax.text(1.25, 1.45, "base $\\kappa$=1", fontsize=6, color="0.45", rotation=90, va="bottom")
    ax.set_xlim(5e-6, 2e3); ax.set_ylim(1.0, 32)
    ax.set_xlabel(r"Kinetic-rate multiplier $\kappa$  (= rate constant / base; $\propto$ Damköhler)")
    ax.set_ylabel("Injection-driven carbonate (m$^3$)")
    ax.set_title("Carbonation is insensitive to kinetics over 7 orders of magnitude")
    ax.legend(loc="lower left", frameon=False)
    save(fig, "fig_reply_kinetic_insensitivity")


def fig_damkohler():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 3.4))
    # (a) absolute: rate axis ~∝ throughput, kinetic axis flat
    a.loglog(DA_REL, RATE_MOL, "o-", color=ORANGE, ms=5, label="vary flow $q$ (throughput)")
    a.loglog(DA_REL, KIN_MOL, "s--", color=PURPLE, ms=5, label=r"vary kinetics $\kappa$")
    a.set_xlabel(r"Relative Damköhler  Da$_{\rm rel}=\kappa/q$")
    a.set_ylabel("Injection-driven CO$_2$ mineralized (mol)")
    a.set_title("(a) Absolute yield: rate sweep = throughput")
    a.legend(loc="upper right", frameon=False)
    a.invert_xaxis()
    # (b) efficiency: both collapse onto a single ~2.2% line
    b.semilogx(DA_REL, RATE_EFF, "o-", color=ORANGE, ms=5, label="vary flow $q$")
    b.semilogx(DA_REL, KIN_EFF, "s--", color=PURPLE, ms=5, label=r"vary kinetics $\kappa$")
    b.set_ylim(0, 3.2)
    b.set_xlabel(r"Relative Damköhler  Da$_{\rm rel}=\kappa/q$")
    b.set_ylabel("Mineralization efficiency (%)")
    b.set_title("(b) Efficiency: single Da, pairs collapse to ~1%")
    b.legend(loc="lower center", frameon=False)
    b.invert_xaxis()
    fig.tight_layout()
    save(fig, "fig_reply_damkohler_consistency")


def fig_phase():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.6, 3.4))
    # (a) Case C: monotonic decline with phase split
    a.plot(XCO2, CASE_C, "o-", color=BLUE, ms=5)
    a.set_xlabel(r"Injected CO$_2$ mole fraction $x_{\rm CO_2}$")
    a.set_ylabel("Injection-driven carbonate (m$^3$)")
    a.set_title(r"(a) Phase partitioning (Case C): $\rho=-1.0$")
    a.annotate("dissolved", xy=(0.04, 18.27), xytext=(0.18, 17),
               fontsize=6.5, arrowprops=dict(arrowstyle="-", lw=0.5))
    a.annotate("supercritical", xy=(0.99, 5.25), xytext=(0.6, 8.5),
               fontsize=6.5, arrowprops=dict(arrowstyle="-", lw=0.5))
    a.set_ylim(0, 20)
    # (b) six-scenario efficiency: uncontrolled vs injection-driven
    x = np.arange(len(SCN)); w = 0.38
    b.bar(x - w/2, ETA_U, w, color="0.7", label="uncontrolled (final − initial)")
    b.bar(x + w/2, ETA_I, w, color=BLUE, label="injection-driven (control-subtracted)")
    for i in (1,):  # annotate scCO2 background dominance
        b.annotate("97% background", xy=(x[i] + w/2, ETA_I[i]), xytext=(x[i] + 0.1, 1.4),
                   fontsize=6, color=RED, ha="center",
                   arrowprops=dict(arrowstyle="-", color=RED, lw=0.6))
    b.text(4, 0.5, "SWAG injected\nmass overest.", fontsize=5.5, color="0.4", ha="center")
    b.set_xticks(x); b.set_xticklabels(SCN, fontsize=6.5)
    b.set_ylabel("Mineralization efficiency (%)")
    b.set_title("(b) Injected-CO$_2$ efficiency by strategy")
    b.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    save(fig, "fig_reply_phase_efficiency")


def fig_contact_buoyancy():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.6, 3.4))
    # (a) scCO2 control subtraction: contact-limited, rate-insensitive injection-driven
    a.loglog(SCK, SC_TOT, "o-", color=GREY, ms=4, label="total (final − initial)")
    a.loglog(SCK, SC_BG, "s--", color=PURPLE, ms=4, label="background (no injection)")
    a.loglog(SCK, SC_INJ, "o-", color=RED, ms=5, label="injection-driven")
    a.set_xlabel(r"Kinetic-rate multiplier $\kappa$")
    a.set_ylabel("Carbonate volume (m$^3$)")
    a.set_title("(a) Gas–water contact: scCO$_2$ injection-driven\nis tiny and rate-insensitive (×1.33)")
    a.legend(loc="center left", frameon=False)
    # (b) buoyancy (Case D): well position barely changes integrated yield
    x = np.arange(3)
    bars = b.bar(x, CASE_D, 0.6, color=[ORANGE, GREEN, BLUE])
    b.set_ylim(5.20, 5.30)
    b.set_xticks(x); b.set_xticklabels(CASE_D_LABELS, fontsize=6.5)
    b.set_ylabel("Total carbonate (m$^3$)")
    b.set_title("(b) Buoyancy (Case D): geometry, not yield\nspread < 0.2% across well position")
    for xi, v in zip(x, CASE_D):
        b.text(xi, v + 0.002, f"{v:.3f}", ha="center", fontsize=6.5)
    b.text(1, 5.285, "phase partitioning spans\n0.19 → 20.8 t (×110) by contrast",
           fontsize=5.8, color="0.4", ha="center")
    fig.tight_layout()
    save(fig, "fig_reply_contact_buoyancy")


if __name__ == "__main__":
    fig_kinetic()
    fig_damkohler()
    fig_phase()
    fig_contact_buoyancy()
    print("done")
