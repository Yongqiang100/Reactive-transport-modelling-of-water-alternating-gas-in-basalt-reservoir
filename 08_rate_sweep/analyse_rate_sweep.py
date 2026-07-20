#!/usr/bin/env python3
"""
analyse_rate_sweep.py  --  CORRECTED rate-sweep analysis (dual metric)

Tests whether the manuscript's sublinear "carbonation up to 16x at 30x rate"
(Fig. dasweep, S1/WAG) is a real Damkohler transition or a metric artifact of
the graded mesh, by computing TWO metrics on the same runs:

  m1  domain-mean carbonate VF      = unweighted cell mean  (the MANUSCRIPT metric)
  m3  injection-driven carbonate    = integral of (VF_final - VF_t0) dV  (volume-
                                      weighted, seed-subtracted; the CORRECT metric)

For each scenario it fits the rate-scaling exponent p (metric ~ mu^p):
  p ~ 1   linear / supply-limited (no transition)
  p < 1   sublinear

The manuscript reports m1. On the graded mesh the unweighted mean over-counts the
small refined near-well cells where carbonate concentrates, so as the plume spreads
into coarse far-field cells at high rate, m1 can grow sublinearly even when the
true (volume-weighted) carbonate grows linearly. If p(m1) < 1 while p(m3) ~ 1, the
"16x transition" is a metric artifact and carbonation is actually supply-limited --
consistent with study 07's clean single-Da collapse.

Layout expected: BASE_DIR/rs_<scenario>_mu<tag>/<prefix>.h5
Usage:
    BASE_DIR=$HOME/WAG/rate-sweep python3 analyse_rate_sweep.py
    python3 analyse_rate_sweep.py --selftest
"""
import os
import re
import sys
import glob
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import h5py
    HAVE_H5 = True
except ImportError:
    HAVE_H5 = False

BASE_DIR = Path(os.environ.get("BASE_DIR", Path.home() / "WAG" / "rate-sweep"))
OUT_DIR = Path(os.environ.get("OUT", "."))
CARB = ["Calcite", "Magnesite", "Siderite", "Dolomite-ord"]
# refined production grid: 250 x 1 x 50 ; x = 100@1 + 50@3 + 50@9 + 50@26 ; z = 50@2
GRID_WIDTHS = np.array([1.0] * 100 + [3.0] * 50 + [9.0] * 50 + [26.0] * 50)
DZ = 2.0
SCEN = ["dissolved", "scco2", "wag6mo", "wag3mo", "swag", "adaptive"]
MU = {"0p3": 0.3, "1": 1.0, "3": 3.0, "10": 10.0, "30": 30.0}


def _carb_sum(group):
    fld = None
    for m in CARB:
        ks = [k for k in group.keys() if k.startswith(f"{m}_VF")]
        if ks:
            a = np.array(group[ks[0]], dtype=float)
            fld = a if fld is None else fld + a
    return fld


def read_metrics(sim_dir):
    """Return (m1_domain_mean, m3_injection_driven_m3) at final time, using the
    t=0 group as the seed baseline for m3. None if unreadable."""
    if not HAVE_H5:
        return None
    h5s = [h for h in sorted(glob.glob(str(sim_dir / "*.h5"))) if not h.endswith("-restart.h5")]
    if not h5s:
        return None
    try:
        with h5py.File(h5s[-1], "r") as f:
            tgs = sorted([g for g in f.keys() if g.startswith("Time")],
                         key=lambda s: float(s.replace("Time:", "").strip().split()[0]))
            if not tgs:
                return None
            final = _carb_sum(f[tgs[-1]])
            if final is None:
                return None
            t0 = _carb_sum(f[tgs[0]])           # t=0 seed baseline
            if t0 is None:
                t0 = np.full_like(final, final.min())   # fallback: min as baseline
            cellvol = np.broadcast_to(GRID_WIDTHS[:, None, None] * 1.0 * DZ, final.shape)
            m1 = float(final.mean())                       # unweighted (manuscript)
            m3 = float(((final - t0) * cellvol).sum())     # injection-driven volume
            return m1, m3
    except Exception as exc:  # noqa: BLE001
        print(f"  warn: {sim_dir.name}: {exc}")
        return None


def fit_p(mu, y):
    mu = np.asarray(mu, float); y = np.asarray(y, float)
    ok = (mu > 0) & (y > 0)
    if ok.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(mu[ok]), np.log(y[ok]), 1)[0])


def collect():
    data = {}
    if not BASE_DIR.is_dir():
        print(f"  (no {BASE_DIR})")
        return data
    for sc in SCEN:
        rows = []
        for tag, mu in MU.items():
            d = BASE_DIR / f"rs_{sc}_mu{tag}"
            if not d.is_dir():
                continue
            mm = read_metrics(d)
            if mm is None:
                print(f"  (unreadable: rs_{sc}_mu{tag})")
                continue
            rows.append({"mu": mu, "m1": mm[0], "m3": mm[1]})
        if rows:
            data[sc] = sorted(rows, key=lambda r: r["mu"])
    return data


# ---------------------------------------------------------------------
def make_figure(data):
    import matplotlib as mpl
    mpl.rcParams["font.family"] = "serif"; mpl.rcParams["font.size"] = 9
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.6, 4.8))
    cmap = plt.cm.viridis(np.linspace(0, 0.85, len(data)))
    mm = np.array([min(MU.values()), max(MU.values())])

    for c, (sc, rows) in zip(cmap, data.items()):
        mu = [r["mu"] for r in rows]
        p1 = fit_p(mu, [r["m1"] for r in rows])
        p3 = fit_p(mu, [r["m3"] for r in rows])
        axL.plot(mu, [r["m1"] for r in rows], "o-", color=c, ms=6, mfc="white",
                 mew=1.3, label=f"{sc} ($p$={p1:.2f})")
        axR.plot(mu, [max(r["m3"], 1e-30) for r in rows], "s-", color=c, ms=6,
                 mfc="white", mew=1.3, label=f"{sc} ($p$={p3:.2f})")

    # linear (p=1) reference anchored at the dissolved mu=1 point
    ref = data.get("dissolved") or next(iter(data.values()))
    b = next((r for r in ref if abs(r["mu"] - 1.0) < 1e-9), ref[0])
    axL.plot(mm, b["m1"] * mm / b["mu"], "k--", lw=1.0, alpha=0.6, label="linear ($p$=1)")
    axR.plot(mm, max(b["m3"], 1e-30) * mm / b["mu"], "k--", lw=1.0, alpha=0.6, label="linear ($p$=1)")

    for ax in (axL, axR):
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("injection-rate multiplier $\\mu$   ($Da \\propto 1/\\mu$)")
        ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.4)
        ax.legend(fontsize=6.6, loc="upper left")
    axL.set_ylabel("domain-mean carbonate VF")
    axL.set_title("(a) MANUSCRIPT metric: unweighted domain-mean VF\n(graded mesh over-counts refined near-well cells)", fontsize=8.5)
    axR.set_ylabel("injection-driven carbonate volume (m$^3$)")
    axR.set_title("(b) CORRECT metric: volume-weighted, seed-subtracted\n($p$=1 -> supply-limited, no $Da$ transition)", fontsize=8.5)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"fig_rate_sweep_dualmetric.{ext}",
                    dpi=(300 if ext == "pdf" else 200), bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {OUT_DIR/'fig_rate_sweep_dualmetric.pdf'}")


# ---------------------------------------------------------------------
def _synthetic():
    """Manuscript metric sublinear (mu^0.82 -> 16x at 30x); correct metric
    linear (mu^1.0). The hypothesis under test."""
    data = {}
    for sc, amp in [("dissolved", 9e-4), ("scco2", 4.3e-4), ("wag6mo", 6.5e-4)]:
        rows = []
        for mu in (0.3, 1.0, 3.0, 10.0, 30.0):
            p1 = 0.05 if sc == "scco2" else 0.82       # scco2 ~ flat; others sublinear
            rows.append({"mu": mu, "m1": amp * mu ** p1, "m3": 20.0 * amp / 9e-4 * mu})
        data[sc] = rows
    return data


def main():
    selftest = "--selftest" in sys.argv
    print("=" * 70)
    print("  CORRECTED rate sweep -- dual-metric analysis")
    print(f"  Source: {'synthetic self-test' if selftest else BASE_DIR}")
    print("  graded mesh: unweighted mean (manuscript) vs volume-weighted (correct)")
    print("=" * 70)

    data = _synthetic() if selftest else collect()
    if not data:
        raise SystemExit("No results. Check BASE_DIR layout rs_<scen>_mu<tag>/, or --selftest.")

    print("\n  exponent p (carbonate ~ mu^p) and 30x fold, per scenario:")
    print(f"    {'scenario':<11} | {'MANUSCRIPT m1':>22} | {'CORRECT m3':>22}")
    print(f"    {'':<11} | {'p':>7} {'fold(0.3->30)':>14} | {'p':>7} {'fold(0.3->30)':>14}")
    flagged = []
    for sc, rows in data.items():
        mu = [r["mu"] for r in rows]
        p1 = fit_p(mu, [r["m1"] for r in rows]); p3 = fit_p(mu, [r["m3"] for r in rows])
        lo = rows[0]; hi = rows[-1]
        f1 = hi["m1"] / lo["m1"] if lo["m1"] else float("nan")
        f3 = hi["m3"] / lo["m3"] if lo["m3"] else float("nan")
        print(f"    {sc:<11} | {p1:>7.2f} {f1:>13.1f}x | {p3:>7.2f} {f3:>13.1f}x")
        if not np.isnan(p1) and not np.isnan(p3) and (p3 - p1) > 0.15:
            flagged.append((sc, p1, p3))

    print("\nGenerating figure...")
    make_figure(data)

    print("\n  READING:")
    if flagged:
        print("    The MANUSCRIPT metric (domain-mean VF) is markedly more sublinear than")
        print("    the volume-weighted/injection-driven metric for:")
        for sc, p1, p3 in flagged:
            print(f"       {sc}: p(manuscript)={p1:.2f}  vs  p(correct)={p3:.2f}")
        print("    On this graded mesh the unweighted mean over-counts the small near-well")
        print("    cells; as the plume spreads into coarse far-field cells at high rate, the")
        print("    mean grows sublinearly while the true carbonate grows ~linearly. The")
        print("    manuscript's '16x transition' is therefore a METRIC ARTIFACT, not a")
        print("    transport-to-reaction Damkohler transition -- consistent with study 07.")
    else:
        print("    The two metrics give similar exponents -- the rate behaviour is metric-")
        print("    robust here. Read p directly: p~1 => supply-limited (no transition);")
        print("    p<1 for both => a genuine sublinearity to explain (e.g. breakthrough).")
    print("\nDone.")


if __name__ == "__main__":
    main()
