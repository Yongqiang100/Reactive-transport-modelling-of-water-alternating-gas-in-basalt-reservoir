#!/usr/bin/env python3
"""
analyse_da_consistency.py  --  Damkohler collapse test

Tests whether carbonation is governed by a SINGLE, consistently-defined
Damkohler number, by varying Da two independent ways and checking collapse:

    Da_rel = kappa / q          (reaction rate / advection rate, normalized)
    rate axis  (kappa=1): Da_rel = 1/q      kinetic axis (q=1): Da_rel = kappa

Observable: carbonation efficiency, made comparable across runs as
    eff_rel = (injection-driven CO2 mineralized) / (CO2 supplied)
with supply proportional to the rate multiplier, and injection-driven
mineralization isolated by subtracting the kappa-matched no-injection control
cell-by-cell (final_injection - final_control), then normalized to baseline.

If a single Da governs, eff_rel is a function of Da_rel ALONE: the rate-axis
and kinetic-axis points fall on one curve, and each matched-Da pair agrees.
If they diverge, the injection-rate "transition" is a supply / residence-time
effect, NOT a Damkohler transition, and the manuscript's scaling is being used
inconsistently.

Usage:
    python3 analyse_da_consistency.py
    BASE_DIR=/path python3 analyse_da_consistency.py
    python3 analyse_da_consistency.py --selftest
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

BASE_DIR = Path(os.environ.get("BASE_DIR", Path.home() / "WAG" / "da-consistency"))
OUT_DIR = Path(os.environ.get("OUT", "."))

# study-07 decks use the refined production grid (NXYZ 250 1 50)
GRID_WIDTHS = [1.0] * 100 + [3.0] * 50 + [9.0] * 50 + [26.0] * 50
DZ = 2.0
MINERAL_PROPS = {                  # molar volume (m^3/mol), CO2 per formula
    "Calcite":      (3.6934e-5, 1),
    "Magnesite":    (2.8018e-5, 1),
    "Siderite":     (2.9378e-5, 1),
    "Dolomite-ord": (6.4365e-5, 2),
}
M_CO2 = 0.04401
RTAG = {"0": 0.0, "1": 1.0, "3": 3.0, "10": 10.0, "30": 30.0}
KTAG = {"1": 1.0, "0p333": 1.0 / 3, "0p1": 0.1, "0p033": 1.0 / 30}


def read_minerals(sim_dir):
    """{mineral: field(NX,1,NZ)} at the final time, or None."""
    if not HAVE_H5:
        return None
    h5s = sorted(glob.glob(str(sim_dir / "*.h5")))
    if not h5s:
        return None
    try:
        with h5py.File(h5s[-1], "r") as f:
            tgs = sorted([g for g in f.keys() if g.startswith("Time")],
                         key=lambda s: float(s.replace("Time:", "").strip().split()[0]))
            if not tgs:
                return None
            g = f[tgs[-1]]
            out = {}
            for m in MINERAL_PROPS:
                ks = [k for k in g.keys() if k.startswith(f"{m}_VF")]
                if ks:
                    out[m] = np.array(g[ks[0]], dtype=float)
            return out or None
    except Exception as exc:  # noqa: BLE001
        print(f"  warn: {sim_dir.name}: {exc}")
        return None


def injection_driven_co2(inj_dir, ctrl_dir):
    """CO2 mineralized by injection = integral of (VF_inj - VF_ctrl) dV per
    mineral -> moles CO2. Control removes the background re-equilibration."""
    mi = read_minerals(inj_dir)
    mc = read_minerals(ctrl_dir)
    if mi is None or mc is None:
        return None
    wx = np.array(GRID_WIDTHS)
    co2_mol = 0.0
    per_mineral = {}
    for m, (Vm, nco2) in MINERAL_PROPS.items():
        if m not in mi:
            continue
        fld_i = mi[m]
        fld_c = mc.get(m, 0.0)
        cellvol = np.broadcast_to(wx[:, None, None] * 1.0 * DZ, fld_i.shape)
        dvol = float(((fld_i - fld_c) * cellvol).sum())   # m^3 mineral (injection-driven)
        per_mineral[m] = dvol
        co2_mol += (dvol / Vm) * nco2
    return {"CO2_mol": co2_mol, "per_mineral": per_mineral}


def parse_run(name):
    m = re.match(r"da_q(\w+?)_k(\w+)$", name)
    if not m:
        return None
    rt, kt = m.group(1), m.group(2)
    if rt not in RTAG or kt not in KTAG:
        return None
    return RTAG[rt], KTAG[kt], kt


def collect():
    """Return {(q, kappa): {'mineral_mol','eff_rel','Da_rel','axis'}}."""
    runs = {}
    if not BASE_DIR.is_dir():
        print(f"  (no {BASE_DIR})")
        return runs
    # locate controls by kappa-tag
    ctrl = {}
    for sub in BASE_DIR.iterdir():
        pr = parse_run(sub.name) if sub.is_dir() else None
        if pr and pr[0] == 0.0:
            ctrl[pr[2]] = sub
    for sub in sorted(BASE_DIR.iterdir()):
        pr = parse_run(sub.name) if sub.is_dir() else None
        if not pr or pr[0] == 0.0:
            continue
        q, kap, kt = pr
        if kt not in ctrl:
            print(f"  (no control da_q0_k{kt} for {sub.name}; skipping)")
            continue
        idc = injection_driven_co2(sub, ctrl[kt])
        if idc is None:
            print(f"  (unreadable: {sub.name})")
            continue
        eff_rel = idc["CO2_mol"] / q                  # supply ~ q
        axis = "rate" if abs(kap - 1.0) < 1e-9 else "kappa"
        runs[(q, kap)] = {"mineral_mol": idc["CO2_mol"], "eff_rel": eff_rel,
                          "Da_rel": kap / q, "axis": axis,
                          "per_mineral": idc["per_mineral"]}
    return runs


# ---------------------------------------------------------------------
def make_figure(runs):
    import matplotlib as mpl
    mpl.rcParams["font.family"] = "serif"; mpl.rcParams["font.size"] = 9
    base = runs.get((1.0, 1.0))
    e0 = base["eff_rel"] if base else 1.0

    rate = sorted([(v["Da_rel"], v["eff_rel"] / e0) for k, v in runs.items()
                   if v["axis"] == "rate" or k == (1.0, 1.0)])
    kap = sorted([(v["Da_rel"], v["eff_rel"] / e0) for k, v in runs.items()
                  if v["axis"] == "kappa" or k == (1.0, 1.0)])

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    if rate:
        x, y = zip(*rate)
        ax.plot(x, y, "o-", color="#1f4e9c", ms=8, mfc="white", mew=1.6,
                label="Da varied via injection rate $q$  ($\\kappa=1$)")
    if kap:
        x, y = zip(*kap)
        ax.plot(x, y, "s--", color="#c0392b", ms=8, mfc="white", mew=1.6,
                label="Da varied via kinetics $\\kappa$  ($q=1$)")
    # matched-Da guide lines
    for dr in (1.0 / 3, 0.1, 1.0 / 30):
        ax.axvline(dr, color="0.8", lw=0.8, ls=":", zorder=0)
    ax.axhline(1.0, color="0.7", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("Damk\u00f6hler number  $Da_{\\mathrm{rel}} = \\kappa/q$  (baseline = 1)")
    ax.set_ylabel("carbonation efficiency (injection-driven), normalized")
    ax.set_title("Damk\u00f6hler collapse test: does a single $Da$ govern carbonation?",
                 fontsize=10)
    ax.grid(True, which="both", ls="-", lw=0.3, alpha=0.25)
    ax.legend(fontsize=8, loc="lower right")
    # annotate the interpretation directly (upper-left empty region)
    ax.annotate("coincide \u2192 single consistent $Da$\n"
                "diverge \u2192 rate axis is supply/residence",
                xy=(0.03, 0.18), xycoords="axes fraction", ha="left", va="bottom",
                fontsize=7.5, color="0.4")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"fig_da_consistency.{ext}",
                    dpi=(300 if ext == "pdf" else 200), bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {OUT_DIR/'fig_da_consistency.pdf'}")


# ---------------------------------------------------------------------
def _synthetic():
    """Encodes the expected outcome: kinetic axis flat (kappa never limits in
    range), rate axis sublinear (16x for 30x -> q^0.8) -> the two DIVERGE."""
    M0 = 3.55e5
    runs = {}
    for q in (1.0, 3.0, 10.0, 30.0):                 # rate axis (kappa=1)
        M = M0 * q ** 0.8
        runs[(q, 1.0)] = {"mineral_mol": M, "eff_rel": M / q,
                          "Da_rel": 1.0 / q, "axis": "rate", "per_mineral": {}}
    for kap in (1.0 / 30, 0.1, 1.0 / 3):             # kinetic axis (q=1)
        M = M0                                        # flat
        runs[(1.0, kap)] = {"mineral_mol": M, "eff_rel": M,
                            "Da_rel": kap, "axis": "kappa", "per_mineral": {}}
    return runs


def main():
    selftest = "--selftest" in sys.argv
    print("=" * 64)
    print("  Damkohler consistency (collapse) test")
    print(f"  Source: {'synthetic self-test' if selftest else BASE_DIR}")
    print("=" * 64)

    runs = _synthetic() if selftest else collect()
    if not runs:
        raise SystemExit("No results. Run the sweep first, or use --selftest.")
    base = runs.get((1.0, 1.0))
    if not base:
        raise SystemExit("Baseline run da_q1_k1 missing; cannot normalize.")
    e0 = base["eff_rel"]

    print("\n  run (q, kappa)      Da_rel    inj-CO2 [mol]   eff_rel (norm)")
    for (q, kap), v in sorted(runs.items(), key=lambda kv: (kv[1]["axis"], kv[1]["Da_rel"])):
        print(f"    q={q:<5g} k={kap:<7.4g}  {v['Da_rel']:>7.3f}  {v['mineral_mol']:>13.3e}  "
              f"{v['eff_rel']/e0:>8.3f}   [{v['axis']}]")

    # Matched-Da pair test: rate vs kinetics at identical Da_rel
    print("\n  MATCHED-Da PAIRS (efficiency via rate vs via kinetics at same Da):")
    print(f"    {'Da_rel':>8} {'eff_rate':>10} {'eff_kappa':>10} {'ratio':>8}  verdict")
    pairs = [(1.0 / 3, (3.0, 1.0), (1.0, 1.0 / 3)),
             (0.1,     (10.0, 1.0), (1.0, 0.1)),
             (1.0 / 30, (30.0, 1.0), (1.0, 1.0 / 30))]
    ratios = []
    for dr, rk, kk in pairs:
        vr, vk = runs.get(rk), runs.get(kk)
        if not (vr and vk):
            continue
        er, ek = vr["eff_rel"] / e0, vk["eff_rel"] / e0
        ratio = er / ek if ek else float("nan")
        ratios.append(ratio)
        verdict = "agree (Da governs)" if 0.8 <= ratio <= 1.25 else "DIVERGE (not Da)"
        print(f"    {dr:>8.3f} {er:>10.3f} {ek:>10.3f} {ratio:>8.2f}  {verdict}")

    print("\nGenerating figure...")
    make_figure(runs)

    if ratios:
        worst = min(ratios)
        print("\n  VERDICT:")
        if all(0.8 <= r <= 1.25 for r in ratios):
            print("    The rate and kinetic axes COLLAPSE onto one curve -> carbonation is")
            print("    governed by a single, consistent Damkohler number. The rate-axis")
            print("    transition is a genuine Da (transport-to-reaction) transition.")
        else:
            print("    The two axes DIVERGE at matched Da (efficiency ratio down to "
                  f"{worst:.2f}).")
            print("    Carbonation is NOT a function of Da = kappa/q alone: reducing Da by")
            print("    raising q lowers efficiency, while reducing Da by lowering kappa does")
            print("    not. The injection-rate 'transition' is a SUPPLY / RESIDENCE-TIME")
            print("    effect, not a Damkohler transition -- the scaling is being applied")
            print("    inconsistently and should be reframed accordingly.")
    print("\nDone.")


if __name__ == "__main__":
    main()
