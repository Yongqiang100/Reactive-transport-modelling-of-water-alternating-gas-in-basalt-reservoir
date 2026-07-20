#!/usr/bin/env python3
"""
analyse_kappa_sweep.py

Reads the global kinetic-rate crossover sweep (28 runs: 14 kappa x 2
endmembers) and produces the central transport-limitation result:

  * carbonate volume fraction at 30 yr as a function of the kinetic-rate
    multiplier kappa, for S1 (dissolved) and S2 (supercritical);
  * the transport-limited plateau value (high-kappa limit);
  * the crossover kappa_crit (carbonation = 50% of plateau), i.e. the
    point where the system leaves the transport-limited regime;
  * the reaction-limited tail slope d(log carb)/d(log kappa) at low kappa
    (-> 1 confirms carbonation ~ kappa, the reaction-limited signature);
  * a figure (fig_kappa_crossover.pdf) with the literature forsterite
    kinetic range shaded, demonstrating that every credible rate sits on
    the plateau.

Because Da ~ kappa at fixed Pe, the kappa axis is a relative-Da axis: the
crossover defines Da ~ O(1), to be compared with the Da ~ 10 transition
found independently along the injection-rate axis.

Usage:
    python3 analyse_kappa_sweep.py
    # BASE_DIR=/path/to/scratch python3 analyse_kappa_sweep.py
    # (self-test rendering with synthetic data, no HDF5 needed:)
    python3 analyse_kappa_sweep.py --selftest
"""

import os
import re
import sys
import glob
import json
import math
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

BASE_DIR = Path(os.environ.get("BASE_DIR", Path.home() / "WAG" / "kinetic-crossover"))
OUT_DIR = Path(os.environ.get("OUT", "."))
CARB_MINERALS = ["Calcite", "Magnesite", "Siderite", "Dolomite-ord"]

# Reference points for context (from the manuscript headline kinetic test).
# Baseline kinetics = Rimstidt et al. (2012) -> kappa = 1.
# Palandri & Kharaka (2004) is ~1500x slower in intrinsic forsterite rate.
RIMSTIDT_KAPPA = 1.0
PALANDRI_KAPPA = 1.0 / 1500.0      # ~6.7e-4 relative to the Rimstidt baseline


# ---------------------------------------------------------------------
# HDF5 reader (same convention as the Da-Pe analysis)
# ---------------------------------------------------------------------
def read_final_carbonate(sim_dir: Path):
    """Mean total carbonate VF at the final time, summed over the four
    carbonate phases. Returns None if no usable HDF5 output is found."""
    if not HAVE_H5:
        return None
    h5s = sorted(glob.glob(str(sim_dir / "*.h5")))
    if not h5s:
        return None
    try:
        with h5py.File(h5s[-1], "r") as f:
            time_grps = sorted(
                [g for g in f.keys() if g.startswith("Time")],
                key=lambda s: float(s.replace("Time:", "").strip().split()[0]),
            )
            if not time_grps:
                return None
            g = f[time_grps[-1]]
            carb_total = 0.0
            found = False
            for mineral in CARB_MINERALS:
                matches = [k for k in g.keys() if k.startswith(f"{mineral}_VF")]
                if matches:
                    carb_total += float(np.array(g[matches[0]]).mean())
                    found = True
            return carb_total if found else None
    except Exception as exc:  # noqa: BLE001
        print(f"  warn: failed to read {sim_dir.name}: {exc}")
        return None


def parse_ktag(tag: str):
    """'k1em3' -> 1e-3, 'k3e0' -> 3.0, 'k1e2' -> 100.0"""
    m = re.match(r"k(\d+)e(m?)(\d+)$", tag)
    if not m:
        return None
    mant = int(m.group(1))
    sign = -1 if m.group(2) == "m" else 1
    exp = int(m.group(3))
    return mant * 10.0 ** (sign * exp)


def collect():
    """Walk BASE_DIR/<scenario>/<ktag>/ and read carbonate VF at 30 yr."""
    results = {"dissolved": [], "scco2": []}
    for scenario in results:
        sdir = BASE_DIR / scenario
        if not sdir.is_dir():
            print(f"  (no {sdir})")
            continue
        for sub in sorted(sdir.iterdir()):
            if not sub.is_dir():
                continue
            kappa = parse_ktag(sub.name)
            if kappa is None:
                continue
            carb = read_final_carbonate(sub)
            if carb is None:
                print(f"  (no usable HDF5 in {scenario}/{sub.name})")
                continue
            results[scenario].append((kappa, carb))
        results[scenario].sort()
    return results


# ---------------------------------------------------------------------
# Regime statistics
# ---------------------------------------------------------------------
def analyse_series(pts):
    """pts: sorted list of (kappa, carb). Returns plateau, kappa_crit,
    tail slope, and a normalized-Da crossover summary."""
    if len(pts) < 3:
        return None
    k = np.array([p[0] for p in pts], float)
    c = np.array([p[1] for p in pts], float)

    # Plateau = mean of the top decade of kappa (transport-limited limit)
    plateau = float(c[k >= k.max() / 10.0].mean())

    # Crossover: largest kappa at which carbonation falls to 50% of plateau,
    # found by linear interpolation in log-kappa.
    half = 0.5 * plateau
    kappa_crit = None
    for i in range(1, len(k)):
        if (c[i - 1] - half) * (c[i] - half) <= 0 and c[i - 1] != c[i]:
            lk = np.interp(half, [c[i - 1], c[i]],
                           [math.log10(k[i - 1]), math.log10(k[i])])
            kappa_crit = 10.0 ** lk
            break

    # Reaction-limited tail slope: log-log slope over the lowest two decades
    lo = k <= (k.min() * 100.0)
    tail_slope = None
    if lo.sum() >= 2 and (c[lo] > 0).all():
        tail_slope = float(np.polyfit(np.log10(k[lo]), np.log10(c[lo]), 1)[0])

    return {
        "plateau": plateau,
        "kappa_crit": kappa_crit,
        "tail_slope": tail_slope,
        "n_points": len(pts),
        "kappa_min": float(k.min()),
        "kappa_max": float(k.max()),
    }


# ---------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------
def make_figure(results, stats):
    import matplotlib as mpl
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.size"] = 9
    mpl.rcParams["axes.linewidth"] = 0.7
    mpl.rcParams["mathtext.fontset"] = "dejavuserif"

    fig, ax = plt.subplots(figsize=(6.2, 4.4))

    styles = {
        "dissolved": dict(color="#1f4e9c", marker="o", label="S1 dissolved ($x_{CO_2}=0.04$)"),
        "scco2":     dict(color="#c0392b", marker="s", label="S2 supercritical ($x_{CO_2}=0.99$)"),
    }

    # Literature forsterite kinetic range (shaded): Palandri (slow) -> Rimstidt (baseline)
    ax.axvspan(PALANDRI_KAPPA, RIMSTIDT_KAPPA, color="0.85", zorder=0,
               label="forsterite rate range\n(Palandri$\\to$Rimstidt, $\\sim$1500$\\times$)")

    all_carb = []
    for scenario, sty in styles.items():
        pts = results.get(scenario, [])
        if not pts:
            continue
        k = np.array([p[0] for p in pts], float)
        c = np.array([p[1] for p in pts], float)
        all_carb.extend(c.tolist())
        ax.plot(k, c, marker=sty["marker"], color=sty["color"], ms=5.5,
                lw=1.4, mfc="white", mec=sty["color"], mew=1.3,
                label=sty["label"], zorder=3)

        st = stats.get(scenario)
        if st:
            # plateau line
            ax.axhline(st["plateau"], color=sty["color"], ls=":", lw=0.8,
                       alpha=0.6, zorder=1)
            # crossover marker
            if st["kappa_crit"]:
                ax.plot([st["kappa_crit"]], [0.5 * st["plateau"]], marker="v",
                        color=sty["color"], ms=8, mfc=sty["color"], zorder=4)

    # ----- choose representation from what the data actually shows -----
    slopes = [stats[s]["tail_slope"] for s in stats
              if stats.get(s) and stats[s].get("tail_slope") is not None]
    crossed = any(stats[s].get("kappa_crit") for s in stats if stats.get(s))
    reaction_limited = crossed or any(sl is not None and sl > 0.5 for sl in slopes)
    ymin = min(all_carb) if all_carb else 1e-4
    ymax = max(all_carb) if all_carb else 1e-3

    if reaction_limited:
        # Full regime diagram: log-log, with the reaction-limited slope-1 guide
        ax.set_yscale("log")
        kref = np.array([1e-5, 1e-2])
        anchor = min(c for c in all_carb if c > 0)
        ax.plot(kref, anchor * (kref / kref[0]), color="0.4", ls="--", lw=0.9, zorder=1)
        ax.annotate("reaction-limited\n(carb $\\propto \\kappa$)", xy=(3e-4, anchor * 30),
                    color="0.35", fontsize=7.5, ha="center")
        ax.annotate("transport-limited plateau", xy=(6.0, ymax * 1.25),
                    color="0.25", fontsize=8, ha="center")
    else:
        # Carbonation is invariant: linear y zoomed to the data; no false guide.
        ax.set_yscale("linear")
        ax.set_ylim(ymin * 0.85, ymax * 1.12)
        ax.annotate("carbonation invariant across the explored range\n"
                    "reaction-limited crossover not reached "
                    "($Da_{\\mathrm{crit}}\\!\\ll$ baseline)",
                    xy=(0.5, 0.95), xycoords="axes fraction",
                    color="0.25", fontsize=8.5, ha="center", va="top")
        # flag the first sign of reaction limitation (a low-kappa dip), if any
        for scenario, sty in styles.items():
            pts = results.get(scenario, [])
            st = stats.get(scenario)
            if pts and st and st.get("plateau") and pts[0][1] < 0.9 * st["plateau"]:
                k0, c0 = pts[0]
                ax.annotate("first sign of\nreaction limitation\n($-${:.0f}% at $\\kappa={:g}$)"
                            .format(100 * (1 - c0 / st["plateau"]), k0),
                            xy=(k0, c0), xytext=(k0 * 4, c0 * 0.82),
                            color=sty["color"], fontsize=7, ha="left", va="top",
                            arrowprops=dict(arrowstyle="->", color=sty["color"], lw=0.8))

    ax.set_xscale("log")
    ax.set_xlabel("kinetic-rate multiplier $\\kappa$   "
                  "(Damk\u00f6hler number $Da \\propto \\kappa$, $Pe$ fixed)")
    ax.set_ylabel(r"carbonate volume fraction at 30 yr")
    ax.set_title(
        "Kinetic sensitivity of carbonation: transport-limited plateau and crossover"
        if reaction_limited else
        "Kinetic insensitivity of carbonation: transport-limited across the explored range",
        fontsize=9.5)
    ax.grid(True, which="major", ls="-", lw=0.3, alpha=0.3)
    ax.grid(True, which="minor", ls="-", lw=0.2, alpha=0.15)
    ax.legend(fontsize=7.2, loc="center right", framealpha=0.95)

    fig.tight_layout()
    out_pdf = OUT_DIR / "fig_kappa_crossover.pdf"
    out_png = OUT_DIR / "fig_kappa_crossover.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(out_png, dpi=200, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  saved {out_pdf}")
    print(f"  saved {out_png}")


# ---------------------------------------------------------------------
# Synthetic self-test (no HDF5): saturation model carb = P*k/(k+k_crit)
# ---------------------------------------------------------------------
def _synthetic():
    kappas = [1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2,
              1e-1, 3e-1, 1.0, 3.0, 10.0, 30.0, 100.0]
    res = {}
    for scenario, P, kc in [("dissolved", 8.93e-4, 1.5e-4), ("scco2", 4.25e-4, 2.2e-4)]:
        res[scenario] = [(k, P * k / (k + kc)) for k in kappas]
    return res


def _load_json(path):
    """Reconstruct the results dict from a saved kappa_results.json."""
    d = json.load(open(path))
    return {s: [(p["kappa"], p["carb_VF"]) for p in d["series"].get(s, [])]
            for s in ("dissolved", "scco2")}


def main():
    selftest = "--selftest" in sys.argv
    from_json = None
    if "--from-json" in sys.argv:
        from_json = sys.argv[sys.argv.index("--from-json") + 1]

    print("=" * 60)
    print("  Kinetic-rate crossover analysis")
    src = ("synthetic self-test" if selftest else
           f"saved JSON: {from_json}" if from_json else str(BASE_DIR))
    print(f"  Source: {src}")
    print("=" * 60)

    if selftest:
        results = _synthetic()
    elif from_json:
        results = _load_json(from_json)
    else:
        results = collect()
    if not any(results.get(s) for s in ("dissolved", "scco2")):
        raise SystemExit("No results found. Run the sweep first, or use --selftest.")

    stats = {}
    for scenario in ("dissolved", "scco2"):
        pts = results.get(scenario, [])
        st = analyse_series(pts)
        stats[scenario] = st
        print(f"\n--- {scenario} ({len(pts)} points) ---")
        if st:
            print(f"  transport-limited plateau : {st['plateau']:.3e}")
            kc = st["kappa_crit"]
            ts = st["tail_slope"]
            if kc:
                print(f"  crossover kappa (50% plateau): {kc:.2e}")
                print(f"  Palandri (kappa={PALANDRI_KAPPA:.1e}) is "
                      f"{'ON the plateau' if PALANDRI_KAPPA > kc else 'in the reaction-limited zone'}")
            else:
                # crossover not bracketed: report the lowest-kappa departure
                k0, c0 = pts[0]
                drop = 100 * (1 - c0 / st["plateau"])
                print(f"  crossover NOT reached: carbonation still "
                      f"{c0/st['plateau']*100:.0f}% of plateau at kappa={k0:g} "
                      f"(-{drop:.0f}% vs plateau); Da_crit lies below the explored range")
            if ts is not None:
                print(f"  low-kappa log-log slope     : {ts:.2f}  "
                      f"(->1 would indicate reaction-limited; ~0 = still transport-limited)")

    payload = {
        "kappa_ladder": sorted({k for s in results.values() for k, _ in s}),
        "series": {s: [{"kappa": k, "carb_VF": c} for k, c in results.get(s, [])]
                   for s in results},
        "stats": stats,
        "reference": {"rimstidt_kappa": RIMSTIDT_KAPPA, "palandri_kappa": PALANDRI_KAPPA},
    }
    out_json = OUT_DIR / "kappa_results.json"
    if not from_json:                       # don't overwrite the source when re-plotting
        with open(out_json, "w") as fh:
            json.dump(payload, fh, indent=2, default=str)
        print(f"\n  saved {out_json}")

    print("\nGenerating figure...")
    make_figure(results, stats)
    print("\nDone.")


if __name__ == "__main__":
    main()
