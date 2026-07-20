#!/usr/bin/env python3
"""
analyse_grid_sweep.py  —  grid-resolution & front-structure verification

Reads the study-06 runs (S1 dissolved) and answers two questions the
baseline mesh could not:

  (1) Is the integrated carbonation grid-converged and Pe-invariant?
      -> mean carbonate VF across the convergence ladder (2/1/0.5 m at D=1e-9)
         and across the Pe sub-sweep (0.5 m at D=1e-9/1e-8/1e-7).

  (2) Does Pe set the WIDTH (sharpness) of the carbonate front?
      -> z-averaged carbonate profile carb(x) on the resolved 0.5 m grid for
         each Pe; the 10-90% leading-edge width is measured and compared.
         If width ~ sqrt(2 D t) (scales with D), Pe controls front sharpness;
         if width is ~constant, morphology is Pe-insensitive too.

Either outcome is decisive: integrated magnitude is expected to stay flat
(confirming the supply-limited result at resolution), while the front-width
behaviour settles the morphology claim that the 2 m mesh could not resolve.

Usage:
    python3 analyse_grid_sweep.py
    # BASE_DIR=/path/to/scratch python3 analyse_grid_sweep.py
    # render-test with synthetic data (no HDF5 needed):
    python3 analyse_grid_sweep.py --selftest
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

BASE_DIR = Path(os.environ.get("BASE_DIR", Path.home() / "WAG" / "grid-resolution"))
OUT_DIR = Path(os.environ.get("OUT", "."))
CARB_MINERALS = ["Calcite", "Magnesite", "Siderite", "Dolomite-ord"]
SECONDS_30YR = 30 * 365.25 * 86400

# Near-well cell widths per grid (must match generate_grid_decks.py)
GRID_WIDTHS = {
    "2m":   [2.0] * 50 + [3.0] * 50 + [9.0] * 50 + [26.0] * 50,
    "1m":   [1.0] * 100 + [3.0] * 50 + [9.0] * 50 + [26.0] * 50,
    "0p5m": [0.5] * 200 + [3.0] * 50 + [9.0] * 50 + [26.0] * 50,
}
CELL_SIZE = {"2m": 2.0, "1m": 1.0, "0p5m": 0.5}
DMULT = {"1x": 1.0e-9, "10x": 1.0e-8, "100x": 1.0e-7}

# Molar volume (m^3/mol) and CO2 per formula unit, for converting carbonate
# volume fraction -> moles -> sequestered CO2 mass.
# NOTE: VF = moles x molar_volume INSIDE PFLOTRAN, so for an exact inversion
# these molar volumes must match the MOLAR_VOLUME entries in hanford.dat. The
# values below are the standard ones; verify against the database with e.g.
#   grep -iA3 'Calcite\|Magnesite\|Siderite\|Dolomite' hanford.dat
MINERAL_PROPS = {
    "Calcite":      (3.6934e-5, 1),   # CaCO3
    "Magnesite":    (2.8018e-5, 1),   # MgCO3
    "Siderite":     (2.9378e-5, 1),   # FeCO3
    "Dolomite-ord": (6.4365e-5, 2),   # CaMg(CO3)2
}
M_CO2 = 0.04401            # kg / mol CO2
FARFIELD_X = 1500.0        # cells beyond this x (m) are unreacted -> give VF0


def x_centers(grid_key):
    w = np.array(GRID_WIDTHS[grid_key])
    edges = np.concatenate([[0.0], np.cumsum(w)])
    return 0.5 * (edges[:-1] + edges[1:])


def read_carb(sim_dir, grid_key):
    """Return (vw_mean, total_vol_m3, uw_mean, profile_x) at the final time.
    vw_mean is the VOLUME-weighted mean carbonate VF (grid-robust on the graded
    mesh); uw_mean is the unweighted cell mean (kept only for reference, since it
    over-weights the refined near-well zone). None on failure."""
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
            field = None
            for m in CARB_MINERALS:
                ks = [k for k in g.keys() if k.startswith(f"{m}_VF")]
                if ks:
                    arr = np.array(g[ks[0]], dtype=float)   # shape (NX, NY, NZ)
                    field = arr if field is None else field + arr
            if field is None:
                return None
            wx = np.array(GRID_WIDTHS[grid_key])            # Delta-x per column
            cellvol = np.broadcast_to(wx[:, None, None] * 1.0 * 2.0, field.shape)
            total_vol = float((field * cellvol).sum())      # m^3 carbonate
            vw_mean = float((field * cellvol).sum() / cellvol.sum())
            uw_mean = float(field.mean())
            profile = field.mean(axis=tuple(range(1, field.ndim)))  # z-avg vs x
            return vw_mean, total_vol, uw_mean, profile
    except Exception as exc:  # noqa: BLE001
        print(f"  warn: {sim_dir.name}: {exc}")
        return None


def carbonate_inventory(sim_dir, grid_key):
    """Net NEW carbonate at the final time, per mineral and as sequestered CO2.

    new volume_m = integral over domain of (VF_m - VF0_m) dV   [m^3 mineral]
    Baseline VF0_m is the true t=0 snapshot (per cell) when the HDF5 contains
    one; otherwise it falls back to the far-field value (x > FARFIELD_X), which
    is unreacted over 30 yr. Converts to CO2 via MINERAL_PROPS.
    Returns dict (incl. 'baseline' tag) or None."""
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
            t_first = float(tgs[0].replace("Time:", "").strip().split()[0])
            g_init = f[tgs[0]] if (len(tgs) >= 2 and t_first < 0.01) else None
            baseline = "t=0 snapshot (per cell)" if g_init is not None \
                else f"far-field x>{FARFIELD_X:g} m (t=0 not in file)"
            wx = np.array(GRID_WIDTHS[grid_key])
            xc = x_centers(grid_key)
            ff = xc > FARFIELD_X                                  # far-field x-mask
            per_mineral = {}
            total_new_vol = 0.0
            total_co2_mol = 0.0
            for m, (Vm, nco2) in MINERAL_PROPS.items():
                ks = [k for k in g.keys() if k.startswith(f"{m}_VF")]
                if not ks:
                    continue
                fld = np.array(g[ks[0]], dtype=float)            # (NX, NY, NZ)
                cellvol = np.broadcast_to(wx[:, None, None] * 1.0 * 2.0, fld.shape)
                if g_init is not None:
                    ki = [k for k in g_init.keys() if k.startswith(f"{m}_VF")]
                    vf0_field = np.array(g_init[ki[0]], dtype=float) if ki else 0.0
                    new_vol = float(((fld - vf0_field) * cellvol).sum())
                    vf0_rep = float(np.mean(vf0_field)) if ki else 0.0
                else:
                    vf0_rep = float(fld[ff, :, :].mean()) if ff.any() else 0.0
                    new_vol = float(((fld - vf0_rep) * cellvol).sum())
                per_mineral[m] = {"new_vol_m3": new_vol, "VF0": vf0_rep}
                total_new_vol += new_vol
                total_co2_mol += (new_vol / Vm) * nco2
            if not per_mineral:
                return None
            return {
                "baseline": baseline,
                "per_mineral": per_mineral,
                "new_carbonate_vol_m3": total_new_vol,
                "CO2_mol": total_co2_mol,
                "CO2_kg": total_co2_mol * M_CO2,
                "CO2_tonnes": total_co2_mol * M_CO2 / 1000.0,
            }
    except Exception as exc:  # noqa: BLE001
        print(f"  warn (inventory): {sim_dir.name}: {exc}")
        return None


def front_width(xc, prof, near_max=220.0):
    """10-90% width of the carbonate leading edge in the near-well zone.
    Returns (width_m, x_front_m, background, peak)."""
    sel = xc <= near_max
    x = xc[sel]
    p = np.asarray(prof)[sel]
    if x.size < 4:
        return None, None, None, None
    background = float(np.median(p[-max(3, x.size // 6):]))  # outer edge of window
    peak = float(p.max())
    if peak - background <= 0:
        return None, None, background, peak
    n = np.clip((p - background) / (peak - background), 0.0, 1.0)
    # profile decreases outward: find x where it crosses 0.9 and 0.1
    def cross(level):
        for i in range(1, len(n)):
            if (n[i - 1] - level) * (n[i] - level) <= 0 and n[i - 1] != n[i]:
                return float(np.interp(level, [n[i], n[i - 1]], [x[i], x[i - 1]]))
        return None
    x90, x10 = cross(0.9), cross(0.1)
    if x90 is None or x10 is None:
        return None, None, background, peak
    return abs(x10 - x90), 0.5 * (x10 + x90), background, peak


def collect():
    runs = {}
    if not BASE_DIR.is_dir():
        print(f"  (no {BASE_DIR})")
        return runs
    for sub in sorted(BASE_DIR.iterdir()):
        m = re.match(r"grid_(\w+?)_D(\w+)$", sub.name)
        if not (sub.is_dir() and m):
            continue
        gkey, dkey = m.group(1), m.group(2)
        if gkey not in GRID_WIDTHS or dkey not in DMULT:
            continue
        res = read_carb(sub, gkey)
        if res is None:
            print(f"  (no usable HDF5 in {sub.name})")
            continue
        vw_mean, total_vol, uw_mean, prof = res
        runs[(gkey, dkey)] = {"mean": vw_mean, "vol": total_vol,
                              "uw_mean": uw_mean, "profile": prof,
                              "inv": carbonate_inventory(sub, gkey)}
    return runs


# ---------------------------------------------------------------------
def make_figure(runs):
    import matplotlib as mpl
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.size"] = 9

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    axA, axB, axC = axes

    # (A) Convergence: mean carb vs near-well cell size at D=1e-9 -----------
    conv = [(CELL_SIZE[g], runs[(g, "1x")]["mean"])
            for g in ("2m", "1m", "0p5m") if (g, "1x") in runs]
    if conv:
        conv.sort(reverse=True)
        cs = [c for c, _ in conv]
        mv = [m for _, m in conv]
        axA.plot(cs, mv, "o-", color="#1f4e9c", ms=7, mfc="white", mew=1.4)
        spread = 100 * (max(mv) - min(mv)) / np.mean(mv)
        axA.set_xlabel("near-well cell size (m)")
        axA.set_ylabel("mean carbonate VF at 30 yr")
        axA.set_title(f"(a) grid convergence (D=1e-9)\nspread {spread:.1f}% across 2$\\to$0.5 m",
                      fontsize=9)
        axA.set_ylim(np.mean(mv) * 0.9, np.mean(mv) * 1.1)
        axA.invert_xaxis()
        axA.grid(True, ls=":", lw=0.5, alpha=0.5)

    # (B) Pe-invariance at the finest grid ---------------------------------
    pe = [(DMULT[d], runs[("0p5m", d)]["mean"])
          for d in ("1x", "10x", "100x") if ("0p5m", d) in runs]
    if pe:
        pe.sort()
        ds = [d for d, _ in pe]
        mv = [m for _, m in pe]
        axB.plot(ds, mv, "s-", color="#c0392b", ms=7, mfc="white", mew=1.4)
        spread = 100 * (max(mv) - min(mv)) / np.mean(mv)
        axB.set_xscale("log")
        axB.set_xlabel("diffusion coefficient $D$ (m$^2$/s)  [$Pe \\propto 1/D$]")
        axB.set_ylabel("mean carbonate VF at 30 yr")
        axB.set_title(f"(b) Pe-invariance at 0.5 m\nspread {spread:.1f}% over $Pe\\sim$10$^5$-10$^7$",
                      fontsize=9)
        axB.set_ylim(np.mean(mv) * 0.9, np.mean(mv) * 1.1)
        axB.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)

    # (C) Front profiles + width vs Pe at the finest grid ------------------
    xc = x_centers("0p5m")
    colors = {"1x": "#1f4e9c", "10x": "#7d3c98", "100x": "#c0392b"}
    widths = {}
    for d in ("1x", "10x", "100x"):
        if ("0p5m", d) not in runs:
            continue
        prof = np.asarray(runs[("0p5m", d)]["profile"])
        sel = xc <= 220.0
        axC.plot(xc[sel], prof[sel], "-", color=colors[d], lw=1.6,
                 label=f"$D$={DMULT[d]:.0e}  ($\\sqrt{{2Dt}}$={math.sqrt(2*DMULT[d]*SECONDS_30YR):.1f} m)")
        w, xf, bg, pk = front_width(xc, prof)
        if w is not None:
            widths[d] = w
    axC.set_xlabel("distance along flow $x$ (m)")
    axC.set_ylabel("z-averaged carbonate VF at 30 yr")
    title = "(c) carbonate front vs $Pe$ (0.5 m grid)"
    if len(widths) >= 2:
        wtxt = ", ".join(f"{d}:{widths[d]:.1f} m" for d in ("1x", "10x", "100x") if d in widths)
        title += f"\n10-90% width  {wtxt}"
    axC.set_title(title, fontsize=9)
    axC.legend(fontsize=7, loc="upper right")
    axC.grid(True, ls=":", lw=0.5, alpha=0.5)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"fig_grid_resolution.{ext}",
                    dpi=(300 if ext == "pdf" else 200), bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {OUT_DIR/'fig_grid_resolution.pdf'}")
    return widths


# ---------------------------------------------------------------------
def _synthetic():
    """Synthetic runs: magnitude grid-converged & Pe-invariant; front width
    set by sqrt(2 D t) (the resolved-diffusion expectation)."""
    runs = {}
    base = 8.93e-4
    seed = 4.0e-4
    # synthetic inventory: ~13 m^3 new carbonate (calcite-dominated), grid-flat
    def fake_inv():
        vol = 13.1
        co2 = (vol / MINERAL_PROPS["Calcite"][0]) * 1.0
        return {"per_mineral": {"Calcite": {"new_vol_m3": vol, "VF0": seed}},
                "new_carbonate_vol_m3": vol, "CO2_mol": co2,
                "CO2_kg": co2 * M_CO2, "CO2_tonnes": co2 * M_CO2 / 1000.0,
                "baseline": "synthetic"}
    for g in ("2m", "1m", "0p5m"):
        # volume-weighted mean is flat (converged); unweighted grows ~+49% as
        # the near-well zone is refined (the graded-mesh artifact)
        runs[(g, "1x")] = {"mean": base,
                           "uw_mean": base * (1 + 0.245 * math.log2(2.0 / CELL_SIZE[g])),
                           "vol": 20.4, "profile": None, "inv": fake_inv()}
    for d in ("1x", "10x", "100x"):
        D = DMULT[d]
        xc = x_centers("0p5m")
        w = math.sqrt(2 * D * SECONDS_30YR)
        xf = 60.0
        prof = seed + (base * 2.2 - seed) * 0.5 * np.array(
            [math.erfc((x - xf) / max(w, 0.6)) for x in xc])
        runs[("0p5m", d)] = {"mean": base, "uw_mean": base * 1.49,
                             "vol": 20.4, "profile": prof, "inv": fake_inv()}
    return runs


def main():
    selftest = "--selftest" in sys.argv
    print("=" * 60)
    print("  Grid-resolution & front-structure verification")
    print(f"  Source: {'synthetic self-test' if selftest else BASE_DIR}")
    print("=" * 60)

    runs = _synthetic() if selftest else collect()
    if not runs:
        raise SystemExit("No results found. Run the sweep first, or use --selftest.")

    # Magnitude robustness (VOLUME-WEIGHTED; unweighted shown only as the
    # graded-mesh artifact it is)
    conv = {g: runs[(g, "1x")]["mean"] for g in ("2m", "1m", "0p5m") if (g, "1x") in runs}
    if len(conv) >= 2:
        vals = list(conv.values())
        print("\n[magnitude] grid convergence at D=1e-9 (volume-weighted mean carbonate VF):")
        for g in ("2m", "1m", "0p5m"):
            if g in conv:
                uw = runs[(g, "1x")].get("uw_mean")
                uws = f"   [unweighted {uw:.4e}]" if uw is not None else ""
                print(f"    {g:>5} (cell {CELL_SIZE[g]} m): {conv[g]:.4e}{uws}")
        spread = (max(vals) - min(vals)) / np.mean(vals)
        print(f"    volume-weighted spread: {100*spread:.2f}% "
              f"-> {'GRID-CONVERGED' if spread < 0.05 else 'NOT converged'}")
        uwv = [runs[(g, "1x")].get("uw_mean") for g in conv if runs[(g, "1x")].get("uw_mean")]
        if len(uwv) >= 2:
            uspread = (max(uwv) - min(uwv)) / np.mean(uwv)
            print(f"    (unweighted cell-mean spread {100*uspread:.0f}% is a graded-mesh "
                  f"artifact: refining the near-well zone adds cells there and re-weights "
                  f"the average; the physical field is unchanged.)")

    pe = {d: runs[("0p5m", d)]["mean"] for d in ("1x", "10x", "100x") if ("0p5m", d) in runs}
    if len(pe) >= 2:
        vals = list(pe.values())
        print("\n[magnitude] Pe-invariance at 0.5 m (volume-weighted):")
        for d in ("1x", "10x", "100x"):
            if d in pe:
                print(f"    D={DMULT[d]:.0e}: {pe[d]:.4e}")
        spread = (max(vals) - min(vals)) / np.mean(vals)
        print(f"    spread: {100*spread:.2f}% "
              f"-> {'Pe-invariant' if spread < 0.05 else 'Pe-DEPENDENT'}")

    print("\nGenerating figure...")
    widths = make_figure(runs)

    if len(widths) >= 2:
        print("\n[morphology] carbonate front 10-90% width vs Pe (0.5 m grid):")
        for d in ("1x", "10x", "100x"):
            if d in widths:
                print(f"    D={DMULT[d]:.0e} (sqrt(2Dt)={math.sqrt(2*DMULT[d]*SECONDS_30YR):.1f} m): "
                      f"width {widths[d]:.1f} m")
        ws = [widths[d] for d in ("1x", "10x", "100x") if d in widths]
        ratio = max(ws) / min(ws)
        if ratio > 2:
            print(f"    front broadens {ratio:.1f}x with D -> Pe DOES control front sharpness "
                  f"(width tracks sqrt(2Dt)); magnitude unaffected.")
        else:
            print(f"    width nearly constant ({ratio:.1f}x) -> morphology also Pe-insensitive "
                  f"in the resolved range.")

    # Extensive headline: net NEW carbonate and sequestered CO2 (grid-converged)
    invs = {g: runs[(g, "1x")].get("inv") for g in ("2m", "1m", "0p5m")
            if (g, "1x") in runs and runs[(g, "1x")].get("inv")}
    if invs:
        fin = invs.get("0p5m") or list(invs.values())[-1]
        print(f"\n[inventory] net NEW carbonate and sequestered CO2 (extensive; "
              f"baseline: {fin.get('baseline', 'n/a')}):")
        print(f"    {'cell':>6} {'new_carb_m3':>13} {'CO2_tonnes':>12}")
        for g in ("2m", "1m", "0p5m"):
            iv = invs.get(g)
            if iv:
                print(f"    {CELL_SIZE[g]:>5}m {iv['new_carbonate_vol_m3']:>13.4e} "
                      f"{iv['CO2_tonnes']:>12.3f}")
        vols = [iv["new_carbonate_vol_m3"] for iv in invs.values()]
        spread = (max(vols) - min(vols)) / np.mean(vols) if np.mean(vols) else 0.0
        print(f"    grid spread: {100*spread:.2f}% "
              f"-> {'GRID-CONVERGED' if spread < 0.05 else 'NOT converged'}")
        print("    per-mineral new volume at finest grid (+ precipitation / - net loss):")
        for m, d in fin["per_mineral"].items():
            tag = "" if d["new_vol_m3"] >= 0 else "   (net dissolution vs baseline)"
            print(f"      {m.split('-')[0]:<10} {d['new_vol_m3']:>+11.4e} m3{tag}")
        print(f"    -> NET {fin['CO2_tonnes']:.2f} t CO2 mineralized per injector over 30 yr "
              f"(use this extensive value, not the domain-mean VF).")
        print("    (CO2 conversion uses standard molar volumes; cross-check against hanford.dat.)")

    payload = {
        "magnitude_convergence": {g: conv.get(g) for g in conv},
        "magnitude_vs_Pe": {DMULT[d]: pe.get(d) for d in pe},
        "front_width_m": {DMULT[d]: widths.get(d) for d in widths},
        "carbonate_inventory": {g: invs.get(g) for g in invs} if invs else {},
    }
    if not selftest:
        with open(OUT_DIR / "grid_results.json", "w") as fh:
            json.dump(payload, fh, indent=2, default=str)
        print(f"\n  saved {OUT_DIR/'grid_results.json'}")
    print("\nDone.")


if __name__ == "__main__":
    main()
