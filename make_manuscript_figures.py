#!/usr/bin/env python3
"""
make_manuscript_figures.py — regenerate the manuscript figure set from the completed
WAG runs, reading the HDF5 output directly with the same validated readers as
compare_to_paper.py (no precomputed JSON summaries needed).

Run on Setonix from the package root with co2conv active:
    source $MYSCRATCH/conda/bin/activate $MYSCRATCH/conda/envs/co2conv
    python3 -c "import matplotlib" || pip install matplotlib
    python3 make_manuscript_figures.py            # writes ./figures/*.pdf (+ .png)

Emits the manuscript's eight figures (matching \\includegraphics filenames) plus the
grid-convergence revision figure:
    fig_domain_schematic         model domain (illustration)
    fig_baseline_comparison      6-panel time series, all six scenarios
    fig_carbonate_breakdown      carbonate-phase split per scenario
    fig_spatial_profiles         radial pH / porosity / carbonate profiles
    fig_gas_saturation_2d        2-D gas-saturation cross-sections (shared scale)
    fig_damkohler_sweep          carbonate VF and dphi vs injection-rate multiplier
    fig_da_sigma_regime          Da-Sigma regime diagram
    fig_kinetic_sensitivity      Palandri/Kharaka vs Rimstidt (auto-discovered runs)
    fig_grid_convergence         metric convergence vs grid + Pe independence (revision)

All carbonate masses come from the carbonate volume fractions; geometry (cell volumes,
near-well / far-field masks, map axes) is read from each run's own HDF5 Coordinates.
"""
import os, sys, glob, re
from pathlib import Path
from collections import defaultdict
import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("h5py missing — activate co2conv:  source $MYSCRATCH/conda/bin/activate $MYSCRATCH/conda/envs/co2conv")
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Rectangle
    from matplotlib.colors import LinearSegmentedColormap, LogNorm
except ImportError:
    sys.exit("matplotlib missing — with co2conv active:  pip install matplotlib")

ROOT = Path(os.environ.get("WAG_ROOT", Path(__file__).resolve().parent))
OUT = ROOT / "figures"; OUT.mkdir(exist_ok=True)

# ---------- validated constants ----------
CARB = ["Calcite", "Magnesite", "Siderite", "Dolomite-ord"]
PRIM = ["Forsterite", "Anorthite", "Diopside"]
GRID_WIDTHS = np.array([1.0] * 100 + [3.0] * 50 + [9.0] * 50 + [26.0] * 50)
DZ = 2.0
NW_X, FF_X = 100.0, 700.0          # near-well <= 100 m ; far-field >= 700 m
INJECT_END = 30.0                   # yr (injection window)
CC = {"Calcite": "#1f77b4", "Magnesite": "#2ca02c", "Siderite": "#ff7f0e", "Dolomite-ord": "#9467bd"}

BASE = [("base_dissolved", "Dissolved (S1)", "#1f77b4", "o"),
        ("base_scco2", "scCO$_2$ (S2)", "#d62728", "s"),
        ("base_wag6mo", "WAG-6mo (S3)", "#2ca02c", "^"),
        ("base_wag3mo", "WAG-3mo (S4)", "#9467bd", "v"),
        ("base_swag", "SWAG (S5)", "#ff7f0e", "D"),
        ("base_adaptive", "Adaptive (S6)", "#17becf", "P")]

# Distinct linestyles so the near-identical WAG curves (S3/S4/S5), which share the
# same carbonate volume fraction, remain visible where they overlap in line plots.
LS_BASE = {"base_wag3mo": "--", "base_swag": ":"}  # others default to solid ("-")
RSC = [("dissolved", "Dissolved", "#1f77b4", "o"), ("scco2", "scCO$_2$", "#d62728", "s"),
       ("wag6mo", "WAG-6mo", "#2ca02c", "^"), ("wag3mo", "WAG-3mo", "#9467bd", "v"),
       ("swag", "SWAG", "#ff7f0e", "D"), ("adaptive", "Adaptive", "#17becf", "P")]
MU = {"0p3": 0.3, "1": 1.0, "3": 3.0, "10": 10.0, "30": 30.0}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 8, "axes.titlesize": 8, "axes.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6, "legend.frameon": False,
    "lines.linewidth": 1.0, "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.05})
LW = 1.0
_made, _skipped = [], []


def _label(ax, t, x=-0.12, y=1.08):
    ax.text(x, y, t, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf"); fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig); _made.append(name); print(f"  wrote figures/{name}.pdf (+ .png)")


# ---------- HDF5 readers ----------
def rdir(study, name): return ROOT / study / "runs" / name
def find_h5(d):
    hs = [h for h in sorted(glob.glob(str(Path(d) / "*.h5"))) if not h.endswith("-restart.h5")]
    return hs[-1] if hs else None
def find_mas(d):
    ms = sorted(glob.glob(str(Path(d) / "*-mas*.dat")))
    return ms[0] if ms else None
def _tg(f):
    return sorted([g for g in f.keys() if g.startswith("Time")],
                  key=lambda s: float(s.replace("Time:", "").strip().split()[0]))
def _tyr(g): return float(g.replace("Time:", "").strip().split()[0])
def get_dset(group, base):
    for k in list(group.keys()):
        if k.startswith(base):
            try:
                return np.array(group[k], dtype=float)
            except Exception:
                continue
    return None
def carb_total(group):
    fld = None
    for m in CARB:
        a = get_dset(group, f"{m}_VF")
        if a is not None:
            fld = a if fld is None else fld + a
    return fld
def edges(f, letter, n):
    for g in ("Coordinates", "Domain", "Grid"):
        if g in f and isinstance(f[g], h5py.Group):
            for k in f[g].keys():
                if k.strip().lower().startswith(letter):
                    a = np.array(f[g][k], dtype=float).ravel()
                    if a.size == n + 1:
                        return a
    return None
def _axes(f, nx, nz):
    xe = edges(f, "x", nx)
    ze = edges(f, "z", nz)
    if xe is None:
        xe = np.concatenate([[0], np.cumsum(GRID_WIDTHS)]) if nx == len(GRID_WIDTHS) else np.arange(nx + 1, dtype=float)
    if ze is None:
        ze = np.arange(nz + 1, dtype=float) * DZ
    return xe, ze
def cell_vol(f, shape):
    nx, ny, nz = (list(shape) + [1, 1, 1])[:3]
    xe, ze = _axes(f, nx, nz)
    ye = edges(f, "y", ny)
    dy = np.diff(ye) if ye is not None else np.array([1.0])
    return np.diff(xe)[:, None, None] * dy[None, :, None] * np.diff(ze)[None, None, :]


def read_run(d):
    """One pass over a run: time series + final spatial profile + final 2-D fields.
    Robust to unreadable individual time groups (skips them)."""
    h5 = find_h5(d)
    if not h5:
        return None
    try:
        f = h5py.File(h5, "r")
    except Exception:
        return None
    with f:
        tg = _tg(f)
        if not tg:
            return None
        # geometry from first readable carbonate field
        shape = None
        for gn in tg:
            try:
                c = carb_total(f[gn])
            except Exception:
                c = None
            if c is not None:
                shape = c.shape; break
        if shape is None:
            return None
        nx, ny, nz = (list(shape) + [1, 1, 1])[:3]
        xe, ze = _axes(f, nx, nz)
        xc = 0.5 * (xe[:-1] + xe[1:])
        nw = xc <= NW_X
        ff = xc >= FF_X
        cv = cell_vol(f, shape)
        o = defaultdict(list); o["t"] = []
        last_ok = None
        for gn in tg:
            try:
                g = f[gn]; c = carb_total(g)
            except Exception:
                continue
            if c is None:
                continue
            o["t"].append(_tyr(gn))
            o["carb_mean"].append(float(c.mean()))
            o["carb_vw"].append(float((c * cv).sum() / cv.sum()))
            o["carb_nw"].append(float(c[nw].mean()))
            ph = get_dset(g, "pH");  ph = get_dset(g, "ph") if ph is None else ph
            o["ph_mean"].append(float(ph.mean()) if ph is not None else np.nan)
            o["ph_nw"].append(float(ph[nw].mean()) if ph is not None else np.nan)
            poro = get_dset(g, "Porosity")
            o["poro_mean"].append(float(poro.mean()) if poro is not None else np.nan)
            o["poro_nw"].append(float(poro[nw].mean()) if poro is not None else np.nan)
            gas = get_dset(g, "Gas_Saturation")
            o["gas_mean"].append(float(gas.mean()) if gas is not None else np.nan)
            fo = get_dset(g, "Forsterite_VF")
            o["fo_mean"].append(float(fo.mean()) if fo is not None else np.nan)
            last_ok = gn
        if not o["t"]:
            return None
        for k in list(o.keys()):
            o[k] = np.array(o[k], dtype=float)
        # final-time spatial profile (mid-depth) + 2-D fields + phase split
        g = f[last_ok]
        kmid = nz // 2
        cfull = carb_total(g)
        prof = {"x": xc, "carb": cfull[:, 0, kmid].copy()}
        ph = get_dset(g, "pH");  ph = get_dset(g, "ph") if ph is None else ph
        if ph is not None:
            prof["ph"] = ph[:, 0, kmid].copy()
        poro = get_dset(g, "Porosity")
        if poro is not None:
            prof["phi"] = poro[:, 0, kmid].copy()
        o2 = {"t": o["t"], "series": o, "prof": prof, "xe": xe, "ze": ze,
              "carb2d": cfull[:, 0, :].copy(),
              "phases": {m: float(np.nanmean(get_dset(g, f"{m}_VF"))) if get_dset(g, f"{m}_VF") is not None else 0.0
                         for m in CARB}}
        gas = get_dset(g, "Gas_Saturation")
        o2["gas2d"] = gas[:, 0, :].copy() if gas is not None else None
        return o2


def well_air_kg(d):
    p = find_mas(d)
    if not p:
        return 0.0
    lines = [ln for ln in open(p).read().splitlines() if ln.strip()]
    cols = [c.strip().strip('"').strip() for c in lines[0].split(",")]
    gi = next((i for i, c in enumerate(cols)
               if "well" in c.lower() and "air" in c.lower() and "kg" in c.lower() and "yr" not in c.lower()), None)
    if gi is None:
        return 0.0
    for ln in reversed(lines[1:]):
        try:
            row = [float(x) for x in ln.split()]
            return abs(row[gi]) if gi < len(row) else 0.0
        except ValueError:
            continue
    return 0.0


def field2d_at(d, base, target_yr):
    """Return (xe, ze, field[nx,nz], t_actual) for `base` at the time group
    nearest target_yr — used for transient fields (e.g. gas saturation, which
    peaks at the end of injection and dissolves during monitoring)."""
    h5 = find_h5(d)
    if not h5:
        return None
    try:
        f = h5py.File(h5, "r")
    except Exception:
        return None
    with f:
        tg = _tg(f)
        if not tg:
            return None
        gn = min(tg, key=lambda g: abs(_tyr(g) - target_yr))
        try:
            a = get_dset(f[gn], base)
        except Exception:
            return None
        if a is None:
            return None
        nx, ny, nz = (list(a.shape) + [1, 1, 1])[:3]
        xe, ze = _axes(f, nx, nz)
        return xe, ze, a[:, 0, :].copy(), _tyr(gn)


# ======================================================================
def fig_domain():
    fig, ax = plt.subplots(figsize=(7.2, 2.8))
    z1, z2, z3, Z = 3, 6, 12, 5
    ax.add_patch(Rectangle((0, 0), z1, Z, fc="#fef2f2", ec="none"))
    ax.add_patch(Rectangle((z1, 0), z2 - z1, Z, fc="#fefce8", ec="none"))
    ax.add_patch(Rectangle((z2, 0), z3 - z2, Z, fc="#eff6ff", ec="none"))
    ax.plot([0, z3, z3, 0, 0], [0, 0, Z, Z, 0], color="#334155", lw=0.8)
    for xb in (z1, z2):
        ax.plot([xb, xb], [0, Z], color="#94a3b8", ls="--", lw=0.5)
    ax.add_patch(Rectangle((0, 1.0), 0.3, 3.0, fc="#fca5a5", ec="#dc2626", lw=1.0))
    akw = dict(arrowstyle="->", mutation_scale=10)
    for z in (2.0, 3.0):
        ax.annotate("", xy=(2.5, z), xytext=(0.4, z), arrowprops={**akw, "color": "#dc2626", "lw": 0.9})
        ax.annotate("", xy=(5.5, z), xytext=(3.2, z), arrowprops={**akw, "color": "#d97706", "lw": 0.8, "alpha": 0.7})
    ax.annotate("", xy=(11.0, 2.5), xytext=(6.3, 2.5), arrowprops={**akw, "color": "#3b82f6", "lw": 0.7, "alpha": 0.4})
    ax.plot([z3, z3], [0, Z], color="#2563eb", lw=2.5)
    ax.text(z1 / 2, Z - 0.3, "Reactive", fontsize=7, ha="center", va="top", color="#dc2626", fontstyle="italic")
    ax.text((z1 + z2) / 2, Z - 0.3, "Transition", fontsize=7, ha="center", va="top", color="#d97706", fontstyle="italic")
    ax.text((z2 + z3) / 2, Z - 0.3, "Near-equilibrium", fontsize=7, ha="center", va="top", color="#2563eb", fontstyle="italic")
    ax.text(z1 / 2, 2.5, "pH 4.5–6", fontsize=6, ha="center", va="center", color="#dc2626", alpha=0.6)
    ax.text((z1 + z2) / 2, 2.5, "pH 6–7.5", fontsize=6, ha="center", va="center", color="#d97706", alpha=0.6)
    ax.text((z2 + z3) / 2, 2.0, "pH 8+, Q/K ≈ 1", fontsize=6, ha="center", va="center", color="#2563eb", alpha=0.6)
    ax.text(z1 / 2, -0.4, "0–100 m\nΔx = 1 m, 100 cells", fontsize=5.5, ha="center", va="top", color="#64748b")
    ax.text((z1 + z2) / 2, -0.4, "100–700 m\nΔx = 3–9 m, 100 cells", fontsize=5.5, ha="center", va="top", color="#64748b")
    ax.text((z2 + z3) / 2, -0.4, "700–2000 m\nΔx = 26 m, 50 cells", fontsize=5.5, ha="center", va="top", color="#64748b")
    ax.text(-0.5, 4.0, "z = 80 m", fontsize=6, ha="right", va="center", color="#64748b")
    ax.text(-0.5, 1.0, "z = 20 m", fontsize=6, ha="right", va="center", color="#64748b")
    ax.text(-0.5, 2.5, "60 °C", fontsize=7.5, ha="right", va="center", color="#dc2626", fontweight="bold")
    ax.text(0.15, Z + 0.3, "Well", fontsize=7, ha="center", va="bottom", fontweight="bold", color="#dc2626")
    ax.text(z3, Z + 0.3, "x = 2 km", fontsize=6, ha="right", va="bottom", color="#64748b")
    ax.text(z3 + 0.3, 2.5, "P = 6 MPa\nOutlet BC", fontsize=7, ha="left", va="center", color="#2563eb", fontweight="bold")
    ax.set_xlim(-1.2, z3 + 1.8); ax.set_ylim(-1.5, Z + 0.8); ax.axis("off")
    save(fig, "fig_domain_schematic")


def fig_baseline():
    data = {}
    for stem, lbl, c, mk in BASE:
        r = read_run(rdir("01_baseline", stem))
        if r:
            data[stem] = r
    if not data:
        _skipped.append("fig_baseline_comparison (no 01_baseline runs)"); return
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.5), sharex=True)
    fig.subplots_adjust(hspace=0.16, wspace=0.42, bottom=0.18)
    panels = [("carb_mean", "Total carbonate VF"), ("ph_mean", "Mean pH"),
              ("poro_nw", "Near-well porosity"), ("fo_mean", "Forsterite VF"),
              ("gas_mean", "Gas saturation"), ("carb_nw", "Near-well carbonate VF")]
    leg = []
    for ax, (key, yl), lab in zip(axes.flat, panels, "abcdef"):
        for stem, lbl, c, mk in BASE:
            if stem not in data:
                continue
            s = data[stem]["series"]; t = s["t"]; y = s[key]
            if np.all(np.isnan(y)):
                continue
            l, = ax.plot(t, y, color=c, lw=LW, ls=LS_BASE.get(stem, "-"), label=lbl)
            if lab == "a":
                leg.append(l)
        ax.axvline(INJECT_END, ls=":", c="0.5", lw=0.7)
        ax.set_ylabel(yl, fontsize=7); _label(ax, lab)
        if lab in "abc":
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("Time (yr)")
        if key in ("carb_mean", "carb_nw", "fo_mean"):
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    fig.legend(handles=leg, ncol=6, fontsize=6, loc="lower center",
               bbox_to_anchor=(0.5, 0.02), frameon=False, columnspacing=1.0, handlelength=1.5)
    save(fig, "fig_baseline_comparison")


def fig_carb_bar():
    labels, cd = [], {m: [] for m in CARB}
    for stem, lbl, c, mk in BASE:
        r = read_run(rdir("01_baseline", stem))
        if not r:
            continue
        labels.append(lbl.split(" (")[0])
        for m in CARB:
            cd[m].append(r["phases"].get(m, 0.0))
    if not labels:
        _skipped.append("fig_carbonate_breakdown (no 01_baseline runs)"); return
    fig, ax = plt.subplots(figsize=(5.5, 2.8)); fig.subplots_adjust(bottom=0.22)
    x = np.arange(len(labels)); bot = np.zeros(len(labels))
    for m in CARB:
        v = np.array(cd[m]); ax.bar(x, v, 0.5, bottom=bot, label=m, color=CC[m], edgecolor="white", linewidth=0.3); bot += v
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=0, ha="center", fontsize=6.5)
    ax.set_ylabel("Volume fraction"); ax.legend(fontsize=5.5, ncol=4, loc="upper right")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    save(fig, "fig_carbonate_breakdown")


def fig_spatial():
    fig, axes = plt.subplots(1, 3, figsize=(8.5, 2.9)); fig.subplots_adjust(wspace=0.34, bottom=0.26)
    leg = []
    any_data = False
    for stem, lbl, c, mk in BASE:
        r = read_run(rdir("01_baseline", stem))
        if not r:
            continue
        any_data = True
        p = r["prof"]; x = p["x"]
        if "ph" in p:
            l, = axes[0].plot(x, p["ph"], color=c, lw=LW, label=lbl); leg.append(l)
        if "phi" in p:
            axes[1].plot(x, p["phi"], color=c, lw=LW)
        if "carb" in p:
            axes[2].plot(x, p["carb"], color=c, lw=LW)
    if not any_data:
        plt.close(fig); _skipped.append("fig_spatial_profiles (no 01_baseline runs)"); return
    for ax in axes:
        ax.set_xlim(0, 200); ax.set_xlabel("Distance from well (m)")
    axes[2].ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    for ax, yl, lab in zip(axes, ["pH", "Porosity", "Carbonate VF"], "abc"):
        ax.set_ylabel(yl); _label(ax, lab, y=1.12)
    fig.legend(handles=leg, ncol=6, fontsize=6, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), frameon=False, columnspacing=1.0, handlelength=1.5)
    save(fig, "fig_spatial_profiles")


def fig_gas2d():
    # Gas saturation is transient: peaks at end of injection (~30 yr) and dissolves
    # during monitoring. The manuscript shows the end-of-injection state, not 100 yr.
    target = INJECT_END
    fields = {}; gmax = 0.0; tshown = None
    for stem, lbl, c, mk in BASE:
        res = field2d_at(rdir("01_baseline", stem), "Gas_Saturation", target)
        if res is not None:
            xe, ze, g2d, ta = res
            fields[stem] = (xe, ze, g2d, lbl)
            gmax = max(gmax, float(np.nanmax(g2d)))
            tshown = ta if tshown is None else tshown
    if not fields:
        _skipped.append("fig_gas_saturation_2d (no Gas_Saturation fields)"); return
    vmax = max(gmax, 1e-3)
    vmin = vmax / 1e3                       # ~3 decades of dynamic range
    norm = LogNorm(vmin=vmin, vmax=vmax)
    cmap = LinearSegmentedColormap.from_list(
        "blue_red", ["#e3f2fd", "#e8f5e9", "#fff9c4", "#ffcc80", "#e53935"], N=256).copy()
    cmap.set_bad("#e3f2fd"); cmap.set_under("#e3f2fd")   # zero / sub-vmin gas -> light blue (no gas)
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.4), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.30, wspace=0.10, bottom=0.13, right=0.87, top=0.88)
    for ax, (stem, lbl, c, mk), lab in zip(axes.flat, BASE, "abcdef"):
        ax.set_facecolor("#e3f2fd")
        if stem in fields:
            xe, ze, sg, _ = fields[stem]
            sgm = np.ma.masked_less_equal(sg.T, 0.0)
            ax.pcolormesh(xe, ze, sgm, cmap=cmap, norm=norm, shading="auto")
        else:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center", color="#94a3b8")
        ax.set_title(lbl, fontsize=7, pad=3); _label(ax, lab, y=1.13)
        ax.axhline(20, color="white", lw=0.6, ls="--", alpha=0.8)
        ax.axhline(80, color="white", lw=0.6, ls="--", alpha=0.8)
        ax.set_xlim(0, 1500); ax.set_ylim(0, 100)  # widened to 1500 m: plume bulk ~800 m, WAG peaks ~1440 m; buoyant gas -> top; no inversion
    fig.text(0.43, 0.04, "Distance from well (m)", fontsize=8)
    fig.text(0.06, 0.5, "Elevation, z (m)", rotation="vertical", fontsize=8)
    cax = fig.add_axes([0.89, 0.13, 0.015, 0.75])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, extend="min")
    cb.set_label("Gas saturation (log scale)", fontsize=7); cb.ax.tick_params(labelsize=6)
    save(fig, "fig_gas_saturation_2d")


def _rate_curve(get_final):
    out = {}
    for sc, lbl, c, mk in RSC:
        pts = []
        for tag, mu in MU.items():
            r = read_run(rdir("08_rate_sweep", f"rs_{sc}_mu{tag}"))
            if r:
                v = get_final(r)
                if v is not None and np.isfinite(v):
                    pts.append((mu, v))
        if len(pts) >= 2:
            pts.sort(); out[sc] = (np.array([p[0] for p in pts]), np.array([p[1] for p in pts]), lbl, c, mk)
    return out


def fig_da_sweep():
    carb = _rate_curve(lambda r: r["series"]["carb_mean"][-1])
    dphi = _rate_curve(lambda r: (r["series"]["poro_mean"][-1] - r["series"]["poro_mean"][0])
                       if np.isfinite(r["series"]["poro_mean"][0]) else None)
    if not carb:
        _skipped.append("fig_damkohler_sweep (no 08_rate_sweep runs)"); return
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9)); fig.subplots_adjust(wspace=0.34, bottom=0.2)
    for sc, (mus, vs, lbl, c, mk) in carb.items():
        axes[0].plot(mus, vs, color=c, lw=LW, marker=mk, ms=4, mfc="none", mec=c, mew=1.0, label=lbl)
    axes[0].set_xscale("log"); axes[0].set_yscale("log")
    axes[0].set_xlabel("Rate multiplier ($Q/Q_0$)"); axes[0].set_ylabel("Carbonate VF (100 yr)")
    _label(axes[0], "a"); axes[0].legend(fontsize=5, ncol=2)
    for sc, (mus, vs, lbl, c, mk) in dphi.items():
        axes[1].plot(mus, vs, color=c, lw=LW, marker=mk, ms=4, mfc="none", mec=c, mew=1.0, label=lbl)
    axes[1].set_xscale("log"); axes[1].axhline(0, color="#94a3b8", lw=0.4)
    axes[1].set_xlabel("Rate multiplier ($Q/Q_0$)"); axes[1].set_ylabel(r"$\Delta\phi\ (\phi_{100}-\phi_0)$")
    _label(axes[1], "b")
    for ax in axes:
        ax.set_xticks([0.3, 1, 3, 10, 30]); ax.set_xticklabels(["0.3", "1", "3", "10", "30"])
    save(fig, "fig_damkohler_sweep")


def fig_da_sigma():
    R_, T_, Vm = 8.314, 333.15, 3.69e-5
    sig3 = 0.7 * 2900 * 9.81 * 600; p0 = 1000 * 9.81 * 600
    kOm = 1e5  # Omega proxy scaling, recalibrated for zero-seeded carbonate VF
    dphi = _rate_curve(lambda r: (r["series"]["poro_mean"][-1] - r["series"]["poro_mean"][0])
                       if np.isfinite(r["series"]["poro_mean"][0]) else None)
    curves = _rate_curve(lambda r: r["series"]["carb_mean"][-1])
    if not curves:
        _skipped.append("fig_da_sigma_regime (no 08_rate_sweep runs)"); return
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2)); fig.subplots_adjust(wspace=0.34, bottom=0.24)
    # panel (a): Da vs porosity change -- direct self-sealing regime map
    for sc, (mus, vs, lbl, c, mk) in dphi.items():
        da = 10.0 / mus
        axes[0].plot(da, vs, color=c, lw=0.6, alpha=0.5)
        axes[0].scatter(da, vs, c=c, s=26, marker=mk, edgecolors="white", linewidths=0.4, label=lbl)
    axes[0].set_xscale("log"); axes[0].invert_xaxis(); axes[0].axhline(0, color="#94a3b8", lw=0.4)
    axes[0].set_xlabel("Damk\u00f6hler number (Da)"); axes[0].set_ylabel(r"$\Delta\phi$ (30 yr)")
    _label(axes[0], "a"); axes[0].legend(fontsize=5, ncol=2)
    # panel (b): Da vs normalised stress ratio Sigma -- crystallization-pressure proxy
    for sc, (mus, vs, lbl, c, mk) in curves.items():
        dav, sigv = [], []
        for mu, tc in zip(mus, vs):
            om = max(1.01, 1 + tc * kOm)
            Pc = (R_ * T_ / Vm) * np.log(om)
            se = sig3 - p0 - 3e6 * min(mu, 5)
            sigv.append(Pc / max(se, 1e5)); dav.append(10.0 / mu)
        axes[1].plot(dav, sigv, color=c, lw=0.6, alpha=0.5)
        axes[1].scatter(dav, sigv, c=c, s=26, marker=mk, edgecolors="white", linewidths=0.4, label=lbl)
    axes[1].set_xscale("log"); axes[1].set_yscale("log"); axes[1].invert_xaxis()
    axes[1].set_xlabel("Damk\u00f6hler number (Da)"); axes[1].set_ylabel(r"Normalised stress ratio ($\Sigma$)")
    _label(axes[1], "b")
    save(fig, "fig_da_sigma_regime")


def _discover_kinetic():
    """Find two run dirs corresponding to the two kinetic parameterizations."""
    pats = [("pk", ("pk", "palandri", "kharaka")), ("rim", ("rim", "rimstidt"))]
    for study in ("05_kinetic_crossover", "04_mechanisms", "03_dape"):
        rd = ROOT / study / "runs"
        if not rd.is_dir():
            continue
        names = [p.name for p in rd.glob("*") if p.is_dir()]
        found = {}
        for key, toks in pats:
            for n in names:
                if any(tok in n.lower() for tok in toks):
                    found[key] = rd / n; break
        if "pk" in found and "rim" in found:
            return found["pk"], found["rim"]
    return None, None


def fig_kin_sens():
    pk_dir, rm_dir = _discover_kinetic()
    if pk_dir is None:
        _skipped.append("fig_kinetic_sensitivity (could not locate Palandri/Rimstidt runs — "
                        "paste `ls 05_kinetic_crossover/runs` and I will wire it up)")
        return
    pk, rm = read_run(pk_dir), read_run(rm_dir)
    if not pk or not rm:
        _skipped.append("fig_kinetic_sensitivity (runs found but unreadable)"); return
    cases = [("Palandri & Kharaka (2004)", pk, "#1f77b4"), ("Rimstidt et al. (2012)", rm, "#d62728")]
    fig = plt.figure(figsize=(7.2, 7.2)); gs = GridSpec(3, 3, figure=fig, hspace=0.2, wspace=0.42)
    fig.subplots_adjust(bottom=0.10)
    leg = []
    def S(d): return d["series"]
    ax = fig.add_subplot(gs[0, 0])
    for lb, d, c in cases:
        l, = ax.plot(S(d)["t"], S(d)["carb_nw"], color=c, lw=LW, label=lb); leg.append(l)
    ax.set_ylabel("Volume fraction"); _label(ax, "a"); ax.tick_params(labelbottom=False)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    axb = fig.add_subplot(gs[0, 1])
    for lb, d, c in cases:
        axb.plot(S(d)["t"], S(d)["fo_mean"], color=c, lw=LW)
    axb.set_ylabel("Forsterite VF"); _label(axb, "b"); axb.tick_params(labelbottom=False)
    axc = fig.add_subplot(gs[0, 2])
    for lb, d, c in cases:
        axc.plot(S(d)["t"], S(d)["ph_nw"], color=c, lw=LW)
    axc.set_ylabel("Near-well pH"); _label(axc, "c"); axc.tick_params(labelbottom=False)
    axd = fig.add_subplot(gs[1, 0])
    for lb, d, c in cases:
        axd.plot(S(d)["t"], S(d)["carb_mean"], color=c, lw=LW)
    axd.set_ylabel("Carbonate VF (mean)"); _label(axd, "d"); axd.tick_params(labelbottom=False)
    axd.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    axe = fig.add_subplot(gs[1, 1])
    for lb, d, c in cases:
        axe.plot(S(d)["t"], S(d)["fo_mean"], color=c, lw=LW)
    axe.set_ylabel("Forsterite VF (mean)"); _label(axe, "e"); axe.tick_params(labelbottom=False)
    axf = fig.add_subplot(gs[1, 2])
    for lb, d, c in cases:
        axf.plot(S(d)["t"], S(d)["poro_nw"], color=c, lw=LW)
    axf.set_ylabel("Near-well porosity"); _label(axf, "f"); axf.tick_params(labelbottom=False)
    axg = fig.add_subplot(gs[2, 0]); axg.set_xlabel("Time (yr)"); axg.set_ylabel("Mean pH"); _label(axg, "g")
    for lb, d, c in cases:
        axg.plot(S(d)["t"], S(d)["ph_mean"], color=c, lw=LW)
    axh = fig.add_subplot(gs[2, 1]); axh.set_xlabel("Time (yr)"); axh.set_ylabel("Gas saturation"); _label(axh, "h")
    for lb, d, c in cases:
        axh.plot(S(d)["t"], S(d)["gas_mean"], color=c, lw=LW)
    axi = fig.add_subplot(gs[2, 2]); axi.axis("off")
    cm = S(cases[0][1])["carb_mean"][-1]; cr = S(cases[1][1])["carb_mean"][-1]
    cd = abs(cr - cm) / max(cm, 1e-12) * 100
    axi.text(0.05, 0.95, f"Rate ratio at pH 8.5:\n  Rimstidt / P&K ≈ 1 500×\n\n"
                         f"Carbonate VF difference\n  (100 yr): {cd:.1f}%\n\n→ Transport-limited\n  (1 − Q/K) → 0",
             transform=axi.transAxes, fontsize=6.5, va="top", fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.4", fc="#f8fafc", ec="#cbd5e1", lw=0.5))
    _label(axi, "i")
    fig.legend(handles=leg, ncol=2, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 0.01),
               frameon=False, columnspacing=1.5, handlelength=2.0)
    save(fig, "fig_kinetic_sensitivity")


def _grid_parse(name):
    cs = next((v for k, v in {"0p5m": 0.5, "1m": 1.0, "2m": 2.0}.items() if k in name), np.nan)
    m = re.search(r"D(\d+)x", name)
    return cs, (int(m.group(1)) if m else np.nan)


def fig_grid():
    rd = ROOT / "06_grid_resolution" / "runs"
    if not rd.is_dir():
        _skipped.append("fig_grid_convergence (no 06_grid_resolution runs)"); return
    data = {}
    for d in sorted(rd.glob("*")):
        r = read_run(d)
        if not r:
            continue
        cs, dm = _grid_parse(d.name)
        data[d.name] = (cs, dm, r["series"]["carb_mean"][-1], r["series"]["carb_vw"][-1])
    if not data:
        _skipped.append("fig_grid_convergence (no readable grid runs)"); return
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.2, 3.0)); fig.subplots_adjust(wspace=0.3, bottom=0.18)
    conv = sorted([(v[0], v[2], v[3]) for v in data.values() if v[1] == 1 and np.isfinite(v[0])])
    if conv:
        cs = [c[0] for c in conv]
        axL.plot(cs, [c[1] for c in conv], "s-", color="#d62728", label="domain-mean")
        axL.plot(cs, [c[2] for c in conv], "o-", color="#1f77b4", label="volume-weighted")
        axL.invert_xaxis(); axL.set_xlabel("Cell size (m)"); axL.set_ylabel("Carbonate VF")
        axL.set_title("(a) Metric convergence vs grid"); axL.legend(loc="center right")
        axL.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    pe = sorted([(v[1], v[2], v[3]) for v in data.values() if v[0] == 0.5 and np.isfinite(v[1])])
    if pe:
        dm = [p[0] for p in pe]
        axR.semilogx(dm, [p[1] for p in pe], "s-", color="#d62728", label="domain-mean")
        axR.semilogx(dm, [p[2] for p in pe], "o-", color="#1f77b4", label="volume-weighted")
        axR.set_ylim(0, max(p[1] for p in pe) * 1.3)
        axR.set_xlabel("Diffusion multiplier (Pe$^{-1}$)"); axR.set_ylabel("Carbonate VF")
        axR.set_title("(b) Diffusion (Pe) independence"); axR.legend(loc="center right")
        axR.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    save(fig, "fig_grid_convergence")


def main():
    print(f"Writing figures to {OUT}/")
    for fn in (fig_domain, fig_baseline, fig_carb_bar, fig_spatial, fig_gas2d,
               fig_da_sweep, fig_da_sigma, fig_kin_sens, fig_grid):
        try:
            fn()
        except Exception as exc:
            _skipped.append(f"{fn.__name__} failed: {exc}")
            print(f"  {fn.__name__} failed: {exc}")
    print("\nSummary:")
    print(f"  produced {len(_made)} figure(s): {', '.join(_made) if _made else '(none)'}")
    if _skipped:
        print("  skipped / issues:")
        for s in _skipped:
            print(f"    - {s}")
    print("Done.")


if __name__ == "__main__":
    main()
