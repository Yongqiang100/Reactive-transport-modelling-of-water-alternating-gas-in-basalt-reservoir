#!/usr/bin/env python3
"""
compare_to_paper.py — compute the headline quantities from the completed WAG
runs and print them next to the values claimed in the manuscript.

Run on Setonix from the package root, with co2conv ACTIVE (provides h5py+numpy):
    source $MYSCRATCH/conda/bin/activate $MYSCRATCH/conda/envs/co2conv
    python3 compare_to_paper.py --inventory    # FIRST: dump HDF5 structure + -mas.dat header
    python3 compare_to_paper.py | tee compare_to_paper.txt

Cell volumes are read from each run's HDF5 geometry (so any grid works), with a
fallback to the 250x1x50 production widths. CO2 ("Air") inventory is read from the
'Global Air Mass in Liquid/Gas Phase' columns of the -mas.dat.

The manuscript numbers below are the ORIGINAL values (30-yr continuous, graded
140x1x25). The runs are the REVISED setup (refined 250x1x50, 100 yr = 30-yr
injection + 70-yr monitoring), so differences are expected.
"""
import os, sys, glob
from pathlib import Path
import numpy as np

try:
    import h5py
    HAVE_H5 = True
except ImportError:
    HAVE_H5 = False

ROOT = Path(os.environ.get("WAG_ROOT", Path(__file__).resolve().parent))

CARB = ["Calcite", "Magnesite", "Siderite", "Dolomite-ord"]
GRID_WIDTHS = np.array([1.0] * 100 + [3.0] * 50 + [9.0] * 50 + [26.0] * 50)  # production x widths
DZ = 2.0
MOLARV = {"Calcite": 3.6934e-5, "Magnesite": 2.8018e-5,
          "Siderite": 2.9378e-5, "Dolomite-ord": 6.4365e-5}     # m^3/mol
CO2_STOICH = {"Calcite": 1, "Magnesite": 1, "Siderite": 1, "Dolomite-ord": 2}
M_CO2 = 0.04401                                                 # kg/mol

PAPER = {"S1_carb_vf": 9.0e-4, "S2_carb_vf": 4.3e-4, "WAG_carb_vf": 6.5e-4,
         "pH": 8.4, "eta_pct": (2.6, 2.9), "S1_cum_CO2_kt": 17.3,
         "rate_ratio_30x": 16.0, "Da_sigma_opt": (3.0, 10.0)}

BASE = {"S1 dissolved": "base_dissolved", "S2 scCO2": "base_scco2",
        "S3 WAG-6mo": "base_wag6mo", "S4 WAG-3mo": "base_wag3mo",
        "S5 SWAG": "base_swag", "S6 adaptive": "base_adaptive"}
SCEN = ["dissolved", "scco2", "wag6mo", "wag3mo", "swag", "adaptive"]
MU = {"0p3": 0.3, "1": 1.0, "3": 3.0, "10": 10.0, "30": 30.0}

# --- injected CO2 from the deck FLOW spec. The reported aqueous-carbon (HCO3-) columns are not in
#     consistent units, so injected is taken from rate x CO2-mole-fraction x molar density x time.
#     gas-phase CO2 is read from the reliable 'well Air [kg]' column at run time; the LIQUID-phase
#     dissolved CO2 (which 'well Air' misses) is computed from each scenario's water-slug spec:
#       (q_l [m^3/s], x_CO2 in liquid, total liquid-injection years within the 30-yr window)
RHO_LIQ = 55556.0                                              # mol/m^3 (water-dominated liquid)
SEC_PER_YR = 365.25 * 86400.0
LIQ_INJ = {"base_dissolved": (1.0e-5, 0.04, 30.0),
           "base_scco2":     (0.0,    0.0,  0.0),
           "base_wag6mo":    (1.0e-5, 0.04, 15.0),
           "base_wag3mo":    (1.0e-5, 0.04, 15.0),
           "base_swag":      (5.0e-6, 0.35, 30.0),             # x=0.35 in liquid is supersaturated -> flagged
           "base_adaptive":  (1.0e-5, 0.04, 11.83)}


def liquid_co2_kg(stem):
    q, x, t = LIQ_INJ.get(stem, (0.0, 0.0, 0.0))
    return q * RHO_LIQ * x * M_CO2 * (t * SEC_PER_YR)


def rdir(study, name): return ROOT / study / "runs" / name
def find_h5(d):
    hs = [h for h in sorted(glob.glob(str(d / "*.h5"))) if not h.endswith("-restart.h5")]
    return hs[-1] if hs else None
def find_mas(d):
    ms = sorted(glob.glob(str(d / "*-mas*.dat")))
    return ms[0] if ms else None
def time_groups(f):
    return sorted([g for g in f.keys() if g.startswith("Time")],
                  key=lambda s: float(s.replace("Time:", "").strip().split()[0]))
def final_time_yr(sim_dir):
    """Last snapshot time (yr) in the run's HDF5, or None."""
    h5 = find_h5(sim_dir)
    if not h5:
        return None
    try:
        with h5py.File(h5, "r") as f:
            tg = time_groups(f)
            return float(tg[-1].replace("Time:", "").strip().split()[0]) if tg else None
    except Exception:  # noqa: BLE001
        return None
def run_status(sim_dir):
    """('DONE'|'TIMEOUT'|'running?'|'no log') from run.log."""
    log = sim_dir / "run.log"
    if not log.exists():
        logs = sorted(glob.glob(str(sim_dir / "*.log")))
        if not logs:
            return "no log"
        log = Path(logs[-1])
    txt = log.read_text(errors="ignore")
    if "Wall Clock Time" in txt:
        return "DONE"
    if "TIME LIMIT" in txt or "CANCELLED" in txt or "DUE TO TIME" in txt:
        return "TIMEOUT"
    return "running?"
def get_dset(group, base):
    for k in group.keys():
        if k.startswith(base):
            return np.array(group[k], dtype=float)
    return None
def carb_total(group):
    fld = None
    for m in CARB:
        a = get_dset(group, f"{m}_VF")
        if a is not None:
            fld = a if fld is None else fld + a
    return fld


def cell_vol_from_file(f, shape):
    """Cell-volume array matching `shape`, read from the HDF5 grid geometry
    (Coordinates/Domain group: X/Y/Z cell-edge arrays). Falls back to the
    250x1x50 production widths, then to uniform (with a warning)."""
    nx, ny, nz = (list(shape) + [1, 1, 1])[:3]
    grp = None
    for g in ("Coordinates", "Domain", "Grid"):
        if g in f and isinstance(f[g], h5py.Group):
            grp = f[g]; break
    if grp is not None:
        def widths(letter, n):
            ks = [k for k in grp.keys() if k.strip().lower().startswith(letter)]
            if not ks:
                return None
            a = np.array(grp[ks[0]], dtype=float).ravel()
            if a.size == n + 1:
                return np.diff(a)          # cell edges -> widths
            if a.size == n:
                return a                   # already per-cell
            return None
        dx, dy, dz = widths("x", nx), widths("y", ny), widths("z", nz)
        if dx is not None and dz is not None:
            if dy is None:
                dy = np.ones(ny)
            vol = dx[:, None, None] * dy[None, :, None] * dz[None, None, :]
            if vol.shape == tuple(shape):
                return vol
    if tuple(shape) == (250, 1, 50):
        return GRID_WIDTHS[:, None, None] * 1.0 * DZ * np.ones(shape)
    print(f"    warn: no grid geometry in HDF5 for shape {tuple(shape)}; uniform weights used")
    return np.ones(shape)


def vf_metrics(sim_dir):
    """(dom_mean_VF, vol_weighted_VF, injection_driven_m3, phi0, phiF, pH) or None."""
    h5 = find_h5(sim_dir)
    if not (HAVE_H5 and h5):
        return None
    try:
        with h5py.File(h5, "r") as f:
            tg = time_groups(f)
            if not tg:
                return None
            g1, g0 = f[tg[-1]], f[tg[0]]
            final = carb_total(g1)
            if final is None:
                return None
            t0 = carb_total(g0)
            if t0 is None:
                t0 = np.zeros_like(final)
            cv = cell_vol_from_file(f, final.shape)
            m1 = float(final.mean())
            vfw = float((final * cv).sum() / cv.sum())
            m3 = float(((final - t0) * cv).sum())
            phi = get_dset(g1, "Porosity"); phi0 = get_dset(g0, "Porosity")
            ph = get_dset(g1, "pH")
            if ph is None:
                ph = get_dset(g1, "ph")
            return (m1, vfw, m3,
                    float(phi0.mean()) if phi0 is not None else None,
                    float(phi.mean()) if phi is not None else None,
                    float(ph.mean()) if ph is not None else None)
    except Exception as exc:  # noqa: BLE001
        print(f"    warn: {sim_dir.name}: {exc}")
        return None


def mineralized_CO2_kg(sim_dir):
    h5 = find_h5(sim_dir)
    if not (HAVE_H5 and h5):
        return None
    try:
        with h5py.File(h5, "r") as f:
            tg = time_groups(f)
            g1, g0 = f[tg[-1]], f[tg[0]]
            tot, cv = 0.0, None
            for m in CARB:
                a1 = get_dset(g1, f"{m}_VF")
                if a1 is None:
                    continue
                a0 = get_dset(g0, f"{m}_VF")
                if a0 is None:
                    a0 = np.zeros_like(a1)
                if cv is None:
                    cv = cell_vol_from_file(f, a1.shape)
                mol = (a1 - a0) * cv / MOLARV[m]
                tot += float((mol * CO2_STOICH[m] * M_CO2).sum())
            return tot
    except Exception:  # noqa: BLE001
        return None


def read_mas(sim_dir):
    """Final CO2 ('Air') inventory: dict(mobile, liq, gas, injected_or_None).
    Columns: 'Global Air Mass in Liquid Phase', '... Gas Phase'. Header is
    comma/quote-delimited; data rows are whitespace-delimited."""
    p = find_mas(sim_dir)
    if not p:
        return None
    try:
        lines = [ln for ln in open(p).read().splitlines() if ln.strip()]
        cols = [c.strip().strip('"').strip() for c in lines[0].split(",")]

        def col(must, exclude=()):
            for i, c in enumerate(cols):
                cl = c.lower()
                if all(s in cl for s in must) and not any(e in cl for e in exclude):
                    return i
            return None
        i_liq = col(["air", "liquid"])                        # Global Air in liquid (mobile, dissolved)
        i_gas = col(["air", "gas"])                           # Global Air in gas (mobile, free)
        i_well = col(["well", "air", "kg"], exclude=["yr"])   # cumulative CO2 injected at well (Air component)
        if i_well is None:
            i_well = col(["injection", "air", "kg"], exclude=["yr"])
        i_out = col(["outlet", "air", "kg"], exclude=["yr"])  # cumulative CO2 vented at outlet (Air)
        # aqueous inorganic carbon (DIC), tracked as the HCO3- component [mol]:
        i_well_dic = col(["well", "hco3", "mol"], exclude=["yr"])
        i_out_dic = col(["outlet", "hco3", "mol"], exclude=["yr"])
        i_glob_dic = col(["global", "hco3", "mol"], exclude=["yr"])
        rows = []
        for ln in lines[1:]:
            try:
                rows.append([float(x) for x in ln.split()])
            except ValueError:
                continue
        if not rows:
            return None
        last, first = rows[-1], rows[0]

        def v(i, row=None):
            row = last if row is None else row
            return float(row[i]) if (i is not None and i < len(row)) else None
        liq, gas = v(i_liq), v(i_gas)
        mobile = None
        if liq is not None or gas is not None:
            mobile = (liq or 0.0) + (gas or 0.0)
        return {"mobile": mobile, "liq": liq, "gas": gas,
                "well": v(i_well), "vented": v(i_out),
                "well_dic": v(i_well_dic), "vented_dic": v(i_out_dic),
                "glob_dic": v(i_glob_dic), "glob_dic0": v(i_glob_dic, first)}
    except Exception as exc:  # noqa: BLE001
        print(f"    warn mas {sim_dir.name}: {exc}")
        return None


def fit_p(mu, y):
    mu = np.asarray(mu, float); y = np.asarray(y, float)
    ok = (mu > 0) & (y > 0)
    return float(np.polyfit(np.log(mu[ok]), np.log(y[ok]), 1)[0]) if ok.sum() >= 2 else float("nan")
def fmt(x, e=False):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "n/a"
    return f"{x:.3e}" if e else f"{x:.3f}"
def fmt_kt(x):
    return "n/a" if x is None else f"{x/1e6:.2f} kt"
def fmt_t(x):
    return "n/a" if x is None else f"{x/1e3:6.1f} t"


def _need_h5():
    print("h5py is not available in this Python.")
    print("Activate co2conv first (this is what provides h5py + numpy):")
    print("  source $MYSCRATCH/conda/bin/activate $MYSCRATCH/conda/envs/co2conv")
    print("If h5py is missing from co2conv: pip install h5py  (or conda install -c conda-forge h5py)")


def inventory():
    for study in ["01_baseline", "06_grid_resolution", "08_rate_sweep"]:
        rd = ROOT / study / "runs"
        if not rd.is_dir():
            continue
        for d in sorted(rd.glob("*")):
            h5 = find_h5(d)
            if not h5:
                continue
            print(f"run.h5 inventory  ({d}/{Path(h5).name})")
            with h5py.File(h5, "r") as f:
                print("  top-level items:")
                for k in f.keys():
                    it = f[k]
                    if isinstance(it, h5py.Group):
                        sub = list(it.keys())
                        ex = sub[:6]
                        print(f"    [group] {k}/  ({len(sub)})  e.g. {ex}")
                    else:
                        print(f"    [dset]  {k}  shape={getattr(it,'shape',None)}")
                tg = time_groups(f)
                print(f"  time groups ({len(tg)}): {tg[0]}  ...  {tg[-1]}")
                cell = None
                ct = carb_total(f[tg[-1]])
                if ct is not None:
                    cv = cell_vol_from_file(f, ct.shape)
                    print(f"  carbonate field shape {ct.shape}; bulk volume from geometry = {cv.sum():.1f} m^3")
            mp = find_mas(d)
            if mp:
                hdr = [c.strip().strip('"').strip() for c in open(mp).readline().split(",")]
                print(f"\n  -mas.dat columns ({len(hdr)}):")
                for i, c in enumerate(hdr):
                    print(f"    [{i:2d}] {c}")
            return
    print("No readable run.h5 found under */runs/. Have the jobs finished?")


def carbon_check():
    """Prove that 'Global HCO3- [mol]' is the TOTAL dissolved-carbon component (DIC) -- i.e. it
    already sums CO2(aq)+HCO3-+CO3--+carbonate complexes -- by cross-checking it against the
    domain integral of the HDF5 Total_HCO3- field over the aqueous (water) volume."""
    for stem in ("base_dissolved", "base_scco2", "base_wag6mo"):
        d = rdir("01_baseline", stem)
        h5, mp = find_h5(d), find_mas(d)
        if not (h5 and mp):
            continue
        print(f"carbon-species check  ({d.name})")
        with h5py.File(h5, "r") as f:
            g = f[time_groups(f)[-1]]
            tots = sorted(k for k in g.keys() if k.startswith("Total_"))
            carbon_basis = [k for k in tots if "hco3" in k.lower() or "co3" in k.lower() or "co2" in k.lower()]
            print(f"  aqueous basis species (Total_*): {len(tots)} total")
            print(f"  carbon-bearing basis species: {carbon_basis or '[none other than HCO3-]'}")
            print("    (CO2(aq), CO3--, CaHCO3+, CaCO3(aq), MgHCO3+ ... are SECONDARY species,")
            print("     summed into the HCO3- component total -- so Total_HCO3- = DIC)")
            dic = get_dset(g, "Total_HCO3-")           # mol/L (analytical/total component)
            por = get_dset(g, "Porosity")
            sg = get_dset(g, "Gas_Saturation")
            if dic is None or por is None:
                print("  (Total_HCO3- or Porosity missing — cannot integrate)\n"); continue
            cv = cell_vol_from_file(f, dic.shape)
            sl = (1.0 - sg) if sg is not None else np.ones_like(dic)
            vwater_L = por * sl * cv * 1000.0          # m^3 -> L
            dic_h5 = float((dic * vwater_L).sum())     # mol
            # diagnostics: actual field magnitudes (to explain any mismatch)
            print(f"  Total_HCO3- [field, mol/L]  min={dic.min():.3e} mean={dic.mean():.3e} max={dic.max():.3e}")
            print(f"  Porosity mean={por.mean():.3f}   Gas_Saturation min/mean/max="
                  f"{(sg.min() if sg is not None else 0):.3f}/{(sg.mean() if sg is not None else 0):.3f}/{(sg.max() if sg is not None else 0):.3f}")
            print(f"  water volume = {vwater_L.sum():.3e} L ; bulk volume = {cv.sum():.1f} m^3")
        lines = [ln for ln in open(mp).read().splitlines() if ln.strip()]
        cols = [c.strip().strip('"').strip() for c in lines[0].split(",")]
        gi = next((i for i, c in enumerate(cols)
                   if "global" in c.lower() and "hco3" in c.lower() and "yr" not in c.lower()), None)
        dic_mb = None
        if gi is not None:
            row = [float(x) for x in lines[-1].split()]
            dic_mb = row[gi] if gi < len(row) else None
        print(f"  DIC from HDF5 Total_HCO3- integral = {dic_h5:.4e} mol")
        if dic_mb is not None:
            print(f"  DIC from mass-balance Global HCO3- = {dic_mb:.4e} mol")
            print(f"  ratio massbal/HDF5 = {dic_mb/dic_h5:.3f}  ->  ~1.0 confirms Global HCO3- = total DIC")
        print()
    print("If the ratios are ~1, the carbon budget already counts ALL aqueous carbon species,")
    print("so DIC = HCO3- component is correct (not just free bicarbonate).")


def rate_status():
    """Report each rate-sweep run's last snapshot time vs the 100 yr target and log
    completion, so incomplete (still-advancing) high-rate runs are visible."""
    target = 100.0
    print(f"Rate-sweep run status (target final time = {target:.0f} yr):")
    print(f"  {'run':<22} {'last_t(yr)':>10} {'status':>10}   note")
    n_incomplete = 0
    for sc in SCEN:
        for tag, mu in MU.items():
            d = rdir("08_rate_sweep", f"rs_{sc}_mu{tag}")
            if not d.is_dir():
                continue
            ft = final_time_yr(d)
            st = run_status(d)
            incomplete = (ft is None) or (ft < target - 1.0)
            if incomplete:
                n_incomplete += 1
            note = "INCOMPLETE -> ratio provisional" if incomplete else ""
            ft_s = f"{ft:.1f}" if ft is not None else "NA"
            print(f"  rs_{sc}_mu{tag:<4} {ft_s:>10} {st:>10}   {note}")
    print(f"\n  {n_incomplete} run(s) not yet at {target:.0f} yr.")
    if n_incomplete:
        print("  Resubmit to let them finish (resume guard skips DONE runs):")
        print("    WALLTIME=23:00:00 NTASKS_PER_RUN=32 ./run_all_setonix.sh   # 08 stragglers")
        print("  Then re-run compare_to_paper.py; the [3] ratios for those scenarios will settle.")


def main():
    if not HAVE_H5:
        _need_h5(); return
    if "--inventory" in sys.argv:
        inventory(); return
    if "--carbon-check" in sys.argv:
        carbon_check(); return
    if "--rate-status" in sys.argv:
        rate_status(); return

    print("=" * 78)
    print("  WAG CO2 — computed results vs manuscript claims")
    print("  (paper = ORIGINAL 30-yr/coarse values; runs = revised 100-yr/refined)")
    print("=" * 78)

    print("\n[1] Baseline (01_baseline) — final total-carbonate volume fraction")
    print(f"  {'scenario':<26} {'paper(dom)':>11} {'computed':>11}   vol-wtd | phi0->phiF | pH")
    for label, stem in BASE.items():
        mm = vf_metrics(rdir("01_baseline", stem))
        if mm is None:
            print(f"  {label:<26} {'(unreadable)':>23}"); continue
        m1, vfw, m3, phi0, phiF, ph = mm
        paper = (PAPER["S1_carb_vf"] if "S1" in label else PAPER["S2_carb_vf"] if "S2" in label
                 else PAPER["WAG_carb_vf"] if "WAG" in label else None)
        pa = f"{paper:.2e}" if paper else "—"
        print(f"  {label:<26} {pa:>11} {fmt(m1, e=True):>11}   "
              f"{fmt(vfw, e=True)} | {fmt(phi0)}->{fmt(phiF)} | {fmt(ph)}")

    print("\n[2] Carbon budget — injected CO2 from the FLOW spec (unit-clean), eta = mineralized / injected")
    print("    gas-phase CO2 = well Air [kg] (reliable); liquid-phase dissolved CO2 = q_l x x_CO2 x rho_liq x t")
    print(f"  {'scenario':<13} {'mineralzd':>10} {'inj_gas':>9} {'inj_liq':>9} {'injected':>9} {'eta':>7}")
    for label, stem in BASE.items():
        d = rdir("01_baseline", stem)
        mas = read_mas(d); minz = mineralized_CO2_kg(d)
        if minz is None:
            print(f"  {label:<13} (unreadable)"); continue
        inj_gas = abs(mas["well"]) if (mas and mas["well"] is not None) else 0.0   # gas-phase CO2 (kg)
        inj_liq = liquid_co2_kg(stem)                                             # liquid-phase dissolved CO2 (kg)
        inj = inj_gas + inj_liq
        eta = (100.0 * minz / inj) if inj > 0 else None
        flag = "  (x_CO2=0.35 liquid supersaturated; injected suspect)" if stem == "base_swag" else ""
        print(f"  {label:<13} {fmt_t(minz)} {fmt_t(inj_gas)} {fmt_t(inj_liq)} {fmt_t(inj)} {fmt(eta):>6}%{flag}")
    print(f"  (paper eta {PAPER['eta_pct'][0]}-{PAPER['eta_pct'][1]} %. Injected is computed from the deck flow")
    print("   condition because the reported aqueous-carbon (HCO3-) columns are not in consistent units;")
    print("   eta must be <=100%, and injected>=mineralized confirms this budget over the HCO3--based one.)")

    print("\n[3] Rate sweep (08_rate_sweep) — carbonate scaling with injection rate")
    print(f"  {'scenario':<12} {'p(dom)':>8} {'p(inj)':>8} {'30x dom':>9} {'30x inj':>9} {'min_t(yr)':>10}")
    any_incomplete = False
    for sc in SCEN:
        mus, m1s, m3s, fts = [], [], [], []
        for tag, mu in MU.items():
            sd = rdir("08_rate_sweep", f"rs_{sc}_mu{tag}")
            mm = vf_metrics(sd)
            if mm:
                mus.append(mu); m1s.append(mm[0]); m3s.append(mm[2])
                ft = final_time_yr(sd)
                if ft is not None:
                    fts.append(ft)
        if len(mus) < 2:
            print(f"  {sc:<12} {'(insufficient)':>8}"); continue
        o = np.argsort(mus); mus = np.array(mus)[o]; m1s = np.array(m1s)[o]; m3s = np.array(m3s)[o]
        def ratio(a):
            b = a[np.isclose(mus, 1.0)]; t = a[np.isclose(mus, 30.0)]
            return float(t[0] / b[0]) if (len(b) and len(t) and b[0] > 0) else float("nan")
        min_t = min(fts) if fts else float("nan")
        incomplete = not (min_t >= 99.0)
        any_incomplete = any_incomplete or incomplete
        flag = " *" if incomplete else ""
        print(f"  {sc:<12} {fit_p(mus, m1s):>8.2f} {fit_p(mus, m3s):>8.2f} "
              f"{ratio(m1s):>9.2f} {ratio(m3s):>9.2f} {min_t:>10.1f}{flag}")
    print(f"  (paper ~{PAPER['rate_ratio_30x']:.0f}x, domain-mean. dom & inj agreeing => genuine, not a metric artifact.)")
    if any_incomplete:
        print("  * min_t < 100 yr: at least one rate point is still advancing -> that ratio is PROVISIONAL.")
        print("    Run  python3 compare_to_paper.py --rate-status  to see which, then resubmit to finish.")

    print("\n[4] Grid resolution (06_grid_resolution) — convergence of the two metrics")
    print(f"  {'grid dir':<28} {'dom-mean VF':>12} {'vol-wtd VF':>12} {'inj-driven m3':>14}")
    rd = ROOT / "06_grid_resolution" / "runs"
    if rd.is_dir():
        for d in sorted(rd.glob("*")):
            mm = vf_metrics(d)
            if mm:
                print(f"  {d.name:<28} {fmt(mm[0], e=True):>12} {fmt(mm[1], e=True):>12} {fmt(mm[2], e=True):>14}")
    print("  (vol-weighted/injection-driven converging while domain-mean drifts = metric artifact.)")

    print("\n" + "=" * 78)
    print("  Finer per-study detail (03 Da-Pe, 04 mechanisms, 05 kappa, 07 Da-consistency):")
    print("  run each analyse_*.py with BASE_DIR=<study>/runs.")
    print("=" * 78)


if __name__ == "__main__":
    main()
