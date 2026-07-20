#!/usr/bin/env python3
"""
extract_two_timepoints.py -- report carbonation metrics at BOTH the end of
injection (t = 30 yr) and the end of monitoring (t = 100 yr) for all six
baseline scenarios, reading the 100-yr HDF5 snapshots (0.5 .. 30 .. 100 yr)
already on disk. NO PFLOTRAN re-run is needed -- the 30-yr snapshot is one of
the OUTPUT TIMES the decks already write.

It reuses the *validated* constants and low-level readers from
analyse_transport_limitation.py (identical carbonate phases, molar volumes,
CO2 stoichiometry, cell geometry, no-injection control, and injected-mass
logic), so the t = 100 column must reproduce the published numbers
(S1 ~ 20.8, S2 ~ 0.19, WAG ~ 10.4, S6 ~ 8.3 t injection-driven; S1 carbonate
VF ~ 8.2e-5). That reproduction is the built-in correctness check; the t = 30
column is the new end-of-injection data for the two-timepoint tables.

Place this file next to analyse_transport_limitation.py in the package root and
run on Setonix with co2conv active:

    source $MYSCRATCH/conda/bin/activate $MYSCRATCH/conda/envs/co2conv
    cd $MYSCRATCH/WAG            # package root (WAG_ROOT)
    python3 extract_two_timepoints.py
    # optional: choose different readout years
    python3 extract_two_timepoints.py --years 30 100

Paste the printed table back and I will reconcile the manuscript with each
number at its actual year.
"""
import argparse
from pathlib import Path
import numpy as np
import h5py

# Reuse the validated module (same directory). Everything below pulls its
# constants and readers so conventions match the published analysis exactly.
try:
    import analyse_transport_limitation as A
except Exception as e:  # pragma: no cover
    raise SystemExit(
        "Could not import analyse_transport_limitation.py -- put this script in "
        "the same directory (the package root) and activate co2conv first.\n"
        f"Import error: {e}"
    )


def _fields_at(d, target_yr):
    """Read the snapshot nearest target_yr from run dir d.

    Returns dict of per-carbonate VF fields at the target snapshot, the initial
    (first-snapshot) fields, cell volume, porosity at target and initial, the
    forsterite VF field, and the actual snapshot years used. None if unreadable.
    """
    h5 = A.find_h5(d)
    if not h5:
        return None
    try:
        f = h5py.File(h5, "r")
    except Exception:
        return None
    with f:
        tg = A._tg(f)
        if not tg:
            return None
        gsel = min(tg, key=lambda g: abs(A._t(g) - target_yr))  # nearest snapshot
        ginit = tg[0]

        def carb_map(gn):
            g = f[gn]
            out = {}
            for m in A.CARB:
                a = A.get_dset(g, f"{m}_VF")
                if a is not None:
                    out[m] = a
            return out

        cur = carb_map(gsel)
        if not cur:
            return None
        ini = carb_map(ginit)
        if not ini:
            ini = {m: np.zeros_like(v) for m, v in cur.items()}
        shape = next(iter(cur.values())).shape
        cv = A.cell_vol(f, shape)
        poro = A.get_dset(f[gsel], "Porosity")
        poro0 = A.get_dset(f[ginit], "Porosity")
        forst = A.get_dset(f[gsel], "Forsterite_VF")
        return {"carb": cur, "carb0": ini, "cv": cv, "poro": poro,
                "poro0": poro0, "forst": forst,
                "t_used": A._t(gsel), "t_init": A._t(ginit)}


def _co2_t(cur, ref, cv):
    """Control-subtracted CO2 mass in tonnes:
    sum_m (cur - ref) * cell_vol / molar_volume * CO2_stoich * M_CO2 / 1e3."""
    s = 0.0
    for m in A.CARB:
        if m in cur:
            rm = ref.get(m, np.zeros_like(cur[m]))
            s += float(((cur[m] - rm) * cv).sum()) / A.MOLARV[m] * A.CO2_STOICH[m]
    return s * A.M_CO2 / 1e3


def _carb_vf_mean(cur, cv):
    """Domain VOLUME-WEIGHTED total carbonate volume fraction (uncontrolled;
    matches the sum-Carb-VF column of the results table). Volume weighting is
    required because the mesh is strongly non-uniform (1 m near-well vs 26 m
    far-field cells) and carbonate concentrates in the small near-well cells,
    so an unweighted mean over cells overstates the field by ~6x."""
    tot = None
    for m in A.CARB:
        if m in cur:
            tot = cur[m] if tot is None else tot + cur[m]
    if tot is None:
        return float("nan")
    return float((tot * cv).sum() / cv.sum())


def _xmask(xmax):
    """Boolean mask over x-cells whose centre lies within xmax metres."""
    xe = np.concatenate([[0.0], np.cumsum(A.GRID_WIDTHS)])
    xc = 0.5 * (xe[:-1] + xe[1:])
    return xc < xmax


def _nearwell(cur, ref, cv, xmax=50.0):
    """(near-well injection-driven CO2 mass [t], near-well fraction of the
    total injection-driven mass) within x < xmax."""
    mask = _xmask(xmax)
    near = 0.0
    tot = 0.0
    for m in A.CARB:
        if m in cur:
            rm = ref.get(m, np.zeros_like(cur[m]))
            d = (cur[m] - rm) * cv / A.MOLARV[m] * A.CO2_STOICH[m]
            near += float(d[mask, :, :].sum())
            tot += float(d.sum())
    near_t = near * A.M_CO2 / 1e3
    frac = (100.0 * near / tot) if tot > 0 else float("nan")
    return near_t, frac


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", type=float, nargs="+", default=[30.0, 100.0],
                    help="readout years (default: 30 100)")
    args = ap.parse_args()
    years = args.years

    # No-injection control (same locations the module uses).
    ctrl = A.runs_dir("07_da_consistency") / "da_q0_k1"
    if not ctrl.is_dir():
        ctrl = A.runs_dir("09_scco2_kappa_controls") / "sk_ctrl_k1e0"
    ctrl_at = {y: _fields_at(ctrl, y) for y in years}

    print("=" * 100)
    print("  TWO-TIMEPOINT CARBONATION SUMMARY  --  end of injection (t=30) vs end of monitoring (t=100)")
    print("  Injected mass is over the 30-yr injection window: injection stops at t=30, so the")
    print("  denominator is IDENTICAL at both readout years -- only the mineralized numerator grows.")
    print("=" * 100)
    for y in years:
        c = ctrl_at[y]
        if c:
            bg = _co2_t(c["carb"], c["carb0"], c["cv"])
            print(f"  no-injection background @ t={c['t_used']:>5.0f} yr : {bg:6.2f} t CO2 "
                  f"(subtracted per mineral to isolate injection-driven yield)")
    print()

    hdr = (f"  {'scenario':<15}{'t[yr]':>6}{'phi':>9}{'dphi[e-4]':>10}"
           f"{'CarbVF[e-5]':>12}{'inj[t]':>8}{'mineral[t]':>11}{'eta[%]':>8}"
           f"{'nearwell[t]':>12}{'nw[%]':>7}")

    for stem, label in A.BASE_SC:
        run = A.runs_dir("01_baseline") / stem
        # injected CO2 mass is time-independent (injection ends at 30 yr)
        injected_t = (A.well_air_kg(run) + A.liquid_co2_kg(stem)) / 1e3
        print("  " + "-" * 98)
        printed_hdr = False
        for y in years:
            fa = _fields_at(run, y)
            c = ctrl_at[y]
            if fa is None or c is None:
                print(f"  {label:<15}{y:>6.0f}   (missing run or snapshot)")
                continue
            if not printed_hdr:
                print(hdr)
                printed_hdr = True
            phi = float(fa["poro"].mean()) if fa["poro"] is not None else float("nan")
            phi0 = float(fa["poro0"].mean()) if fa["poro0"] is not None else 0.15
            dphi = (phi - phi0) * 1e4
            cvf = _carb_vf_mean(fa["carb"], fa["cv"]) * 1e5
            mineral_t = _co2_t(fa["carb"], c["carb"], fa["cv"])       # injection-driven
            nw_t, nw_frac = _nearwell(fa["carb"], c["carb"], fa["cv"])
            eta = (100.0 * mineral_t / injected_t) if injected_t > 0 else float("nan")
            print(f"  {label:<15}{fa['t_used']:>6.0f}{phi:>9.5f}{dphi:>10.2f}"
                  f"{cvf:>12.2f}{injected_t:>8.1f}{mineral_t:>11.2f}{eta:>8.2f}"
                  f"{nw_t:>12.2f}{nw_frac:>7.1f}")

    print("\n  CHECK: the t=100 mineral[t] column should reproduce the published")
    print("         injection-driven masses (S1~20.8, S2~0.19, WAG~10.4, S6~8.3 t) and")
    print("         S1 CarbVF~8.2e-5. If it does, the t=30 column is trustworthy.")
    print("  Columns: phi=mean porosity; dphi=phi-phi0; CarbVF=domain-mean total carbonate VF")
    print("  (uncontrolled, matches the results table); mineral=control-subtracted injection-")
    print("  driven CO2 mass; eta=mineral/injected; nearwell=injection-driven mass within 50 m")
    print("  of the well and its % of the domain total.")


if __name__ == "__main__":
    main()
