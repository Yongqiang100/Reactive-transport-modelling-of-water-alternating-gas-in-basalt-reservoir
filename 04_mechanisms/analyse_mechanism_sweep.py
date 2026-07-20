#!/usr/bin/env python3
"""
analyse_mechanism_sweep.py

Reads Case C/D/E mechanism sweep results and produces three diagnostic
figures plus a unified mechanism summary figure.

For each case, extracts:
  - Final-time mean carbonate volume fraction (the primary outcome)
  - Gas-water contact metric (Σ|∇Sg| in the near-wellbore region)
  - Buoyancy override fraction (gas mass above mid-height)
  - Mean gas saturation

Produces:
  fig_caseC_phase.pdf       — carbonate vs mole fraction
  fig_caseD_buoyancy.pdf    — carbonate vs injection elevation
  fig_caseE_contact.pdf     — carbonate vs gas residual saturation
  fig_mechanisms_summary.pdf — 3-panel composite figure
  mechanisms_results.json   — raw data

Plus three cross-section figures showing gas saturation maps for each case.

Usage:
    BASE_DIR=~/WAG OUT=~/WAG/figures python3 analyse_mechanism_sweep.py
"""

import os
import glob
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

try:
    import h5py
except ImportError:
    print("ERROR: h5py required")
    raise SystemExit(1)

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 8

BASE = Path(os.environ.get('BASE_DIR', os.path.expanduser('~/WAG')))
OUT  = Path(os.environ.get('OUT', os.path.expanduser('~/WAG/figures')))
OUT.mkdir(parents=True, exist_ok=True)

# Subdirectory → parameter value map for each case
CASE_C_MAP = {  # mole fraction
    'phase_xfrac0p20': 0.20,
    'phase_xfrac0p40': 0.40,
    'phase_xfrac0p60': 0.60,
    'phase_xfrac0p80': 0.80,
}
CASE_D_MAP = {  # injection z-position (center of well)
    'buoyancy_top':    15.0,    # z=0-30, center=15
    'buoyancy_bottom': 85.0,    # z=70-100, center=85
}
CASE_E_MAP = {  # gas residual saturation
    'contact_Sgres0p05': 0.05,
    'contact_Sgres0p20': 0.20,
    'contact_Sgres0p30': 0.30,
}

# Baselines from the Da-Pe Case A scco2_kappa1p0 run already on disk
# These will be auto-detected if available
BASELINE_VALUES = {
    'mole_fraction_baseline': (0.99, None),   # to be filled from existing S2 run
    'elevation_baseline': (50.0, None),       # mid-formation (z=20-80, center=50)
    'sgres_baseline': (0.10, None),
}


# ---------------------------------------------------------------------
# HDF5 reader
# ---------------------------------------------------------------------
def read_final_state(sim_dir):
    h5s = sorted(glob.glob(str(sim_dir / "*.h5")))
    h5s = [p for p in h5s if 'hanford' not in os.path.basename(p).lower()]
    if not h5s:
        return None
    try:
        with h5py.File(h5s[-1], 'r') as f:
            def time_value(name):
                try:
                    s = name.replace('Time:', '').strip().split()[0]
                    return float(s)
                except (ValueError, IndexError):
                    return -1.0
            time_grps = sorted([g for g in f.keys() if g.startswith('Time')],
                               key=time_value)
            if not time_grps:
                return None
            g = f[time_grps[-1]]

            # Gas saturation (Nx, Ny=1, Nz) → 2D
            if 'Gas_Saturation' not in g:
                return None
            gs = np.array(g['Gas_Saturation'])
            gs_2d = gs[:, 0, :]

            # Carbonate VF
            carb_total = 0.0
            for mineral in ['Calcite', 'Magnesite', 'Siderite', 'Dolomite-ord']:
                matches = [k for k in g.keys() if k.startswith(f'{mineral}_VF')]
                if matches:
                    carb_total += float(np.array(g[matches[0]]).mean())

            return {
                'gas_sat_2d': gs_2d,
                'carb_total': carb_total,
            }
    except Exception as e:
        print(f"  warn reading {sim_dir}: {e}")
        return None


def gas_water_contact(gas_sat_2d, x_max_idx=95):
    """Σ|∇Sg| over the near-wellbore region (first 95 cells, ~0-500m)."""
    region = gas_sat_2d[:x_max_idx, :]
    gx, gz = np.gradient(region)
    return float(np.sqrt(gx**2 + gz**2).sum())


def buoyancy_override(gas_sat_2d):
    """Fraction of gas mass above the formation mid-height."""
    upper = gas_sat_2d[:, :12]   # z indices 0-11 = upper half
    lower = gas_sat_2d[:, 13:]   # z indices 13-24 = lower half
    u = float(upper.sum())
    l = float(lower.sum())
    return u / (u + l + 1e-20) if (u + l) > 1e-20 else 0.0


def mean_gas_sat(gas_sat_2d, x_max_idx=95):
    return float(gas_sat_2d[:x_max_idx, :].mean())


# ---------------------------------------------------------------------
# Collect results
# ---------------------------------------------------------------------
def collect():
    mech_dir = BASE / 'mechanisms'
    if not mech_dir.exists():
        print(f"ERROR: {mech_dir} does not exist")
        return None

    results = {'C': [], 'D': [], 'E': []}

    for case_letter, case_map, key_name in [
        ('C', CASE_C_MAP, 'mole_fraction'),
        ('D', CASE_D_MAP, 'z_center'),
        ('E', CASE_E_MAP, 'gas_residual_sat'),
    ]:
        case_root = mech_dir / f'case{case_letter}'
        if not case_root.exists():
            print(f"  (no {case_root})")
            continue
        for subname, value in case_map.items():
            sub = case_root / subname
            if not sub.is_dir():
                continue
            state = read_final_state(sub)
            if state is None:
                print(f"  (no usable HDF5 in {sub.name})")
                continue
            gs2d = state['gas_sat_2d']
            results[case_letter].append({
                key_name: value,
                'label': subname,
                'carb_total': state['carb_total'],
                'contact': gas_water_contact(gs2d),
                'override': buoyancy_override(gs2d),
                'mean_gs': mean_gas_sat(gs2d),
                '_gs2d': gs2d,
            })

    # Also load baseline (Case A scco2_kappa1p0) if available
    baseline_sub = BASE / 'suiteA' / 'scco2_kappa1p0'
    if baseline_sub.exists():
        state = read_final_state(baseline_sub)
        if state is not None:
            BASELINE_VALUES['mole_fraction_baseline'] = (0.99, state['carb_total'])
            BASELINE_VALUES['elevation_baseline'] = (50.0, state['carb_total'])
            BASELINE_VALUES['sgres_baseline'] = (0.10, state['carb_total'])
            BASELINE_VALUES['_baseline_gs2d'] = state['gas_sat_2d']

    # Also load S1 baseline (dissolved) for Case C reference
    s1_baseline = BASE / 'suiteA' / 'dissolved_kappa1p0'
    if s1_baseline.exists():
        state = read_final_state(s1_baseline)
        if state is not None:
            BASELINE_VALUES['s1_dissolved'] = (0.04, state['carb_total'])

    return results


# ---------------------------------------------------------------------
# Individual case figures
# ---------------------------------------------------------------------
def fig_caseC(results):
    """Phase partitioning: carbonate vs mole fraction."""
    pts = sorted(results['C'], key=lambda r: r['mole_fraction'])
    if not pts:
        return
    fig, ax = plt.subplots(figsize=(4.5, 3.0))
    fig.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.18)

    # Add S1 (dissolved) and S2 (scco2) endpoints if available
    extra_x, extra_y = [], []
    if BASELINE_VALUES.get('s1_dissolved'):
        x, y = BASELINE_VALUES['s1_dissolved']
        if y is not None:
            extra_x.append(x); extra_y.append(y)
    s2x, s2y = BASELINE_VALUES.get('mole_fraction_baseline', (None, None))
    if s2y is not None:
        extra_x.append(s2x); extra_y.append(s2y)

    x_swept = [r['mole_fraction'] for r in pts]
    y_swept = [r['carb_total'] for r in pts]
    all_x = extra_x[:1] + x_swept + (extra_x[1:] if len(extra_x) > 1 else [])
    all_y = extra_y[:1] + y_swept + (extra_y[1:] if len(extra_y) > 1 else [])

    # Plot the continuum
    ax.plot(all_x, np.array(all_y)*1e4, marker='o', color='#9c27b0',
            markersize=6, lw=1.5, mec='white', mew=0.6, zorder=3)

    # Highlight endpoints if present
    if BASELINE_VALUES.get('s1_dissolved') and BASELINE_VALUES['s1_dissolved'][1] is not None:
        ax.scatter([0.04], [BASELINE_VALUES['s1_dissolved'][1]*1e4],
                   s=90, marker='o', color='#1f6fb4', edgecolor='black',
                   lw=0.6, zorder=4, label='S1 dissolved (baseline)')
    if s2y is not None:
        ax.scatter([s2x], [s2y*1e4], s=90, marker='s', color='#c92e1d',
                   edgecolor='black', lw=0.6, zorder=4,
                   label='S2 scCO$_2$ (baseline)')

    ax.set_xlabel(r'Injection CO$_2$ mole fraction $x_{CO_2}$', fontsize=9)
    ax.set_ylabel(r'Carbonate VF at 30 yr ($\times 10^{-4}$)', fontsize=9)
    ax.set_title('Case C: phase partitioning', fontsize=9)
    ax.set_xlim(-0.05, 1.05)
    ax.legend(loc='lower left', fontsize=7, frameon=False)
    ax.tick_params(labelsize=8)
    plt.savefig(OUT / 'fig_caseC_phase.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f"  ✓ {OUT / 'fig_caseC_phase.pdf'}")


def fig_caseD(results):
    """Buoyancy override: carbonate vs injection elevation."""
    pts = list(results['D'])
    if not pts:
        return
    # Build a 3-point dataset including baseline
    baseline_y = BASELINE_VALUES.get('elevation_baseline', (None, None))[1]
    if baseline_y is not None:
        pts.append({
            'z_center': 50.0,
            'label': 'baseline (middle)',
            'carb_total': baseline_y,
        })
    pts = sorted(pts, key=lambda r: r['z_center'])

    fig, ax = plt.subplots(figsize=(4.5, 3.0))
    fig.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.18)
    z = [r['z_center'] for r in pts]
    y = [r['carb_total'] for r in pts]
    ax.plot(z, np.array(y)*1e4, marker='^', color='#2e7d32',
            markersize=7, lw=1.5, mec='white', mew=0.6)

    # Annotate
    for r in pts:
        ax.annotate(r['label'].replace('buoyancy_', '').replace('_', ' '),
                    xy=(r['z_center'], r['carb_total']*1e4),
                    xytext=(8, 4), textcoords='offset points', fontsize=7)

    ax.set_xlabel('Injection-well center depth (m)', fontsize=9)
    ax.set_ylabel(r'Carbonate VF at 30 yr ($\times 10^{-4}$)', fontsize=9)
    ax.set_title('Case D: buoyancy override (S2 scCO$_2$)', fontsize=9)
    ax.tick_params(labelsize=8)
    plt.savefig(OUT / 'fig_caseD_buoyancy.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f"  ✓ {OUT / 'fig_caseD_buoyancy.pdf'}")


def fig_caseE(results):
    """Gas-water contact: carbonate vs gas residual saturation."""
    pts = list(results['E'])
    if not pts:
        return
    baseline_y = BASELINE_VALUES.get('sgres_baseline', (None, None))[1]
    if baseline_y is not None:
        pts.append({
            'gas_residual_sat': 0.10,
            'label': 'baseline',
            'carb_total': baseline_y,
        })
    pts = sorted(pts, key=lambda r: r['gas_residual_sat'])

    fig, ax = plt.subplots(figsize=(4.5, 3.0))
    fig.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.18)
    x = [r['gas_residual_sat'] for r in pts]
    y = [r['carb_total'] for r in pts]
    ax.plot(x, np.array(y)*1e4, marker='D', color='#f57c00',
            markersize=6, lw=1.5, mec='white', mew=0.6)
    ax.set_xlabel(r'Gas residual saturation $S_{gr}$', fontsize=9)
    ax.set_ylabel(r'Carbonate VF at 30 yr ($\times 10^{-4}$)', fontsize=9)
    ax.set_title('Case E: gas-water contact (S2 scCO$_2$)', fontsize=9)
    ax.tick_params(labelsize=8)
    plt.savefig(OUT / 'fig_caseE_contact.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f"  ✓ {OUT / 'fig_caseE_contact.pdf'}")


# ---------------------------------------------------------------------
# Summary panel
# ---------------------------------------------------------------------
def fig_summary(results):
    """Three-panel composite: phase, buoyancy, contact."""
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 2.8))
    fig.subplots_adjust(left=0.07, right=0.97, top=0.88, bottom=0.20, wspace=0.35)

    # Panel a: Case C
    ax = axes[0]
    pts = sorted(results['C'], key=lambda r: r['mole_fraction'])
    if BASELINE_VALUES.get('s1_dissolved') and BASELINE_VALUES['s1_dissolved'][1] is not None:
        s1 = BASELINE_VALUES['s1_dissolved']
    else:
        s1 = (0.04, None)
    s2 = BASELINE_VALUES.get('mole_fraction_baseline', (0.99, None))

    all_x = ([s1[0]] if s1[1] is not None else []) + [r['mole_fraction'] for r in pts] + ([s2[0]] if s2[1] is not None else [])
    all_y = ([s1[1]] if s1[1] is not None else []) + [r['carb_total'] for r in pts] + ([s2[1]] if s2[1] is not None else [])
    if all_x:
        ax.plot(all_x, np.array(all_y)*1e4, marker='o', color='#9c27b0',
                markersize=5, lw=1.2, mec='white', mew=0.5)
    ax.set_xlabel(r'CO$_2$ mole fraction', fontsize=8)
    ax.set_ylabel(r'Carbonate VF ($\times 10^{-4}$)', fontsize=8)
    ax.set_title('(a) Phase partitioning', fontsize=8.5)
    ax.tick_params(labelsize=7)

    # Panel b: Case D
    ax = axes[1]
    pts = list(results['D'])
    if BASELINE_VALUES.get('elevation_baseline', (None, None))[1] is not None:
        pts.append({'z_center': 50.0, 'carb_total':
                    BASELINE_VALUES['elevation_baseline'][1]})
    pts = sorted(pts, key=lambda r: r['z_center'])
    if pts:
        ax.plot([r['z_center'] for r in pts],
                np.array([r['carb_total'] for r in pts])*1e4,
                marker='^', color='#2e7d32', markersize=6, lw=1.2,
                mec='white', mew=0.5)
    ax.set_xlabel('Well center depth (m)', fontsize=8)
    ax.set_ylabel(r'Carbonate VF ($\times 10^{-4}$)', fontsize=8)
    ax.set_title('(b) Buoyancy override', fontsize=8.5)
    ax.tick_params(labelsize=7)

    # Panel c: Case E
    ax = axes[2]
    pts = list(results['E'])
    if BASELINE_VALUES.get('sgres_baseline', (None, None))[1] is not None:
        pts.append({'gas_residual_sat': 0.10,
                    'carb_total': BASELINE_VALUES['sgres_baseline'][1]})
    pts = sorted(pts, key=lambda r: r['gas_residual_sat'])
    if pts:
        ax.plot([r['gas_residual_sat'] for r in pts],
                np.array([r['carb_total'] for r in pts])*1e4,
                marker='D', color='#f57c00', markersize=5, lw=1.2,
                mec='white', mew=0.5)
    ax.set_xlabel(r'Gas residual sat. $S_{gr}$', fontsize=8)
    ax.set_ylabel(r'Carbonate VF ($\times 10^{-4}$)', fontsize=8)
    ax.set_title('(c) Gas-water contact', fontsize=8.5)
    ax.tick_params(labelsize=7)

    fig.suptitle('Three fundamental mechanisms controlling scCO$_2$ carbonation',
                 fontsize=9, y=0.99)
    plt.savefig(OUT / 'fig_mechanisms_summary.pdf', dpi=300,
                bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f"  ✓ {OUT / 'fig_mechanisms_summary.pdf'}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
if __name__ == '__main__':
    print("="*60)
    print("  Mechanism Sweep Analysis (Cases C, D)")
    print("="*60)
    print(f"  Base: {BASE}")
    print(f"  Out:  {OUT}")
    print()

    results = collect()
    if results is None:
        raise SystemExit(1)

    # Print summary table
    print("\n--- Case C: Phase Partitioning ---")
    print(f"{'mole_frac':<12} {'carb_VF':<14} {'contact':<10} "
          f"{'override':<10} {'mean_Sg':<10}")
    for r in sorted(results['C'], key=lambda x: x['mole_fraction']):
        print(f"{r['mole_fraction']:<12.2f} {r['carb_total']:<14.4e} "
              f"{r['contact']:<10.1f} {r['override']:<10.4f} "
              f"{r['mean_gs']:<10.4f}")

    print("\n--- Case D: Buoyancy Override ---")
    print(f"{'z_center':<12} {'carb_VF':<14} {'contact':<10} "
          f"{'override':<10} {'mean_Sg':<10}")
    for r in sorted(results['D'], key=lambda x: x['z_center']):
        print(f"{r['z_center']:<12.1f} {r['carb_total']:<14.4e} "
              f"{r['contact']:<10.1f} {r['override']:<10.4f} "
              f"{r['mean_gs']:<10.4f}")

    # Case E (gas-residual-saturation sweep) was dropped from the study;
    # only print/plot it if such runs happen to be present.
    if results.get('E'):
        print("\n--- Case E: Gas-Water Contact ---")
        print(f"{'Sgres':<12} {'carb_VF':<14} {'contact':<10} "
              f"{'override':<10} {'mean_Sg':<10}")
        for r in sorted(results['E'], key=lambda x: x['gas_residual_sat']):
            print(f"{r['gas_residual_sat']:<12.2f} {r['carb_total']:<14.4e} "
                  f"{r['contact']:<10.1f} {r['override']:<10.4f} "
                  f"{r['mean_gs']:<10.4f}")

    # Save JSON (strip the gas-saturation arrays)
    json_safe = {}
    for k, v in results.items():
        json_safe[k] = [{kk: vv for kk, vv in r.items() if kk != '_gs2d'}
                        for r in v]
    json_safe['baselines'] = {k: v for k, v in BASELINE_VALUES.items()
                              if not k.startswith('_')}
    json_out = OUT / 'mechanisms_results.json'
    with open(json_out, 'w') as f:
        json.dump(json_safe, f, indent=2, default=str)
    print(f"\n  ✓ {json_out}")

    print("\nGenerating figures...")
    fig_caseC(results)
    fig_caseD(results)
    if results.get('E'):
        fig_caseE(results)
        fig_summary(results)
    print("\nDone.")
