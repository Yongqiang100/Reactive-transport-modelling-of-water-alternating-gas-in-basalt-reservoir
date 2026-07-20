#!/usr/bin/env python3
"""
analyse_dape_sweep.py

Reads results from Suites A and B, computes the diagnostic
Da and Pe values, and produces three diagnostic figures:

  fig_suiteA_kinetic_vs_rate.pdf:
      Side-by-side comparison of the original rate sweep
      vs the kinetic-rate sweep at fixed flow. If the curves
      overlap when plotted against Da, this proves Da governs
      the behaviour (not Pe).

  fig_suiteB_diffusion.pdf:
      Carbonation outcome as a function of D/D₀ at fixed flow.
      If the curve is flat, diffusion is irrelevant across the
      explored Pe range.

  fig_dape_diagnostic.pdf:
      Unified Da-Pe scatter overlay showing where each sweep
      lives in (Da, Pe) space.

Outputs go to ./figures/
"""

import os
import glob
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

try:
    import h5py
    HAS_H5 = True
except ImportError:
    HAS_H5 = False
    print("WARNING: h5py not available; will fall back to JSON if present")

# Paths
BASE = Path(os.environ.get('BASE_DIR', os.path.expanduser('~/WAG/DaPe-disentangling')))
OUT  = Path(os.environ.get('OUT', os.path.expanduser('~/WAG/figures')))
OUT.mkdir(parents=True, exist_ok=True)

# Constants for Da, Pe estimation
D_REF = 1.0e-9       # m²/s reference diffusion
L_CHAR = 100.0       # m characteristic length
Q_BASE = 1.0e-5      # m/s base Darcy velocity at the 20 m well

# Map of run identifiers → (suite, κ, D)
def parse_label(suite, name):
    """
    suite: 'suiteA' or 'suiteB'
    name: subdirectory name, e.g. 'dissolved_kappa1p0' or 'scco2_D10p0'
    Returns (suite_letter, kappa, D).
    """
    if suite == 'suiteA':
        # Format: <scenario>_kappa<value>
        if 'kappa' not in name:
            return ('A', None, None)
        kappa_tag = name.split('kappa')[-1]
        try:
            kappa = float(kappa_tag.replace('p', '.'))
        except ValueError:
            return ('A', None, None)
        return ('A', kappa, D_REF)
    if suite == 'suiteB':
        # Format: <scenario>_D<value>
        # Use rsplit so '_D' matches the suffix-D and not 'dissolved'
        if '_D' not in name:
            return ('B', None, None)
        D_tag = name.rsplit('_D', 1)[-1]
        try:
            D = float(D_tag.replace('p', '.')) * 1.0e-9
        except ValueError:
            return ('B', None, None)
        return ('B', 1.0, D)
    return (None, None, None)

def da_pe(kappa, D, q=Q_BASE, L=L_CHAR):
    """Compute relative Da and Pe; return (None, None) if inputs invalid."""
    if kappa is None or D is None:
        return None, None
    return kappa, q * L / D

# ---------------------------------------------------------------------
# HDF5 reader for carbonate VF at final time
# ---------------------------------------------------------------------
def read_final_carbonate(sim_dir):
    """Read the final-time mean carbonate volume fraction from PFLOTRAN HDF5 output.

    Time-group format: 'Time:  3.00000E+01 y'  (note two spaces after colon)
    Mineral VF dataset format: 'Calcite_VF [m^3 mnrl_m^3 bulk]'
    Sums calcite + magnesite + siderite + dolomite-ord domain-mean VF.
    """
    if not HAS_H5:
        return None
    h5s = sorted(glob.glob(str(sim_dir / "*.h5")))
    h5s = [p for p in h5s if 'hanford' not in os.path.basename(p).lower()]
    if not h5s:
        return None
    try:
        with h5py.File(h5s[-1], 'r') as f:
            # Parse time-group names of form "Time:  3.00000E+01 y"
            def time_value(grp_name):
                try:
                    # Strip 'Time:' prefix and ' y' suffix, parse the number
                    s = grp_name.replace('Time:', '').strip()
                    # Drop trailing unit (' y', ' yr', etc.)
                    s = s.split()[0]
                    return float(s)
                except (ValueError, IndexError):
                    return -1.0

            time_grps = sorted(
                [g for g in f.keys() if g.startswith('Time')],
                key=time_value
            )
            if not time_grps:
                return None

            # Diagnostic: print which group we pick (helpful when debugging)
            final_grp = time_grps[-1]
            g = f[final_grp]

            carb_total = 0.0
            minerals_found = []
            for mineral in ['Calcite', 'Magnesite', 'Siderite', 'Dolomite-ord']:
                # The actual key includes units, e.g. 'Calcite_VF [m^3 mnrl_m^3 bulk]'
                # Find any dataset whose name starts with '<Mineral>_VF'
                matches = [k for k in g.keys() if k.startswith(f'{mineral}_VF')]
                if not matches:
                    continue
                key = matches[0]
                arr = np.array(g[key])
                carb_total += float(arr.mean())
                minerals_found.append(mineral)

            if not minerals_found:
                print(f"  warn: no carbonate VF datasets found in {sim_dir.name} group '{final_grp}'")
                return None
            return carb_total
    except Exception as e:
        print(f"  warn reading {sim_dir}: {e}")
        return None

# ---------------------------------------------------------------------
# Build results table
# ---------------------------------------------------------------------
def collect_results():
    """Walk both suite directories and tabulate the results."""
    results = []
    for suite in ['suiteA', 'suiteB']:
        suite_dir = BASE / suite
        if not suite_dir.exists():
            print(f"  (skipping nonexistent {suite_dir})")
            continue
        for sub in sorted(suite_dir.iterdir()):
            if not sub.is_dir(): continue
            label = sub.name
            # Identify scenario from subdirectory name
            scenario = None
            for scenario_key in ('dissolved', 'scco2'):
                if scenario_key in label:
                    scenario = scenario_key
                    break
            if scenario is None:
                continue
            sweep, kappa, D = parse_label(suite, label)
            if kappa is None or D is None:
                print(f"  (skipping unparseable label: {label})")
                continue
            carb = read_final_carbonate(sub)
            da, pe = da_pe(kappa, D)
            results.append({
                'suite': sweep,
                'scenario': scenario,
                'label': label,
                'kappa': kappa,
                'D': D,
                'Da_rel': da,
                'Pe': pe,
                'carbonate_VF': carb,
            })
    return results

# ---------------------------------------------------------------------
# Generate three diagnostic figures
# ---------------------------------------------------------------------
def fig_suiteA(results):
    """Original rate sweep vs kinetic sweep — both expressed as Da.
    If they overlap, Da governs.

    NOTE: we plot SUITE A only here; the original rate-sweep data must
    be loaded separately if you want the overlay (next function).
    """
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    fig.subplots_adjust(wspace=0.32, bottom=0.20, right=0.96, left=0.10)

    colours = {'dissolved': '#1f77b4', 'scco2': '#d62728'}
    markers = {'dissolved': 'o', 'scco2': 's'}

    suiteA = [r for r in results if r['suite'] == 'A' and r['carbonate_VF'] is not None]

    for ax, x_var, x_label in zip(axes,
                                  ['kappa', 'Da_rel'],
                                  [r'$\kappa$ (kinetic-rate multiplier)',
                                   r'$Da_{\mathrm{rel}}$ ($\propto \kappa$)']):
        for scenario in ['dissolved', 'scco2']:
            pts = [(r[x_var], r['carbonate_VF'])
                   for r in suiteA if r['scenario'] == scenario]
            if not pts:
                continue
            pts.sort()
            xs, ys = zip(*pts)
            ax.plot(xs, ys, marker=markers[scenario], color=colours[scenario],
                    markersize=5, lw=1.2, mec='white', mew=0.5,
                    label='Dissolved (S1)' if scenario=='dissolved' else 'scCO$_2$ (S2)')
        ax.set_xscale('log')
        ax.set_xlabel(x_label, fontsize=8)
        ax.set_ylabel('Carbonate VF at 30 yr', fontsize=8)
        ax.tick_params(labelsize=7)

    axes[0].legend(fontsize=6, loc='best', frameon=False)
    fig.suptitle('Suite A — Kinetic rate sweep at fixed flow (constant Pe)',
                 fontsize=9, y=0.99)

    out_file = OUT / 'fig_suiteA_kinetic.pdf'
    fig.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ {out_file}")

def fig_suiteB(results):
    """Carbonation vs D/D₀ at fixed flow. Flat curve = diffusion irrelevant."""
    fig, ax = plt.subplots(figsize=(4.5, 3.0))
    fig.subplots_adjust(bottom=0.18, left=0.16, right=0.95, top=0.93)

    colours = {'dissolved': '#1f77b4', 'scco2': '#d62728'}
    markers = {'dissolved': 'o', 'scco2': 's'}

    suiteB = [r for r in results if r['suite'] == 'B' and r['carbonate_VF'] is not None]

    for scenario in ['dissolved', 'scco2']:
        pts = [(r['D'] / D_REF, r['carbonate_VF'])
               for r in suiteB if r['scenario'] == scenario]
        if not pts:
            continue
        pts.sort()
        xs, ys = zip(*pts)
        ax.plot(xs, ys, marker=markers[scenario], color=colours[scenario],
                markersize=5, lw=1.2, mec='white', mew=0.5,
                label='Dissolved (S1)' if scenario=='dissolved' else 'scCO$_2$ (S2)')

    ax.set_xscale('log')
    ax.set_xlabel(r'$D/D_0$ (diffusion-coefficient multiplier)', fontsize=8)
    ax.set_ylabel('Carbonate VF at 30 yr', fontsize=8)
    ax.legend(fontsize=6, loc='best', frameon=False)
    ax.tick_params(labelsize=7)
    ax.set_title('Suite B — Diffusion coefficient at fixed flow (constant Da)',
                 fontsize=8.5)

    out_file = OUT / 'fig_suiteB_diffusion.pdf'
    fig.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ {out_file}")

def fig_dape_overlay(results):
    """Locate each run in (Da, Pe) space.

    For Suite A (κ varies, D fixed): each point shifts Da only → vertical
        line in (Da, Pe) at constant Pe.
    For Suite B (D varies, κ fixed): each point shifts Pe only → horizontal
        line in (Da, Pe) at constant Da.
    The original rate sweep would move along a diagonal (both Da and Pe
        change with rate).
    """
    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    fig.subplots_adjust(bottom=0.15, left=0.15)

    colA = {'dissolved': '#1f77b4', 'scco2': '#d62728'}

    for r in results:
        if r['carbonate_VF'] is None:
            continue
        m = 'o' if r['scenario'] == 'dissolved' else 's'
        c = colA[r['scenario']]
        fc = c if r['suite'] == 'A' else 'none'
        ax.scatter(r['Da_rel'], r['Pe'], marker=m, s=50, c=fc if fc!='none' else 'white',
                   edgecolor=c, linewidths=1.2,
                   label=f"Suite {r['suite']}: {r['scenario']}")

    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'$Da_{\mathrm{rel}}$ (kinetic-rate scaling)', fontsize=8)
    ax.set_ylabel(r'$Pe = q L / D$', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title('Da–Pe disentangling experimental design', fontsize=9)

    # Deduplicate legend
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), fontsize=6, loc='best',
              frameon=False, ncol=2)
    ax.axhline(1, color='#cccccc', lw=0.5)
    ax.axvline(1, color='#cccccc', lw=0.5)
    ax.text(0.95, 0.05, 'Pe ≫ 1 region\n(advection-dominated)',
            transform=ax.transAxes, fontsize=6, ha='right', color='#666666')

    out_file = OUT / 'fig_dape_design.pdf'
    fig.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ {out_file}")

# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------
if __name__ == '__main__':
    print("="*60)
    print("  Da-Pe Disentangling Analysis")
    print("="*60)
    print(f"  Base: {BASE}")
    print(f"  Out:  {OUT}")
    print()

    results = collect_results()
    print(f"\nCollected {len(results)} simulation results.")

    # Dump to JSON for downstream
    json_out = OUT / 'dape_results.json'
    with open(json_out, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"  ✓ {json_out}")

    # Print a tidy table
    print(f"\n{'Suite':<6} {'Scenario':<10} {'κ':<8} {'D':<12} {'Da_rel':<10} {'Pe':<12} {'carb_VF':<12}")
    for r in results:
        carb = f"{r['carbonate_VF']:.4e}" if r['carbonate_VF'] is not None else 'N/A'
        print(f"{r['suite']:<6} {r['scenario']:<10} {r['kappa']:<8} "
              f"{r['D']:<12.2e} {r['Da_rel']:<10.3f} {r['Pe']:<12.2e} {carb}")

    print("\nGenerating figures...")
    fig_suiteA(results)
    fig_suiteB(results)
    fig_dape_overlay(results)
    print("\nDone.")
