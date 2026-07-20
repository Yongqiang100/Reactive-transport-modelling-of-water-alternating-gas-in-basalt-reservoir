#!/usr/bin/env python3
"""
Post-process PFLOTRAN WAG results — reads time-stamped HDF5 groups.
Generates comparison plots and summary table.

Run: python3 postprocess_v2.py
"""
import os, glob, json
import numpy as np
import h5py

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

BASE = os.path.dirname(os.path.abspath(__file__))

SCENARIOS = [
    ("scenario1_dissolved", "S1: Dissolved CO₂",  "#2563eb"),
    ("scenario2_scco2",     "S2: scCO₂",          "#dc2626"),
    ("scenario3_wag6mo",    "S3: WAG-6mo",         "#16a34a"),
    ("scenario4_wag3mo",    "S4: WAG-3mo",         "#9333ea"),
    ("scenario5_swag",      "S5: SWAG",            "#ea580c"),
    ("scenario6_adaptive",  "S6: Adaptive",        "#0891b2"),
]

# Carbonate minerals — these are the CO2 storage products
CARBONATES = ['Calcite', 'Magnesite', 'Siderite', 'Dolomite-ord']
# Primary silicates — these dissolve
PRIMARIES = ['Forsterite', 'Anorthite', 'Diopside']
# Secondary non-carbonates
SECONDARIES = ['Kaolinite', 'SiO2(am)']

def read_h5(scenario_dir):
    """Read all time groups from the PFLOTRAN HDF5 file."""
    h5_files = sorted(glob.glob(os.path.join(scenario_dir, "pflotran*.h5")))
    if not h5_files:
        h5_files = sorted(glob.glob(os.path.join(scenario_dir, "*.h5")))
    # Filter out hanford.dat
    h5_files = [f for f in h5_files if 'hanford' not in f]
    if not h5_files:
        return None

    data = {
        'times_yr': [],
        'porosity_mean': [], 'porosity_nearwell': [], 'porosity_farfield': [],
        'perm_mean': [], 'perm_nearwell': [], 'perm_farfield': [],
        'ph_mean': [], 'ph_nearwell': [], 'ph_farfield': [],
        'gas_sat_mean': [],
        'pressure_mean': [],
    }
    # Mineral volume fractions (domain-averaged)
    for m in CARBONATES + PRIMARIES + SECONDARIES:
        data[f'{m}_vf_mean'] = []
        data[f'{m}_vf_nearwell'] = []

    for h5f in h5_files:
        with h5py.File(h5f, 'r') as f:
            # Find all time groups
            time_groups = sorted(
                [g for g in f.keys() if g.startswith('Time')],
                key=lambda g: f[g].attrs.get('Time (s)', [0])[0]
                    if 'Time (s)' in f[g].attrs
                    else float(g.split(':')[1].strip().split()[0])
            )

            for gname in time_groups:
                grp = f[gname]

                # Extract time in years
                if 'Time (s)' in grp.attrs:
                    t_s = grp.attrs['Time (s)']
                    t_yr = float(t_s) / (365.25 * 86400) if np.isscalar(t_s) else float(t_s[0]) / (365.25 * 86400)
                else:
                    # Parse from group name: "Time:  9.00000E+00 y"
                    parts = gname.split(':')[1].strip().split()
                    t_yr = float(parts[0])

                data['times_yr'].append(t_yr)

                # Grid is (50, 1, 25) = (x, y, z)
                # Near-wellbore: first 10 x-cells (0-100m)
                # Far-field: cells 20-50 (200-500m)
                nw = slice(0, 10)
                ff = slice(20, 50)

                # Porosity
                key = 'Porosity'
                if key in grp:
                    arr = grp[key][:]
                    data['porosity_mean'].append(np.mean(arr))
                    data['porosity_nearwell'].append(np.mean(arr[nw, :, :]))
                    data['porosity_farfield'].append(np.mean(arr[ff, :, :]))

                # Permeability
                key = 'Permeability_X [m^2]'
                if key not in grp:
                    key = 'Permeability_X'
                if key in grp:
                    arr = grp[key][:]
                    data['perm_mean'].append(np.mean(arr))
                    data['perm_nearwell'].append(np.mean(arr[nw, :, :]))
                    data['perm_farfield'].append(np.mean(arr[ff, :, :]))

                # pH
                if 'pH' in grp or 'ph' in grp:
                    key = 'pH' if 'pH' in grp else 'ph'
                    arr = grp[key][:]
                    data['ph_mean'].append(np.mean(arr))
                    data['ph_nearwell'].append(np.mean(arr[nw, :, :]))
                    data['ph_farfield'].append(np.mean(arr[ff, :, :]))

                # Gas saturation
                for key in ['Gas_Saturation', 'Gas Saturation']:
                    if key in grp:
                        data['gas_sat_mean'].append(np.mean(grp[key][:]))
                        break

                # Pressure
                for key in ['Liquid_Pressure [Pa]', 'Liquid_Pressure']:
                    if key in grp:
                        data['pressure_mean'].append(np.mean(grp[key][:]))
                        break

                # Mineral volume fractions
                for m in CARBONATES + PRIMARIES + SECONDARIES:
                    for suffix in [f'{m}_VF [m^3 mnrl_m^3 bulk]', f'{m}_VF']:
                        if suffix in grp:
                            arr = grp[suffix][:]
                            data[f'{m}_vf_mean'].append(np.mean(arr))
                            data[f'{m}_vf_nearwell'].append(np.mean(arr[nw, :, :]))
                            break
                    else:
                        data[f'{m}_vf_mean'].append(np.nan)
                        data[f'{m}_vf_nearwell'].append(np.nan)

    # Convert to numpy
    for k in data:
        data[k] = np.array(data[k])

    return data


def read_mass_balance(scenario_dir):
    """Read mass balance .dat file."""
    mas_files = glob.glob(os.path.join(scenario_dir, "*-mas*.dat"))
    if not mas_files:
        return None
    try:
        with open(mas_files[0]) as fh:
            header = fh.readline().strip().replace('"', '').split(',')
        raw = np.genfromtxt(mas_files[0], skip_header=1, delimiter=',')
        return {'columns': [h.strip() for h in header], 'data': raw}
    except:
        try:
            with open(mas_files[0]) as fh:
                header = fh.readline().strip().replace('"', '').split()
            raw = np.loadtxt(mas_files[0], skiprows=1)
            return {'columns': header, 'data': raw}
        except Exception as e:
            print(f"  Warning: {e}")
            return None


# ============================================================
# Read all scenarios
# ============================================================
print("=" * 65)
print("  PFLOTRAN WAG Results — Post-Processing")
print("=" * 65)

all_data = {}
all_mass = {}

for dirname, label, color in SCENARIOS:
    sdir = os.path.join(BASE, dirname)
    if not os.path.isdir(sdir):
        print(f"  ✘ {dirname}/ not found")
        continue

    print(f"\n  Reading {label}...")
    d = read_h5(sdir)
    if d and len(d['times_yr']) > 0:
        all_data[dirname] = d
        print(f"    {len(d['times_yr'])} time snapshots, "
              f"t = {d['times_yr'][0]:.1f} – {d['times_yr'][-1]:.1f} yr")
        if len(d['porosity_mean']) > 0:
            print(f"    Porosity: {d['porosity_mean'][0]:.4f} → {d['porosity_mean'][-1]:.4f}")
        if len(d['ph_mean']) > 0:
            print(f"    pH: {d['ph_mean'][0]:.2f} → {d['ph_mean'][-1]:.2f}")
        # Total carbonate volume fraction
        carb_total = sum(d[f'{m}_vf_mean'] for m in CARBONATES)
        if len(carb_total) > 0 and not np.all(np.isnan(carb_total)):
            print(f"    Total carbonate VF: {np.nanmean(carb_total[0]):.5f} → {np.nanmean(carb_total[-1]):.5f}")
    else:
        print(f"    No time-series data found")

    m = read_mass_balance(sdir)
    if m:
        all_mass[dirname] = m

# ============================================================
# Summary table
# ============================================================
print("\n" + "=" * 65)
print("  SCENARIO COMPARISON (Final State)")
print("=" * 65)
print(f"{'Scenario':<22} {'φ_final':>8} {'k/k₀':>8} {'pH':>6} "
      f"{'Carb VF':>9} {'Calcite':>9} {'Mgsite':>9}")
print("-" * 75)

summary = {}
for dirname, label, color in SCENARIOS:
    if dirname not in all_data:
        continue
    d = all_data[dirname]
    row = {'label': label}

    phi_f = d['porosity_mean'][-1] if len(d['porosity_mean']) else np.nan
    phi_0 = d['porosity_mean'][0] if len(d['porosity_mean']) else 0.15
    k_f = d['perm_mean'][-1] if len(d['perm_mean']) else np.nan
    k_0 = d['perm_mean'][0] if len(d['perm_mean']) else 1e-13
    k_ratio = k_f / k_0 if k_0 > 0 else np.nan
    ph_f = d['ph_mean'][-1] if len(d['ph_mean']) else np.nan

    carb_vf = sum(d[f'{m}_vf_mean'][-1] for m in CARBONATES
                  if len(d[f'{m}_vf_mean']) > 0 and not np.isnan(d[f'{m}_vf_mean'][-1]))
    calc_vf = d['Calcite_vf_mean'][-1] if len(d['Calcite_vf_mean']) > 0 else np.nan
    mgs_vf = d['Magnesite_vf_mean'][-1] if len(d['Magnesite_vf_mean']) > 0 else np.nan

    row.update({
        'porosity_final': float(phi_f),
        'porosity_initial': float(phi_0),
        'perm_ratio': float(k_ratio),
        'ph_final': float(ph_f),
        'carbonate_vf': float(carb_vf),
        'calcite_vf': float(calc_vf),
        'magnesite_vf': float(mgs_vf),
    })
    summary[dirname] = row

    print(f"{label:<22} {phi_f:>8.5f} {k_ratio:>8.4f} {ph_f:>6.2f} "
          f"{carb_vf:>9.6f} {calc_vf:>9.6f} {mgs_vf:>9.6f}")

# ============================================================
# Plots
# ============================================================
print("\nGenerating figures...")

fig = plt.figure(figsize=(18, 14))
gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)
fig.suptitle('WAG CO₂ Mineralization in Basalt — PFLOTRAN Results', fontsize=15, fontweight='bold', y=0.98)

# 1. Total carbonate VF vs time
ax1 = fig.add_subplot(gs[0, 0])
for dirname, label, color in SCENARIOS:
    if dirname not in all_data: continue
    d = all_data[dirname]
    carb = sum(d[f'{m}_vf_mean'] for m in CARBONATES)
    if len(carb) == len(d['times_yr']):
        ax1.plot(d['times_yr'], carb, color=color, label=label, lw=2)
ax1.set_xlabel('Time (yr)')
ax1.set_ylabel('Total carbonate VF')
ax1.set_title('Carbonate Precipitation')
ax1.legend(fontsize=7)
ax1.grid(True, alpha=0.3)

# 2. Porosity vs time
ax2 = fig.add_subplot(gs[0, 1])
for dirname, label, color in SCENARIOS:
    if dirname not in all_data: continue
    d = all_data[dirname]
    if len(d['porosity_mean']) == len(d['times_yr']):
        ax2.plot(d['times_yr'], d['porosity_mean'], color=color, label=label, lw=2)
ax2.set_xlabel('Time (yr)')
ax2.set_ylabel('Porosity (mean)')
ax2.set_title('Porosity Evolution')
ax2.legend(fontsize=7)
ax2.grid(True, alpha=0.3)

# 3. Permeability ratio vs time
ax3 = fig.add_subplot(gs[0, 2])
for dirname, label, color in SCENARIOS:
    if dirname not in all_data: continue
    d = all_data[dirname]
    if len(d['perm_mean']) == len(d['times_yr']) and d['perm_mean'][0] > 0:
        ax3.plot(d['times_yr'], d['perm_mean'] / d['perm_mean'][0],
                 color=color, label=label, lw=2)
ax3.set_xlabel('Time (yr)')
ax3.set_ylabel('k / k₀')
ax3.set_title('Permeability Ratio')
ax3.legend(fontsize=7)
ax3.grid(True, alpha=0.3)
ax3.set_yscale('log')

# 4. pH vs time
ax4 = fig.add_subplot(gs[1, 0])
for dirname, label, color in SCENARIOS:
    if dirname not in all_data: continue
    d = all_data[dirname]
    if len(d['ph_mean']) == len(d['times_yr']):
        ax4.plot(d['times_yr'], d['ph_mean'], color=color, label=label, lw=2)
ax4.set_xlabel('Time (yr)')
ax4.set_ylabel('pH')
ax4.set_title('pH Evolution')
ax4.legend(fontsize=7)
ax4.grid(True, alpha=0.3)

# 5. Individual carbonate minerals (stacked for S1)
ax5 = fig.add_subplot(gs[1, 1])
if 'scenario1_dissolved' in all_data:
    d = all_data['scenario1_dissolved']
    t = d['times_yr']
    bottom = np.zeros(len(t))
    mineral_colors = {'Calcite': '#2563eb', 'Magnesite': '#16a34a',
                      'Siderite': '#ea580c', 'Dolomite-ord': '#9333ea'}
    for m in CARBONATES:
        vf = d[f'{m}_vf_mean']
        if len(vf) == len(t) and not np.all(np.isnan(vf)):
            vf_clean = np.nan_to_num(vf)
            ax5.fill_between(t, bottom, bottom + vf_clean,
                           alpha=0.7, label=m, color=mineral_colors.get(m, '#999'))
            bottom += vf_clean
ax5.set_xlabel('Time (yr)')
ax5.set_ylabel('Volume Fraction')
ax5.set_title('Carbonate Minerals (S1: Dissolved)')
ax5.legend(fontsize=7)
ax5.grid(True, alpha=0.3)

# 6. Primary mineral dissolution
ax6 = fig.add_subplot(gs[1, 2])
for dirname, label, color in SCENARIOS:
    if dirname not in all_data: continue
    d = all_data[dirname]
    forst = d['Forsterite_vf_mean']
    if len(forst) == len(d['times_yr']) and not np.all(np.isnan(forst)):
        ax6.plot(d['times_yr'], forst, color=color, label=label, lw=2)
ax6.set_xlabel('Time (yr)')
ax6.set_ylabel('Forsterite VF')
ax6.set_title('Olivine Dissolution')
ax6.legend(fontsize=7)
ax6.grid(True, alpha=0.3)

# 7. Near-wellbore vs far-field porosity (S1)
ax7 = fig.add_subplot(gs[2, 0])
if 'scenario1_dissolved' in all_data:
    d = all_data['scenario1_dissolved']
    t = d['times_yr']
    if len(d['porosity_nearwell']) == len(t):
        ax7.plot(t, d['porosity_nearwell'], 'b-', lw=2, label='Near-wellbore')
    if len(d['porosity_farfield']) == len(t):
        ax7.plot(t, d['porosity_farfield'], 'r--', lw=2, label='Far-field')
ax7.set_xlabel('Time (yr)')
ax7.set_ylabel('Porosity')
ax7.set_title('Spatial Porosity (S1: Dissolved)')
ax7.legend(fontsize=8)
ax7.grid(True, alpha=0.3)

# 8. Gas saturation
ax8 = fig.add_subplot(gs[2, 1])
for dirname, label, color in SCENARIOS:
    if dirname not in all_data: continue
    d = all_data[dirname]
    if len(d['gas_sat_mean']) == len(d['times_yr']):
        ax8.plot(d['times_yr'], d['gas_sat_mean'], color=color, label=label, lw=2)
ax8.set_xlabel('Time (yr)')
ax8.set_ylabel('Gas Saturation')
ax8.set_title('Free CO₂ Phase')
ax8.legend(fontsize=7)
ax8.grid(True, alpha=0.3)

# 9. Final carbonate VF bar chart
ax9 = fig.add_subplot(gs[2, 2])
labels_plot = []
carb_finals = []
colors_plot = []
for dirname, label, color in SCENARIOS:
    if dirname not in all_data: continue
    d = all_data[dirname]
    carb = sum(d[f'{m}_vf_mean'][-1] for m in CARBONATES
               if len(d[f'{m}_vf_mean']) > 0 and not np.isnan(d[f'{m}_vf_mean'][-1]))
    labels_plot.append(label)
    carb_finals.append(carb)
    colors_plot.append(color)
bars = ax9.bar(labels_plot, carb_finals, color=colors_plot, alpha=0.85)
ax9.set_ylabel('Total Carbonate VF')
ax9.set_title('Final Carbonate (30 yr)')
ax9.tick_params(axis='x', rotation=30)
for bar, val in zip(bars, carb_finals):
    ax9.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
             f'{val:.5f}', ha='center', va='bottom', fontsize=8)
ax9.grid(True, alpha=0.3, axis='y')

plot_path = os.path.join(BASE, "wag_results_comparison.png")
plt.savefig(plot_path, dpi=150, bbox_inches='tight')
print(f"  ✔ Saved {plot_path}")

# ============================================================
# Export JSON
# ============================================================
out_json = os.path.join(BASE, "pflotran_results.json")
export = {}
for dirname, label, color in SCENARIOS:
    if dirname not in all_data: continue
    d = all_data[dirname]
    export[dirname] = {
        'label': label,
        'times_yr': d['times_yr'].tolist(),
        'porosity_mean': d['porosity_mean'].tolist(),
        'perm_mean': d['perm_mean'].tolist(),
        'ph_mean': d['ph_mean'].tolist(),
    }
    for m in CARBONATES:
        export[dirname][f'{m}_vf'] = d[f'{m}_vf_mean'].tolist()

with open(out_json, 'w') as f:
    json.dump(export, f, indent=2)
print(f"  ✔ Saved {out_json}")

print("\n" + "=" * 65)
print("  Done! Check wag_results_comparison.png")
print("=" * 65)
