#!/usr/bin/env python3
"""
plot_spatial_fields.py

Direct visual comparison of the spatial fields for Cases D and E.

For each case it builds a 3-row x N-column panel:
    Row 1: Gas saturation     (the CAUSE that varies across runs)
    Row 2: Net carbonate VF   (above the 4e-4 seed)  -- the EFFECT
    Row 3: pH                 (the geochemical state)

Each row shares a color scale across all columns, and the colorbar is
annotated with the actual data range. This makes the comparison honest:
  - the gas row visibly changes across columns (gas redistributes)
  - the carbonate and pH rows look uniform, and the range annotation shows
    the variation is sub-percent / sub-0.1-pH-unit (negligible).

Region shown: near-well, x < 300 m.

Produces fig_fields_caseD.png/.pdf and fig_fields_caseE.png/.pdf.

Usage:
    BASE_DIR=~/WAG OUT=~/WAG/figures python3 plot_spatial_fields.py
"""

import os
import glob
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

DX = np.array([2.0]*50 + [10.0]*40 + [30.0]*50)
DZ = np.array([4.0]*25)
X_EDGES = np.concatenate([[0], np.cumsum(DX)])
Z_EDGES = np.concatenate([[0], np.cumsum(DZ)])
X_CENTERS = np.cumsum(DX) - DX/2.0

CARB_MINERALS = ['Calcite', 'Magnesite', 'Siderite', 'Dolomite-ord']
TOTAL_SEED = 4.0e-4
X_VIEW_MAX = 300.0

CASE_D_RUNS = [
    ('z=15 m\n(well 0-30)',   'mechanisms/caseD/buoyancy_top'),
    ('z=50 m\n(baseline)',    'suiteA/scco2_kappa1p0'),
    ('z=85 m\n(well 70-100)', 'mechanisms/caseD/buoyancy_bottom'),
]
CASE_E_RUNS = [
    ('Sgr=0.02',             'mechanisms/caseE/contact_Sgres0p02'),
    ('Sgr=0.05',             'mechanisms/caseE/contact_Sgres0p05'),
    ('Sgr=0.10\n(baseline)', 'suiteA/scco2_kappa1p0'),
    ('Sgr=0.20',             'mechanisms/caseE/contact_Sgres0p20'),
    ('Sgr=0.30',             'mechanisms/caseE/contact_Sgres0p30'),
    ('Sgr=0.40',             'mechanisms/caseE/contact_Sgres0p40'),
    ('Sgr=0.50',             'mechanisms/caseE/contact_Sgres0p50'),
]


def read_fields(rel):
    """Return (gas_2d, net_carb_2d, ph_2d) at final time, shape (Nx, Nz)."""
    sim_dir = BASE / rel
    h5s = sorted(glob.glob(str(sim_dir / "*.h5")))
    h5s = [p for p in h5s if 'hanford' not in os.path.basename(p).lower()]
    if not h5s:
        print(f"  WARN: no HDF5 in {sim_dir}")
        return None, None, None
    with h5py.File(h5s[-1], 'r') as f:
        def tv(name):
            try:
                return float(name.replace('Time:', '').strip().split()[0])
            except (ValueError, IndexError):
                return -1.0
        tgrps = sorted([g for g in f.keys() if g.startswith('Time')], key=tv)
        g = f[tgrps[-1]]
        gas = np.array(g['Gas_Saturation'])[:, 0, :] if 'Gas_Saturation' in g else None
        carb = None
        for m in CARB_MINERALS:
            mk = [k for k in g.keys() if k.startswith(f'{m}_VF')]
            if mk:
                arr = np.array(g[mk[0]])[:, 0, :]
                carb = arr if carb is None else carb + arr
        if carb is not None:
            carb = carb - TOTAL_SEED
        ph = None
        pk = [k for k in g.keys() if k.strip().lower() == 'ph' or k.startswith('pH')]
        if pk:
            ph = np.array(g[pk[0]])[:, 0, :]
    return gas, carb, ph


def make_case_figure(runs, case_letter, data_override=None):
    if data_override is not None:
        data = data_override
    else:
        data = []
        for label, rel in runs:
            gas, carb, ph = read_fields(rel)
            if gas is None:
                continue
            data.append({'label': label, 'gas': gas, 'carb': carb, 'ph': ph})
    if not data:
        print(f"  No data for Case {case_letter}")
        return

    n = len(data)
    x_idx = int(np.searchsorted(X_CENTERS, X_VIEW_MAX))
    xe = X_EDGES[:x_idx+1]
    ze = Z_EDGES

    gas_max = max(d['gas'][:x_idx].max() for d in data)
    carb_lo = min(d['carb'][:x_idx].min() for d in data)
    carb_hi = max(d['carb'][:x_idx].max() for d in data)
    ph_lo = min(d['ph'][:x_idx].min() for d in data)
    ph_hi = max(d['ph'][:x_idx].max() for d in data)
    # Guard against zero-width color ranges
    if carb_hi - carb_lo < 1e-30:
        carb_hi = carb_lo + 1e-12
    if ph_hi - ph_lo < 1e-9:
        ph_hi = ph_lo + 1e-3

    fig, axes = plt.subplots(3, n, figsize=(2.4*n + 1.2, 6.8),
                             sharex=True, sharey=True)
    if n == 1:
        axes = axes.reshape(3, 1)
    fig.subplots_adjust(left=0.16, right=0.86, top=0.87, bottom=0.08,
                        hspace=0.22, wspace=0.10)

    row_specs = [
        ('Gas saturation', 'gas', 'viridis', 0.0, gas_max,
         f'0 \u2013 {gas_max:.2f}'),
        ('Net carbonate VF', 'carb', 'BuPu', carb_lo, carb_hi,
         f'{carb_lo*1e5:.3f}\u2013{carb_hi*1e5:.3f}\n\u00d710\u207b\u2075'),
        ('pH', 'ph', 'cividis', ph_lo, ph_hi,
         f'{ph_lo:.3f}\u2013{ph_hi:.3f}'),
    ]

    ims = []
    for row, (rlabel, key, cmap, vmin, vmax, rng) in enumerate(row_specs):
        im_row = None
        for col, d in enumerate(data):
            ax = axes[row, col]
            field = d[key][:x_idx, :].T
            im_row = ax.pcolormesh(xe, ze, field, cmap=cmap,
                                   vmin=vmin, vmax=vmax, shading='flat')
            # NOTE: no invert_yaxis — PFLOTRAN z increases upward (z=0 base,
            # z=100 top), so the buoyant gas cap correctly appears at the TOP.
            if row == 0:
                ax.set_title(d['label'], fontsize=8)
            if row == 2:
                ax.set_xlabel('Distance (m)', fontsize=8)
            if col == 0:
                ax.set_ylabel('Elevation (m)', fontsize=7.5)
            ax.tick_params(labelsize=6.5)
        ims.append((row, im_row, rlabel, rng))

    # Robust placement using actual axis positions
    fig.canvas.draw()
    for row, im_row, rlabel, rng in ims:
        pos_last = axes[row, -1].get_position()
        pos_first = axes[row, 0].get_position()
        cax = fig.add_axes([pos_last.x1 + 0.012, pos_last.y0,
                            0.015, pos_last.height])
        cb = fig.colorbar(im_row, cax=cax)
        cb.ax.tick_params(labelsize=6.5)
        ycen = (pos_first.y0 + pos_first.y1) / 2.0
        fig.text(0.035, ycen, f'{rlabel}\n({rng})',
                 rotation=90, va='center', ha='center',
                 fontsize=8, fontweight='bold')

    fig.suptitle(
        f'Case {case_letter}: gas saturation varies across runs, '
        f'carbonate and pH do not',
        fontsize=10, y=0.985)

    for ext in ('png', 'pdf'):
        out = OUT / f'fig_fields_case{case_letter}.{ext}'
        fig.savefig(out, dpi=(150 if ext == 'png' else 300),
                    bbox_inches='tight', pad_inches=0.06)
    plt.close(fig)
    print(f"  \u2713 fig_fields_case{case_letter}.png / .pdf")
    print(f"    Gas saturation : 0 \u2013 {gas_max:.3f}  (varies clearly across runs)")
    sp = 100*(carb_hi-carb_lo)/carb_hi if carb_hi else 0
    print(f"    Net carbonate  : {carb_lo*1e5:.4f} \u2013 {carb_hi*1e5:.4f} \u00d710\u207b\u2075  "
          f"(spread {sp:.3f}% of value)")
    print(f"    pH             : {ph_lo:.3f} \u2013 {ph_hi:.3f}  "
          f"(spread {ph_hi-ph_lo:.3f} pH units)")


if __name__ == '__main__':
    print("="*70)
    print("  Spatial Field Visualization \u2014 Cases D and E")
    print("="*70)
    print(f"  Base: {BASE}")
    print(f"  Showing near-well region x < {X_VIEW_MAX:.0f} m\n")

    print("Case D:")
    make_case_figure(CASE_D_RUNS, 'D')
    print("\nCase E:")
    make_case_figure(CASE_E_RUNS, 'E')

    print(f"""
{'='*70}
  HOW TO READ THESE FIGURES
{'='*70}
  Three rows share a color scale across columns:

    Row 1 (Gas saturation): the CAUSE. This visibly changes across
      columns \u2014 gas sits at different depths (D) or fills more pore
      space (E). Range ~0 to 0.7.

    Row 2 (Net carbonate): the EFFECT. Columns look essentially
      IDENTICAL; the range annotation shows sub-percent variation.

    Row 3 (pH): uniform at ~8.4 (alkaline); sub-0.1-pH-unit spread.

  Visual takeaway: the free-gas phase redistributes substantially, but
  carbonation and pH are insensitive to it \u2014 the null result for
  buoyancy (D) and contact (E), shown directly.

  Send the two PNGs back and we can view them together.
""")
    print("Done.")
