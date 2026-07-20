#!/usr/bin/env python3
"""Regenerate Figure 9 (single-panel Da-Sigma physical-regime diagram) with the two moved
corner labels. ALL file reading is delegated to make_manuscript_figures.py -- this script only
computes Sigma and draws the figure, so the data path stays identical to the main pipeline.

  * 'Transport-limited (self-sealing risk)'  moved right  (x 0.03 -> 0.11)
  * 'Reaction-limited (low carbonation)'      moved up     (y 0.03 -> 0.12)

Run from the WAG dir (so `import make_manuscript_figures` resolves and ROOT points at the data):
    cd $MYSCRATCH/WAG && "$PY" fig_da_sigma_regime.py
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

import make_manuscript_figures as mmf   # reuse ROOT, RSC, MU, read_run, rdir, _rate_curve, save, OUT

def main():
    R_, T_, Vm = 8.314, 333.15, 3.69e-5
    sig3 = 0.7 * 2900 * 9.81 * 600; p0 = 1000 * 9.81 * 600
    KOM = 1e5   # Omega proxy scaling -- matches mmf.fig_da_sigma (recalibrated for zero-seeded
                # carbonate VF). Set to 1e4 for the older (pre-zero-seeding) Sigma scale.

    # --- file reading via the main pipeline: total carbonate VF at final time per (scenario, rate) ---
    curves = mmf._rate_curve(lambda r: r["series"]["carb_mean"][-1])

    fig, ax = plt.subplots(figsize=(4.5, 4.2)); fig.subplots_adjust(bottom=0.22)
    n = 0
    for sc, (mus, vs, lbl, c, mk) in curves.items():
        dav, sigv = [], []
        for mu, tc in zip(mus, vs):
            om = max(1.01, 1 + tc * KOM)
            Pc = (R_ * T_ / Vm) * np.log(om)
            se = sig3 - p0 - 3e6 * min(mu, 5)
            sigv.append(Pc / max(se, 1e5)); dav.append(10.0 / mu)
        ax.plot(dav, sigv, color=c, ls='-', lw=0.6, alpha=0.4, zorder=3)
        ax.scatter(dav, sigv, c=c, s=28, marker=mk, edgecolors='white',
                   linewidths=0.4, label=lbl, zorder=5)
        n += len(dav)

    ax.axvline(3.0, color='#e2e8f0', ls='-', lw=0.5, zorder=1)
    ax.axhline(800, color='#e2e8f0', ls='-', lw=0.5, zorder=1)

    bbox_kw = dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85)
    kw = dict(fontsize=5.5, fontstyle='italic', zorder=6)
    ax.text(0.11, 0.97, 'Transport-limited\n(self-sealing risk)', transform=ax.transAxes,
            va='top', color='#dc2626', bbox=bbox_kw, **kw)                     # moved right
    ax.text(0.97, 0.97, 'Distributed\nprecipitation', transform=ax.transAxes,
            va='top', ha='right', color='#16a34a', bbox=bbox_kw, **kw)
    ax.text(0.03, 0.03, 'Clogging risk', transform=ax.transAxes,
            color='#d97706', bbox=bbox_kw, **kw)
    ax.text(0.97, 0.12, 'Reaction-limited\n(low carbonation)', transform=ax.transAxes,
            ha='right', color='#64748b', bbox=bbox_kw, **kw)                   # moved up

    # arrow, with its label rotated parallel to it and sitting just above the line
    ax1f, ay1f, ax2f, ay2f = 0.52, 0.68, 0.75, 0.42    # tail -> head (axes fraction)
    ax.annotate('', xy=(ax2f, ay2f), xytext=(ax1f, ay1f),
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='#94a3b8', lw=0.7))
    (PX1, PY1) = ax.transAxes.transform((ax1f, ay1f))
    (PX2, PY2) = ax.transAxes.transform((ax2f, ay2f))
    ang = np.degrees(np.arctan2(PY2 - PY1, PX2 - PX1))
    ax.text((ax1f + ax2f) / 2, (ay1f + ay2f) / 2, 'Increasing injection rate',
            transform=ax.transAxes, rotation=ang, rotation_mode='anchor',
            ha='center', va='bottom', fontsize=5.5, color='#94a3b8')

    ax.set_xscale('log')
    ax.set_xlabel('Damk\u00f6hler number (Da)')
    ax.set_ylabel('Normalised stress ratio (\u03a3)')
    if n:
        ax.legend(fontsize=5.5, ncol=6, loc='lower center', bbox_to_anchor=(0.5, -0.22),
                  frameon=False, columnspacing=0.8, handletextpad=0.3, markerscale=1.0)
    else:
        print("  WARNING: mmf._rate_curve returned no runs -- empty plot (check 08_rate_sweep/runs)")

    mmf.save(fig, "fig_da_sigma_regime")   # writes figures/fig_da_sigma_regime.pdf (+ .png)
    print(f"  points plotted: {n}")

if __name__ == "__main__":
    main()
