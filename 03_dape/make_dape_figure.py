#!/usr/bin/env python3
"""
make_dape_figure.py

Generates fig_dape_disentangling.pdf — the Da-Pe disentangling figure
that responds to peer review concerns about whether the rate-sweep
trends are governed by Da or by associated changes in Pe.

The figure has two panels:
  Panel (a) Case A: carbonate VF vs kinetic-rate multiplier κ
            (Da varies, Pe held fixed)
  Panel (b) Case B: carbonate VF vs Péclet number Pe
            (Pe varies, Da held fixed)

Data values are from 18 PFLOTRAN simulations on hpc01, varying κ
(intrinsic rate constant multiplier) for Case A and varying D
(molecular diffusion coefficient) for Case B.

Usage:
    python3 make_dape_figure.py

Outputs (in current working directory):
    fig_dape_disentangling.pdf
    fig_dape_disentangling.png

Dependencies:
    numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# ---------------------------------------------------------------------
# Plot style — journal-ready, serif typeface
# ---------------------------------------------------------------------
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 8
mpl.rcParams['axes.linewidth'] = 0.6
mpl.rcParams['xtick.major.width'] = 0.6
mpl.rcParams['ytick.major.width'] = 0.6
mpl.rcParams['xtick.major.size'] = 3
mpl.rcParams['ytick.major.size'] = 3

# ---------------------------------------------------------------------
# Simulation results (carbonate volume fraction at 30 years)
# ---------------------------------------------------------------------

# Case A: kinetic-rate sweep at fixed flow (varies Da, holds Pe)
caseA_dissolved = {
    0.1:  8.7890e-04,
    0.3:  9.0064e-04,
    1.0:  8.9321e-04,
    3.0:  8.9466e-04,
    10.0: 9.0891e-04,
}
caseA_scco2 = {
    0.1:  4.1677e-04,
    0.3:  4.3591e-04,
    1.0:  4.2546e-04,
    3.0:  4.2540e-04,
    10.0: 4.0000e-04,
}

# Case B: diffusion-coefficient sweep at fixed flow + kinetics
# (varies Pe, holds Da)
caseB_dissolved = {
    1.0e-10: 8.9321e-04,  # Pe = 1e7
    1.0e-9:  8.9321e-04,  # Pe = 1e6
    1.0e-8:  8.9321e-04,  # Pe = 1e5
    1.0e-7:  8.9320e-04,  # Pe = 1e4
}
caseB_scco2 = {
    1.0e-10: 4.2546e-04,
    1.0e-9:  4.2546e-04,
    1.0e-8:  4.2546e-04,
    1.0e-7:  4.2546e-04,
}

# Convert D to Pe (Q_BASE = base Darcy velocity, L_CHAR = characteristic length)
Q_BASE = 1.0e-5
L_CHAR = 100.0
def pe(D):
    return Q_BASE * L_CHAR / D


# ---------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
fig.subplots_adjust(left=0.09, right=0.97, top=0.88, bottom=0.18, wspace=0.30)

dissolved_color = '#1f6fb4'
scco2_color     = '#c92e1d'

# ===== Panel (a): Case A — Da varies, Pe fixed =====
ax = axes[0]
k_d = sorted(caseA_dissolved.keys())
v_d = [caseA_dissolved[k] for k in k_d]
k_s = sorted(caseA_scco2.keys())
v_s = [caseA_scco2[k] for k in k_s]

ax.semilogx(k_d, np.array(v_d) * 1e4,
            marker='o', color=dissolved_color, markersize=5.5,
            lw=1.2, mec='white', mew=0.6,
            label='Dissolved (S1)')
ax.semilogx(k_s, np.array(v_s) * 1e4,
            marker='s', color=scco2_color, markersize=5.5,
            lw=1.2, mec='white', mew=0.6,
            label='Supercritical (S2)')

# Reference lines at baseline (κ = 1) values
ax.axhline(caseA_dissolved[1.0] * 1e4, color=dissolved_color,
           lw=0.5, ls=':', alpha=0.5)
ax.axhline(caseA_scco2[1.0] * 1e4, color=scco2_color,
           lw=0.5, ls=':', alpha=0.5)

ax.set_xlabel(r'Kinetic-rate multiplier $\kappa$', fontsize=8.5)
ax.set_ylabel(r'Carbonate volume fraction at 30 yr ($\times 10^{-4}$)', fontsize=8.5)
ax.set_ylim(2.5, 11)
ax.set_xticks([0.1, 0.3, 1, 3, 10])
ax.set_xticklabels(['0.1', '0.3', '1', '3', '10'])
ax.legend(loc='upper left', fontsize=7.5, frameon=False)
ax.text(0.96, 0.95, '(a)', transform=ax.transAxes, fontsize=9,
        fontweight='bold', va='top', ha='right')
ax.set_title(r'Case A: $Da$ varies, $Pe$ fixed', fontsize=8.5, pad=4)

# Flatness annotations
ax.text(0.96, 0.62, r'$\pm 3.4\%$ over 100$\times$ in $\kappa$',
        transform=ax.transAxes, fontsize=6.5, ha='right', va='center',
        color=dissolved_color)
ax.text(0.96, 0.27, r'$\pm 9.0\%$ over 100$\times$ in $\kappa$',
        transform=ax.transAxes, fontsize=6.5, ha='right', va='center',
        color=scco2_color)


# ===== Panel (b): Case B — Pe varies, Da fixed =====
ax = axes[1]
D_d = sorted(caseB_dissolved.keys())
pe_d = [pe(D) for D in D_d]
v_bd = [caseB_dissolved[D] for D in D_d]
D_s = sorted(caseB_scco2.keys())
pe_s = [pe(D) for D in D_s]
v_bs = [caseB_scco2[D] for D in D_s]

order_d = np.argsort(pe_d)
pe_d_sorted = [pe_d[i] for i in order_d]
v_bd_sorted = [v_bd[i] for i in order_d]
order_s = np.argsort(pe_s)
pe_s_sorted = [pe_s[i] for i in order_s]
v_bs_sorted = [v_bs[i] for i in order_s]

ax.semilogx(pe_d_sorted, np.array(v_bd_sorted) * 1e4,
            marker='o', color=dissolved_color, markersize=5.5,
            lw=1.2, mec='white', mew=0.6,
            label='Dissolved (S1)')
ax.semilogx(pe_s_sorted, np.array(v_bs_sorted) * 1e4,
            marker='s', color=scco2_color, markersize=5.5,
            lw=1.2, mec='white', mew=0.6,
            label='Supercritical (S2)')

# Reference lines at the Case A baseline (visual continuity with panel a)
ax.axhline(caseA_dissolved[1.0] * 1e4, color=dissolved_color,
           lw=0.5, ls=':', alpha=0.5)
ax.axhline(caseA_scco2[1.0] * 1e4, color=scco2_color,
           lw=0.5, ls=':', alpha=0.5)

ax.set_xlabel('Péclet number $Pe = qL/D$', fontsize=8.5)
ax.set_ylabel(r'Carbonate volume fraction at 30 yr ($\times 10^{-4}$)', fontsize=8.5)
ax.set_ylim(2.5, 11)
ax.set_xticks([1e4, 1e5, 1e6, 1e7])
ax.legend(loc='upper left', fontsize=7.5, frameon=False)
ax.text(0.96, 0.95, '(b)', transform=ax.transAxes, fontsize=9,
        fontweight='bold', va='top', ha='right')
ax.set_title(r'Case B: $Pe$ varies, $Da$ fixed', fontsize=8.5, pad=4)

# Annotation
ax.text(0.96, 0.50,
        'Both curves constant\nto 4 sig. figs.\nover 1000$\\times$ in $D$',
        transform=ax.transAxes, fontsize=6.5, ha='right', va='center',
        color='#444444')

# Shade the Pe range spanned by the original rate sweep
ax.axvspan(1e5, 3e7, alpha=0.06, color='gray', zorder=0)
ax.text(0.5, 0.05, r'shaded: $Pe$ range of rate sweep',
        transform=ax.transAxes, fontsize=5.5, ha='center', va='bottom',
        color='gray', style='italic')


# ---------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------
plt.savefig('fig_dape_disentangling.pdf', dpi=300,
            bbox_inches='tight', pad_inches=0.05)
plt.savefig('fig_dape_disentangling.png', dpi=200,
            bbox_inches='tight', pad_inches=0.05)
print('Saved: fig_dape_disentangling.pdf and fig_dape_disentangling.png')
