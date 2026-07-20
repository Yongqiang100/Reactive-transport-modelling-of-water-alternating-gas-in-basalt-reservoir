#!/usr/bin/env python3
"""
Master Figure Generator — Nature-style minimal aesthetics.

WAG CO₂ Mineralization in Basalt (Journal of Hydrology)
Reads: ~/WAG/, ~/WAG/WAG-damkohler/, ~/WAG/WAG-KineticSensitivity/
Outputs: ~/WAG/figures/ (PDF + PNG)
"""
import os, glob, json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap

# ── Paths ──────────────────────────────────────────────────
HOME = os.path.expanduser("~")
WAG  = os.path.join(HOME, "WAG")
DA   = os.path.join(WAG, "WAG-damkohler")
KS   = os.path.join(WAG, "WAG-KineticSensitivity")
OUT  = os.path.join(WAG, "figures"); os.makedirs(OUT, exist_ok=True)

# ── Scenario metadata ─────────────────────────────────────
SC = [
    ("scenario1_dissolved", "Dissolved",    "#1f77b4"),
    ("scenario2_scco2",     "scCO₂",       "#d62728"),
    ("scenario3_wag6mo",    "WAG-6 mo",     "#2ca02c"),
    ("scenario4_wag3mo",    "WAG-3 mo",     "#9467bd"),
    ("scenario5_swag",      "SWAG",         "#ff7f0e"),
    ("scenario6_adaptive",  "Adaptive WAG", "#17becf"),
]
DA_LEVELS = ["da_high","da_base","da_medium","da_medlow","da_low"]
DA_MULT   = [0.3, 1.0, 3.0, 10.0, 30.0]
CARB = ['Calcite','Magnesite','Siderite','Dolomite-ord']
PRIM = ['Forsterite','Anorthite','Diopside']
ALL_MIN = CARB + PRIM + ['Kaolinite','SiO2(am)']

# ── Nature-style rcParams ─────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica','Arial','DejaVu Sans'],
    'font.size': 7,
    'axes.labelsize': 8,
    'axes.titlesize': 8,
    'axes.linewidth': 0.6,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 6,
    'legend.frameon': False,
    'legend.handlelength': 1.5,
    'lines.linewidth': 1.0,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})
LW = 1.0
LW_THICK = 1.4

def _label(ax, text, x=-0.12, y=1.08):
    ax.text(x, y, text, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top')

def _save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p + '.pdf')
    fig.savefig(p + '.png', dpi=150)
    plt.close(fig)
    print(f"  ✔ {p}.pdf  +  .png")


# ── HDF5 reader ───────────────────────────────────────────
try:
    import h5py; HAS_H5 = True
except ImportError:
    HAS_H5 = False; print("WARNING: h5py unavailable — JSON only")

def read_h5(sdir):
    if not HAS_H5: return None
    h5s = sorted([f for f in glob.glob(os.path.join(sdir,"*.h5"))
                  if 'hanford' not in os.path.basename(f)])
    if not h5s: return None
    d = {k: [] for k in ['times_yr','porosity_mean','porosity_nw','porosity_ff',
                          'ph_mean','ph_nw','ph_ff','perm_mean','gas_sat_mean']}
    for m in ALL_MIN:
        for s in ['_vf_mean','_vf_nw','_vf_ff']:
            d[f'{m}{s}'] = []
    d['spatial'] = None
    nw, ff = slice(0,25), slice(50,None)
    for h5f in h5s:
        try:
            with h5py.File(h5f,'r') as f:
                tgs = sorted([g for g in f.keys() if g.startswith('Time')],
                             key=lambda g: float(g.split(':')[1].strip().split()[0]))
                for gn in tgs:
                    g = f[gn]; t = float(gn.split(':')[1].strip().split()[0])
                    d['times_yr'].append(t)
                    if 'Porosity' in g:
                        a=g['Porosity'][:]
                        d['porosity_mean'].append(np.mean(a))
                        d['porosity_nw'].append(np.mean(a[nw,:,:]))
                        d['porosity_ff'].append(np.mean(a[ff,:,:]))
                    for pk in ['Permeability_X [m^2]','Permeability_X']:
                        if pk in g: d['perm_mean'].append(np.mean(g[pk][:])); break
                    ph_ok=False
                    for pk in ['pH','ph','PH']:
                        if pk in g:
                            a=g[pk][:]
                            d['ph_mean'].append(np.mean(a))
                            d['ph_nw'].append(np.mean(a[nw,:,:]))
                            d['ph_ff'].append(np.mean(a[ff,:,:])); ph_ok=True; break
                    if not ph_ok:
                        for k2 in g.keys():
                            if 'pH' in k2:
                                a=g[k2][:]
                                if a.ndim>=2:
                                    d['ph_mean'].append(np.mean(a))
                                    d['ph_nw'].append(np.mean(a[nw,:,:]))
                                    d['ph_ff'].append(np.mean(a[ff,:,:])); break
                    for gk in ['Gas_Saturation','Gas Saturation']:
                        if gk in g: d['gas_sat_mean'].append(np.mean(g[gk][:])); break
                    for m in ALL_MIN:
                        ok=False
                        for sf in [f'{m}_VF [m^3 mnrl_m^3 bulk]',f'{m}_VF']:
                            if sf in g:
                                a=g[sf][:]
                                d[f'{m}_vf_mean'].append(np.mean(a))
                                d[f'{m}_vf_nw'].append(np.mean(a[nw,:,:]))
                                d[f'{m}_vf_ff'].append(np.mean(a[ff,:,:])); ok=True; break
                        if not ok:
                            d[f'{m}_vf_mean'].append(np.nan)
                            d[f'{m}_vf_nw'].append(np.nan)
                            d[f'{m}_vf_ff'].append(np.nan)
                    if gn==tgs[-1]:
                        sp={'x': np.arange(g['Porosity'].shape[0]) if 'Porosity' in g else np.arange(140)}
                        if 'Porosity' in g: sp['phi']=g['Porosity'][:,0,g['Porosity'].shape[2]//2]
                        for pk in ['pH','ph']:
                            if pk in g: sp['ph']=g[pk][:,0,g[pk].shape[2]//2]; break
                        for m in CARB:
                            for sf in [f'{m}_VF [m^3 mnrl_m^3 bulk]',f'{m}_VF']:
                                if sf in g: sp[f'{m}']=g[sf][:,0,g[sf].shape[2]//2]; break
                        d['spatial']=sp
        except Exception as e:
            print(f"  warn {h5f}: {e}")
    for k in d:
        if k!='spatial' and isinstance(d[k],list): d[k]=np.array(d[k])
    return d

def load_json(fp):
    if os.path.exists(fp):
        with open(fp) as f: return json.load(f)
    return None

def _carb_sum(d, suffix='_vf_mean'):
    t=d.get('times_yr',[])
    c=np.zeros(len(t))
    for m in CARB:
        a=d.get(f'{m}{suffix}',np.zeros(len(t)))
        if len(a)==len(t): c+=np.nan_to_num(a)
    return c


# ══════════════════════════════════════════════════════════
#  FIGURE: Domain Schematic
# ══════════════════════════════════════════════════════════
def fig_domain():
    print("  Domain schematic...")
    fig, ax = plt.subplots(figsize=(7.2, 2.8))

    z1, z2, z3 = 3, 6, 12
    Z = 5

    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((0,0), z1, Z, fc='#fef2f2', ec='none', zorder=0))
    ax.add_patch(Rectangle((z1,0), z2-z1, Z, fc='#fefce8', ec='none', zorder=0))
    ax.add_patch(Rectangle((z2,0), z3-z2, Z, fc='#eff6ff', ec='none', zorder=0))

    ax.plot([0,z3,z3,0,0],[0,0,Z,Z,0], color='#334155', lw=0.8, zorder=3)

    for xb in [z1, z2]:
        ax.plot([xb,xb],[0,Z], color='#94a3b8', ls='--', lw=0.5, zorder=2)

    well_bot, well_top = 1.0, 4.0
    ax.add_patch(Rectangle((0, well_bot), 0.3, well_top-well_bot,
                            fc='#fca5a5', ec='#dc2626', lw=1.0, zorder=4))

    akw = dict(arrowstyle='->', mutation_scale=10)
    for z in [2.0, 3.0]:
        ax.annotate('', xy=(2.5, z), xytext=(0.4, z),
                    arrowprops={**akw, 'color':'#dc2626', 'lw':0.9}, zorder=5)
    for z in [2.0, 3.0]:
        ax.annotate('', xy=(5.5, z), xytext=(3.2, z),
                    arrowprops={**akw, 'color':'#d97706', 'lw':0.8, 'alpha':0.7}, zorder=5)
    ax.annotate('', xy=(11.0, 2.5), xytext=(6.3, 2.5),
                arrowprops={**akw, 'color':'#3b82f6', 'lw':0.7, 'alpha':0.4}, zorder=5)

    ax.plot([z3, z3], [0, Z], color='#2563eb', lw=2.5, zorder=4)

    ax.text(z1/2, Z-0.3, 'Reactive', fontsize=7, ha='center', va='top',
            color='#dc2626', fontstyle='italic')
    ax.text((z1+z2)/2, Z-0.3, 'Transition', fontsize=7, ha='center', va='top',
            color='#d97706', fontstyle='italic')
    ax.text((z2+z3)/2, Z-0.3, 'Near-equilibrium', fontsize=7, ha='center', va='top',
            color='#2563eb', fontstyle='italic')

    ax.text(z1/2, 2.5, 'pH 4.5–6', fontsize=6, ha='center', va='center',
            color='#dc2626', alpha=0.6)
    ax.text((z1+z2)/2, 2.5, 'pH 6–7.5', fontsize=6, ha='center', va='center',
            color='#d97706', alpha=0.6)
    ax.text((z2+z3)/2, 2.5-0.5, 'pH 8+, Q/K ≈ 1', fontsize=6, ha='center', va='center',
            color='#2563eb', alpha=0.6)

    ax.text(z1/2, -0.4, '0–100 m\nΔx = 2 m, 50 cells', fontsize=5.5,
            ha='center', va='top', color='#64748b')
    ax.text((z1+z2)/2, -0.4, '100–500 m\nΔx = 10 m, 40 cells', fontsize=5.5,
            ha='center', va='top', color='#64748b')
    ax.text((z2+z3)/2, -0.4, '500–2000 m\nΔx = 30 m, 50 cells', fontsize=5.5,
            ha='center', va='top', color='#64748b')

    ax.text(-0.5, well_top, 'z = 80 m', fontsize=6, ha='right', va='center', color='#64748b')
    ax.text(-0.5, well_bot, 'z = 20 m', fontsize=6, ha='right', va='center', color='#64748b')
    ax.text(-0.5, 2.5, '60 °C', fontsize=7.5, ha='right', va='center',
            color='#dc2626', fontweight='bold')

    ax.text(0.15, Z+0.3, 'Well', fontsize=7, ha='center', va='bottom',
            fontweight='bold', color='#dc2626')
    ax.text(z3, Z+0.3, 'x = 2 km', fontsize=6,
            ha='right', va='bottom', color='#64748b')

    ax.text(z3+0.3, 2.5, 'P = 6 MPa\nOutlet BC', fontsize=7, ha='left', va='center',
            color='#2563eb', fontweight='bold')

    ax.set_xlim(-1.2, z3+1.8)
    ax.set_ylim(-1.5, Z+0.8)
    ax.axis('off')
    _save(fig, 'fig_domain_schematic')


# ══════════════════════════════════════════════════════════
#  FIGURE: Kinetic Sensitivity (9-panel)
# ══════════════════════════════════════════════════════════
def fig_kin_sens():
    print("  Kinetic sensitivity...")
    js = load_json(os.path.join(KS,"kinetic_sensitivity.json"))
    pk_h5 = read_h5(os.path.join(KS,"kinetics_pk2004"))
    rm_h5 = read_h5(os.path.join(KS,"kinetics_rimstidt2012"))
    pk, rm = {}, {}
    if js:
        for jk, tgt in [("kinetics_pk2004",pk),("kinetics_rimstidt2012",rm)]:
            if jk in js:
                jd=js[jk]
                tgt['t']=np.array(jd['times_yr'])
                for f in ['porosity_nw','ph_nw']:
                    if f in jd: tgt[f]=np.array(jd[f])
                for m in CARB+PRIM:
                    k2=f'{m}_vf_nw'
                    if k2 in jd: tgt[k2]=np.array(jd[k2])
                if 'metrics' in jd: tgt['metrics']=jd['metrics']
    for tgt,h5 in [(pk,pk_h5),(rm,rm_h5)]:
        if tgt and h5:
            nt = len(tgt.get('t',[]))
            for f in ['porosity_mean','ph_mean']:
                if f in h5 and len(h5[f])==nt: tgt[f]=h5[f]
            for m in CARB+PRIM:
                k2=f'{m}_vf_mean'
                if k2 in h5 and len(h5.get(k2,[]))==nt: tgt[k2]=h5[k2]

    cases = [("Palandri & Kharaka (2004)", pk, "#1f77b4"),
             ("Rimstidt et al. (2012)",     rm, "#d62728")]

    fig = plt.figure(figsize=(7.2, 7.2))
    gs = GridSpec(3,3,figure=fig,hspace=0.18,wspace=0.40)
    fig.subplots_adjust(bottom=0.10)

    def _safe(ax, d, key, **kw):
        t = d.get('t', [])
        y = d.get(key, [])
        if len(y) == len(t) > 0:
            return ax.plot(t, y, **kw)
        return []

    def _carb_nw(d):
        t = d.get('t', [])
        c = np.zeros(len(t))
        for m in CARB:
            a = d.get(f'{m}_vf_nw', np.zeros(len(t)))
            if len(a) == len(t): c += np.nan_to_num(a)
        return t, c

    def _carb_avg(d):
        t = d.get('t', [])
        c = np.zeros(len(t))
        for m in CARB:
            a = d.get(f'{m}_vf_mean', [])
            if len(a) == len(t):
                c += np.nan_to_num(a)
            else:
                a2 = d.get(f'{m}_vf_nw', np.zeros(len(t)))
                if len(a2) == len(t): c += np.nan_to_num(a2)
        return t, c

    legend_lines = []

    ax_a = fig.add_subplot(gs[0,0]); ax_b = fig.add_subplot(gs[0,1], sharex=ax_a); ax_c = fig.add_subplot(gs[0,2], sharex=ax_a)
    for lb,d,c in cases:
        t,y = _carb_nw(d)
        if len(t) > 0:
            l, = ax_a.plot(t,y,color=c,lw=LW,label=lb)
            legend_lines.append(l)
    ax_a.set_ylabel('Volume fraction'); _label(ax_a,'a')
    ax_a.tick_params(labelbottom=False)

    for lb,d,c in cases: _safe(ax_b,d,'Forsterite_vf_nw',color=c,lw=LW)
    ax_b.set_ylabel('Volume fraction'); _label(ax_b,'b')
    ax_b.tick_params(labelbottom=False)

    for lb,d,c in cases: _safe(ax_c,d,'ph_nw',color=c,lw=LW)
    ax_c.set_ylabel('pH'); _label(ax_c,'c')
    ax_c.tick_params(labelbottom=False)

    ax_d = fig.add_subplot(gs[1,0], sharex=ax_a); ax_e = fig.add_subplot(gs[1,1], sharex=ax_a); ax_f = fig.add_subplot(gs[1,2], sharex=ax_a)
    for lb,d,c in cases:
        t,y = _carb_avg(d)
        if len(t) > 0: ax_d.plot(t,y,color=c,lw=LW)
    ax_d.set_ylabel('Volume fraction'); _label(ax_d,'d')
    ax_d.ticklabel_format(axis='y',style='sci',scilimits=(-3,3))
    ax_d.tick_params(labelbottom=False)

    for lb,d,c in cases:
        t = d.get('t', [])
        y = d.get('Forsterite_vf_mean', d.get('Forsterite_vf_nw', []))
        if len(y) == len(t) > 0: ax_e.plot(t,y,color=c,lw=LW)
    ax_e.set_ylabel('Volume fraction'); _label(ax_e,'e')
    ax_e.tick_params(labelbottom=False)

    for lb,d,c in cases: _safe(ax_f,d,'porosity_nw',color=c,lw=LW)
    ax_f.set_ylabel('Porosity'); _label(ax_f,'f')
    ax_f.tick_params(labelbottom=False)

    ax_g = fig.add_subplot(gs[2,0], sharex=ax_a); ax_h = fig.add_subplot(gs[2,1], sharex=ax_a)
    for lb,d,c in cases: _safe(ax_g,d,'Magnesite_vf_nw',color=c,lw=LW)
    ax_g.set_ylabel('Volume fraction'); ax_g.set_xlabel('Time (yr)')
    ax_g.ticklabel_format(axis='y',style='sci',scilimits=(-3,3)); _label(ax_g,'g')

    for lb,d,c in cases: _safe(ax_h,d,'Dolomite-ord_vf_nw',color=c,lw=LW)
    ax_h.set_ylabel('Volume fraction'); ax_h.set_xlabel('Time (yr)')
    ax_h.ticklabel_format(axis='y',style='sci',scilimits=(-3,3)); _label(ax_h,'h')

    ax_i = fig.add_subplot(gs[2,2]); ax_i.axis('off')
    pk_m=pk.get('metrics',{}); rm_m=rm.get('metrics',{})
    cd=abs(rm_m.get('carb_nw',0)-pk_m.get('carb_nw',1))/max(pk_m.get('carb_nw',1),1e-12)*100
    fd=abs(rm_m.get('fo_dissolved_pct',0)-pk_m.get('fo_dissolved_pct',0))
    txt=(f"Rate ratio at pH 8.5:\n  Rimstidt / P&K ≈ 1 500×\n\n"
         f"Output difference (30 yr):\n  Carbonate VF  {cd:.1f}%\n"
         f"  Forsterite    {fd:.2f}% pts\n\n"
         f"→ Transport-limited\n  (1 − Q/K) → 0")
    ax_i.text(0.05,0.95,txt,transform=ax_i.transAxes,fontsize=6.5,va='top',
              fontfamily='monospace',
              bbox=dict(boxstyle='round,pad=0.4',fc='#f8fafc',ec='#cbd5e1',lw=0.5))
    _label(ax_i,'i')

    fig.legend(handles=legend_lines, ncol=2, fontsize=7,
               loc='lower center', bbox_to_anchor=(0.5, 0.01),
               frameon=False, columnspacing=1.5, handlelength=2.0)
    _save(fig,'fig_kinetic_sensitivity')


# ══════════════════════════════════════════════════════════
#  FIGURE: Baseline Comparison (6-panel)
# ══════════════════════════════════════════════════════════
def fig_baseline():
    print("  Baseline comparison...")
    data={}
    for dn,lb,c in SC:
        d=read_h5(os.path.join(WAG,dn))
        if d and len(d.get('times_yr',[]))>0: data[dn]=d

    fig,axes=plt.subplots(2,3,figsize=(7.2,4.5), sharex=True)
    fig.subplots_adjust(hspace=0.15,wspace=0.40, bottom=0.18)

    configs=[
        (lambda d: _carb_sum(d),'Total carbonate VF'),
        ('ph_mean','pH'),
        ('porosity_nw','Near-well porosity'),
        ('Forsterite_vf_mean','Forsterite VF'),
        ('gas_sat_mean','Gas saturation'),
        (lambda d: _carb_sum(d,'_vf_nw'),'Near-well carbonate VF'),
    ]
    labels='abcdef'
    lines_for_legend = []

    for ax,(field,yl),lab in zip(axes.flat,configs,labels):
        for dn,lb,c in SC:
            if dn not in data: continue
            d=data[dn]; t=d['times_yr']
            y=field(d) if callable(field) else d.get(field,[])
            if len(y)==len(t):
                l, = ax.plot(t,y,color=c,lw=LW,label=lb)
                if lab=='a': lines_for_legend.append(l)
        ax.set_ylabel(yl, fontsize=7)
        _label(ax,lab)
        if lab in 'def':
            ax.set_xlabel('Time (yr)')
        else:
            ax.tick_params(labelbottom=False)

    fig.legend(handles=lines_for_legend, ncol=6, fontsize=6,
               loc='lower center', bbox_to_anchor=(0.5, 0.02),
               frameon=False, columnspacing=1.0, handlelength=1.5)
    _save(fig,'fig_baseline_comparison')


# ══════════════════════════════════════════════════════════
#  FIGURE: Spatial Profiles
# ══════════════════════════════════════════════════════════
def fig_spatial():
    print("  Spatial profiles...")
    fig,axes=plt.subplots(1,3,figsize=(8.5,2.8))
    fig.subplots_adjust(wspace=0.35, bottom=0.25)

    lines_for_legend = []
    for dn,lb,c in SC:
        d=read_h5(os.path.join(WAG,dn))
        if d is None or d.get('spatial') is None: continue
        sp=d['spatial']; x=sp.get('x',np.arange(140))
        if 'ph' in sp and len(sp['ph'])==len(x):
            l, = axes[0].plot(x,sp['ph'],color=c,lw=LW,label=lb)
            lines_for_legend.append(l)
        if 'phi' in sp and len(sp['phi'])==len(x):
            axes[1].plot(x,sp['phi'],color=c,lw=LW)
        ct=np.zeros(len(x))
        for m in CARB:
            if m in sp and len(sp[m])==len(x): ct+=sp[m]
        if np.any(ct>0): axes[2].plot(x,ct,color=c,lw=LW)

    # Zoom panels b and c to near-wellbore region
    axes[1].set_xlim(0, 30)
    axes[2].set_xlim(0, 20)

    for ax,yl,lab in zip(axes,['pH','Porosity','Carbonate VF'],'abc'):
        ax.set_ylabel(yl)
        _label(ax, lab, y=1.12)

    fig.text(0.5, 0.13, 'Cell index', fontsize=8)

    fig.legend(handles=lines_for_legend, ncol=6, fontsize=6,
               loc='lower center', bbox_to_anchor=(0.5, 0.0),
               frameon=False, columnspacing=1.0, handlelength=1.5)
    _save(fig,'fig_spatial_profiles')


# ══════════════════════════════════════════════════════════
#  FIGURE: Carbonate Breakdown
# ══════════════════════════════════════════════════════════
def fig_carb_bar():
    print("  Carbonate breakdown...")
    fig,ax=plt.subplots(figsize=(5.5,2.8))
    fig.subplots_adjust(bottom=0.22)
    labels=[]; cd={m:[] for m in CARB}
    cc={'Calcite':'#1f77b4','Magnesite':'#2ca02c','Siderite':'#ff7f0e','Dolomite-ord':'#9467bd'}

    for dn,lb,c in SC:
        d=read_h5(os.path.join(WAG,dn))
        if d is None: continue
        labels.append(lb)
        for m in CARB:
            v=d.get(f'{m}_vf_mean',[])
            cd[m].append(float(v[-1]) if len(v)>0 and not np.isnan(v[-1]) else 0)

    x=np.arange(len(labels)); w=0.5; bot=np.zeros(len(labels))
    for m in CARB:
        v=np.array(cd[m])
        ax.bar(x,v,w,bottom=bot,label=m,color=cc[m],edgecolor='white',linewidth=0.3)
        bot+=v
    ax.set_xticks(x); ax.set_xticklabels(labels,rotation=0,ha='center',fontsize=6.5)
    ax.set_ylabel('Volume fraction'); ax.legend(fontsize=5.5,ncol=4,loc='upper right')
    ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,3))
    _save(fig,'fig_carbonate_breakdown')


# ══════════════════════════════════════════════════════════
#  FIGURE: Damköhler Sweep (2-panel)
# ══════════════════════════════════════════════════════════
def fig_da_sweep():
    print("  Damköhler sweep...")
    da_json=load_json(os.path.join(DA,"damkohler_results.json"))
    fig,axes=plt.subplots(1,2,figsize=(7.2,2.8))
    fig.subplots_adjust(wspace=0.35)

    markers = ['o','s','D','^','v','P']

    for ax_i,(ylabel,get_val) in enumerate([
        ('Carbonate VF at 30 yr', lambda d: sum(float(d[f'{m}_vf_mean'][-1])
            for m in CARB if len(d.get(f'{m}_vf_mean',[]))>0
            and not np.isnan(d[f'{m}_vf_mean'][-1]))),
        ('Δφ (φ₃₀ − φ₀)', lambda d: float(d['porosity_mean'][-1]-d['porosity_mean'][0])
            if len(d.get('porosity_mean',[]))>1 else None)]):
        ax=axes[ax_i]
        for (dn,lb,c),mk in zip(SC,markers):
            pts=[]
            for dl,mu in zip(DA_LEVELS,DA_MULT):
                d=read_h5(os.path.join(DA,dl,dn))
                if d and len(d.get('times_yr',[]))>0:
                    v=get_val(d)
                    if v is not None: pts.append((mu,v))
                elif da_json:
                    k=f"{dl}/{dn}"
                    if k in da_json and ax_i==0:
                        pts.append((mu,da_json[k].get('total_carbonate_vf',0)))
            if pts:
                ms,vs=zip(*pts)
                ax.plot(ms,vs,color=c,lw=LW,marker=mk,markersize=4,
                        markeredgecolor='white',markeredgewidth=0.4,label=lb)
        ax.set_xscale('log'); ax.set_xlabel('Rate multiplier (Q/Q₀)')
        ax.set_ylabel(ylabel); _label(ax,'ab'[ax_i])
        if ax_i==1: ax.axhline(0,color='#94a3b8',lw=0.4)
    axes[0].legend(fontsize=5,ncol=2)
    _save(fig,'fig_damkohler_sweep')


# ══════════════════════════════════════════════════════════
#  FIGURE: Da–Σ Regime Diagram
# ══════════════════════════════════════════════════════════
def fig_da_sigma():
    print("  Da–Σ regime diagram...")
    R_=8.314; T_=333.15; Vm=3.69e-5
    sig3=0.7*2900*9.81*600; p0=1000*9.81*600

    markers = ['o','s','D','^','v','P']

    fig,ax=plt.subplots(figsize=(4.5,4.2))
    fig.subplots_adjust(bottom=0.22)

    for (dn,lb,c), mk in zip(SC, markers):
        dav,sigv=[],[]
        for dl,mu in zip(DA_LEVELS,DA_MULT):
            d=read_h5(os.path.join(DA,dl,dn))
            if d is None or len(d.get('times_yr',[]))==0: continue
            tc=sum(float(d[f'{m}_vf_mean'][-1]) for m in CARB
                   if len(d.get(f'{m}_vf_mean',[]))>0 and not np.isnan(d[f'{m}_vf_mean'][-1]))
            da=10./mu
            om=max(1.01,1+tc*1e4)
            Pc=(R_*T_/Vm)*np.log(om)
            dp=3e6*min(mu,5); se=sig3-p0-dp
            S=Pc/max(se,1e5)
            dav.append(da); sigv.append(S)
        if dav:
            ax.plot(dav, sigv, color=c, ls='-', lw=0.6, alpha=0.4, zorder=3)
            ax.scatter(dav, sigv, c=c, s=28, marker=mk, edgecolors='white',
                       linewidths=0.4, label=lb, zorder=5)

    ax.axvline(3.0, color='#e2e8f0', ls='-', lw=0.5, zorder=1)
    ax.axhline(800, color='#e2e8f0', ls='-', lw=0.5, zorder=1)

    bbox_kw = dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85)
    kw = dict(fontsize=5.5, fontstyle='italic', zorder=6)
    ax.text(0.03, 0.97, 'Transport-limited\n(self-sealing risk)', transform=ax.transAxes,
            va='top', color='#dc2626', bbox=bbox_kw, **kw)
    ax.text(0.97, 0.97, 'Distributed\nprecipitation', transform=ax.transAxes,
            va='top', ha='right', color='#16a34a', bbox=bbox_kw, **kw)
    ax.text(0.03, 0.03, 'Clogging risk', transform=ax.transAxes,
            color='#d97706', bbox=bbox_kw, **kw)
    ax.text(0.97, 0.03, 'Reaction-limited\n(low carbonation)', transform=ax.transAxes,
            ha='right', color='#64748b', bbox=bbox_kw, **kw)

    ax.annotate('Increasing \ninjection rate', xy=(0.75, 0.42), xytext=(0.52, 0.68),
                xycoords='axes fraction', textcoords='axes fraction',
                fontsize=5.5, color='#94a3b8', ha='center',
                arrowprops=dict(arrowstyle='->', color='#94a3b8', lw=0.7))

    ax.set_xscale('log')
    ax.set_xlabel('Damköhler number (Da)')
    ax.set_ylabel('Normalised stress ratio (Σ)')

    ax.legend(fontsize=5.5, ncol=6, loc='lower center',
              bbox_to_anchor=(0.5, -0.22), frameon=False,
              columnspacing=0.8, handletextpad=0.3, markerscale=1.0)

    _save(fig,'fig_da_sigma_regime')


# ══════════════════════════════════════════════════════════
#  FIGURE: Gas Saturation Cross-Sections (gravity segregation)
# ══════════════════════════════════════════════════════════
def read_2d_field(sdir, field_names, fallback_shape=(140,25)):
    """Read a 2D x-z cross-section of a field from the last timestep of the last HDF5 file."""
    if not HAS_H5: return None
    h5s = sorted([f for f in glob.glob(os.path.join(sdir,"*.h5"))
                  if 'hanford' not in os.path.basename(f)])
    if not h5s: return None
    try:
        with h5py.File(h5s[-1],'r') as f:
            tgs = sorted([g for g in f.keys() if g.startswith('Time')],
                         key=lambda g: float(g.split(':')[1].strip().split()[0]))
            if not tgs: return None
            g = f[tgs[-1]]
            for fn in field_names:
                if fn in g:
                    arr = g[fn][:]  # shape: (nx, ny, nz)
                    return arr[:,0,:]  # x-z slice (ny=1 for 2D)
    except Exception as e:
        print(f"  warn read_2d_field {sdir}: {e}")
    return None

def fig_gas_saturation_2d():
    print("  Gas saturation 2D cross-sections...")

    # Grid geometry: non-uniform radial
    dx = np.concatenate([np.full(50, 2.0), np.full(40, 10.0), np.full(50, 30.0)])
    x_edges = np.concatenate([[0], np.cumsum(dx)])
    x_centers = 0.5*(x_edges[:-1] + x_edges[1:])

    dz = np.full(25, 4.0)
    z_edges = np.concatenate([[0], np.cumsum(dz)])

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.2), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.29, wspace=0.08, bottom=0.14, right=0.87, top=0.92)

    vmin, vmax = 0, 0.15
    # Custom colormap: light blue (zero gas) → yellow → red (high gas)
    from matplotlib.colors import LinearSegmentedColormap as LSC
    cmap = LSC.from_list('blue_red', ['#e3f2fd','#e8f5e9','#fff9c4','#ffcc80','#e53935'], N=256)
    im = None

    for ax, (dn, lb, c), lab in zip(axes.flat, SC, 'abcdef'):
        sg = read_2d_field(os.path.join(WAG, dn),
                           ['Gas_Saturation', 'Gas Saturation'])
        if sg is not None:
            im = ax.pcolormesh(x_edges, z_edges, sg.T,
                               cmap=cmap, vmin=vmin, vmax=vmax, shading='flat')
        else:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=7, color='#94a3b8')

        ax.set_title(lb, fontsize=7, pad=3)
        _label(ax, lab, y=1.15)

        # Injection interval markers
        ax.axhline(20, color='white', lw=0.6, ls='--', alpha=0.8)
        ax.axhline(80, color='white', lw=0.6, ls='--', alpha=0.8)

        # Zoom to near-wellbore region where features are visible
        ax.set_xlim(0, 500)
        ax.set_ylim(0, 100)

        if lab in 'ad':
            pass  # y-ticks visible from sharey
        if lab in 'def':
            pass  # x-ticks visible from sharex

    # Shared axis labels
    fig.text(0.43, 0.06, 'Distance from well (m)', fontsize=8)
    fig.text(0.07, 0.5, 'Depth (m)', rotation='vertical', fontsize=8)

    # Shared colorbar
    cbar_ax = fig.add_axes([0.89, 0.14, 0.015, 0.78])
    cbar = fig.colorbar(plt.cm.ScalarMappable(
        norm=plt.Normalize(vmin, vmax), cmap=cmap), cax=cbar_ax)
    cbar.set_label('Gas saturation', fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    _save(fig, 'fig_gas_saturation_2d')


# ══════════════════════════════════════════════════════════
#  TABLE: Export data for LaTeX
# ══════════════════════════════════════════════════════════
def export_table():
    print("  Exporting table data...")
    rows=[]
    for dn,lb,c in SC:
        d=read_h5(os.path.join(WAG,dn))
        if d is None: rows.append({'label':lb,'status':'no data'}); continue
        r={'label':lb}
        if len(d.get('porosity_mean',[]))>0:
            r['phi0']=f"{d['porosity_mean'][0]:.5f}"
            r['phi30']=f"{d['porosity_mean'][-1]:.5f}"
        if len(d.get('ph_mean',[]))>0: r['ph']=f"{d['ph_mean'][-1]:.2f}"
        tc=sum(float(d[f'{m}_vf_mean'][-1]) for m in CARB
               if len(d.get(f'{m}_vf_mean',[]))>0 and not np.isnan(d[f'{m}_vf_mean'][-1]))
        r['carb']=f"{tc:.6f}"
        rows.append(r)
    p=os.path.join(OUT,"table_data.json")
    with open(p,'w') as f: json.dump(rows,f,indent=2)
    print(f"  ✔ {p}")


# ══════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════
if __name__=="__main__":
    print("="*50)
    print("  Nature-style Figure Generator")
    print("="*50)
    for d,n in [(WAG,"Baseline"),(DA,"Da sweep"),(KS,"Kinetics")]:
        print(f"  {'✔' if os.path.isdir(d) else '✘'} {n}: {d}")
    print()

    fig_domain()
    fig_kin_sens()
    fig_baseline()
    fig_spatial()
    fig_carb_bar()
    fig_gas_saturation_2d()
    fig_da_sweep()
    fig_da_sigma()
    export_table()

    print(f"\n  All figures → {OUT}/")
