"""Figures 5 and 6 of the paper, drawn from the MATLAB results (matlab/results/*.mat).
All quantities are computed in MATLAB (postprocess.m). This script only draws them.

Usage: python fig_from_matlab.py [results_dir] [output_dir]
Defaults: ../../matlab/results and ../ relative to this file."""
import os, sys, numpy as np, scipy.io as sio, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__))
res = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, '..', '..', 'matlab', 'results')
out = os.path.join(sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, '..'), '')
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Liberation Serif', 'Times New Roman', 'DejaVu Serif'],
                     'font.size': 9, 'axes.linewidth': 0.6, 'mathtext.fontset': 'stix', 'pdf.fonttype': 42, 'ps.fonttype': 42})
C1, C2, C3, C4 = '#1f4e79', '#c0504d', '#7f7f7f', '#e39d2c'
L = lambda f: sio.loadmat(f'{res}/{f}.mat', squeeze_me=True, variable_names=[
    'P', 'Pv', 'Pw', 'R', 'S_il', 'h', 'B', 'U', 's_v', 's_w', 'contact_nodes', 'cum_kkt', 'cum_line', 'He', 'N_edge', 'halfBk2'])


def rot(Q, ang): c, s = np.cos(ang), np.sin(ang); return Q @ np.array([[c, s], [-s, c]])


def walls(Q, R):   # lumen outline for drawing: offset the centerline by R along the node normal
    e = np.diff(Q, axis=0); tt = np.vstack([e[:1], 0.5 * (e[1:] + e[:-1]), e[-1:]])
    tt /= np.linalg.norm(tt, axis=1)[:, None]; n = np.stack([-tt[:, 1], tt[:, 0]], 1)
    return Q + R[:, None] * n, Q - R[:, None] * n


d = L('base'); P, Pv, Pw, R = d['P'], d['Pv'], d['Pw'], d['R']
T = P[-1] - P[-2]; ang = np.pi / 2 - np.arctan2(T[1], T[0])
P_, Pv_, Pw_ = rot(P, ang), rot(Pv, ang), rot(Pw, ang)
J = np.atleast_1d(d['contact_nodes']).astype(int)          # wire node numbers j (0 = sheath)
fig, ax = plt.subplots(1, 3, figsize=(7.2, 4.2), gridspec_kw={'width_ratios': [0.9, 1.25, 1.25]})
a = ax[0]
for Q, col, lab in [(P_, C3, 'State 1'), (Pv_, C1, 'State 2')]:
    w1, w2 = walls(Q, R); a.plot(w1[:, 0], w1[:, 1], color=col, lw=0.8, label=lab); a.plot(w2[:, 0], w2[:, 1], color=col, lw=0.8)
a.plot(Pw_[:, 0], Pw_[:, 1], color=C2, lw=1.2, label='wire')
a.plot(Pw_[J, 0], Pw_[J, 1], 'o', ms=2.2, color='k', label='contact')
a.annotate('tip', Pw_[-1], xytext=(6, 0), textcoords='offset points', fontsize=8, va='center')
a.annotate('sheath', Pw_[0], xytext=(8, -4), textcoords='offset points', fontsize=8)
a.set_aspect('equal'); a.axis('off'); a.legend(loc='upper center', fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.02), ncol=2)
a.set_title('A', loc='left', fontweight='bold')
a = ax[1]; s = d['s_w']
a.plot(s, d['cum_kkt'], color=C1, lw=1.4, label='KKT multipliers')
a.plot(s, d['cum_line'], color=C2, lw=1.0, ls='--', label=r'$-B(\kappa^{\prime\prime}+\frac{1}{2}\kappa^3)$, $H=0$')
a.axvline(float(d['S_il']), color=C3, lw=0.6, ls=':'); a.text(float(d['S_il']) + 3, 0.05 * d['cum_kkt'][-1], 'end of iliac', fontsize=7, color=C3)
a.set_xlabel('arc length along the wire $s$ (mm)'); a.set_ylabel('cumulative wall reaction (N)')
a.legend(fontsize=7, frameon=False, loc='upper left'); a.set_title('B', loc='left', fontweight='bold')
a = ax[2]; h = float(d['h']); se = (np.arange(len(d['He'])) + 0.5) * h
a.plot(se[1:], -d['halfBk2'][1:], color=C2, lw=1.0, label=r'$-\frac{1}{2}B\kappa^2$')
a.plot(se[1:], d['N_edge'][1:], color=C1, lw=0.9, ls='--', label='axial force $N$ from multipliers')
a.plot(se[1:], d['He'][1:], color='k', lw=0.9, label='material force $H$ from multipliers')
a.axhline(0, color=C3, lw=0.5); a.set_xlabel('arc length along the wire $s$ (mm)'); a.set_ylabel('force (N)')
a.legend(fontsize=7, frameon=False, loc='lower right', bbox_to_anchor=(1.03, 0.02)); a.set_title('C', loc='left', fontweight='bold')
for a in ax[1:]: a.spines[['top', 'right']].set_visible(False)
fig.tight_layout(); fig.savefig(out + 'fig_equilibrium.pdf'); fig.savefig(out + 'fig_equilibrium.png', dpi=200)

fig, a = plt.subplots(1, 1, figsize=(3.9, 2.7))
runs = [('tip60', '60'), ('tip90', '90'), ('base', '120'), ('tip150', '150'), ('tip180', '180')]
sty = ['-', '--', '-.', ':', '-']; cols = [C3, C2, C1, C4, 'k']; lws = [2.2, 1.2, 1.2, 1.2, 0.6]
for (f, lab), ls, c, lw in zip(runs, sty, cols, lws):
    dd = L(f); a.plot(dd['s_v'], dd['U'], ls=ls, color=c, lw=lw, label=lab + ' mm')
a.axvline(float(d['S_il']), color=C3, lw=0.5, ls=':')
a.set_xlabel(r'arc length along the vessel $\sigma$ (mm)'); a.set_ylabel('vessel displacement (mm)')
a.legend(fontsize=6.5, frameon=False, title='tip in the aorta', title_fontsize=6.5)
a.spines[['top', 'right']].set_visible(False)
fig.tight_layout(); fig.savefig(out + 'fig_length.pdf'); fig.savefig(out + 'fig_length.png', dpi=200)
print('ok')
