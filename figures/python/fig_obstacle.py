"""Figure: the equilibrium as an obstacle problem with a deformable obstacle.
Panel A is a schematic of the classical problem. Panel B draws the MATLAB equilibrium (matlab/results/base.mat);
contact nodes and forces are read from that file, not recomputed.

Usage: python fig_obstacle.py [output_dir] [results_dir]"""
import os, sys, numpy as np, matplotlib, scipy.io as sio
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, '..'), '')
res = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, '..', '..', 'matlab', 'results')
plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'serif',
                     'font.serif': ['Liberation Serif', 'Times New Roman', 'DejaVu Serif'],
                     'mathtext.fontset': 'stix', 'font.size': 8})
C1, C2, CG, CY, CB = '#1f4e79', '#c0392b', '0.55', '#e0a800', '#8a5a00'

fig = plt.figure(figsize=(7.2, 3.9))
# ---------- A: string pinned at both ends, pushed up by a rigid obstacle ----------
c = 1.0
for k, (h, lab) in enumerate([(0.08, 'barely touching'), (0.35, 'pushed higher'), (0.8, 'higher still')]):
    ax = fig.add_axes([0.02 + 0.33 * k, 0.56, 0.30, 0.38]); ax.axis('off')
    a = 1 - np.sqrt(1 - h / c)                       # free boundary, from tangency
    xs = np.linspace(-1.05, 1.05, 400); psi = h - c * xs ** 2
    ax.plot(xs, psi, color=C2, lw=0.9)
    xa = np.linspace(-a, a, 50)
    ax.plot([-1, -a], [0, h - c * a * a], color=C1, lw=1.4); ax.plot([a, 1], [h - c * a * a, 0], color=C1, lw=1.4)
    ax.plot(xa, h - c * xa ** 2, color=CY, lw=3.0, solid_capstyle='butt')
    ax.plot([-1, 1], [0, 0], 'o', color='k', ms=3)
    ax.plot([-a, a], [h - c * a * a] * 2, 'o', color=CB, ms=3)
    for xi in np.linspace(-a, a, 3 if a > 0.15 else 1):
        yi = h - c * xi ** 2
        ax.add_patch(FancyArrowPatch((xi, yi - 0.22), (xi, yi - 0.01), arrowstyle='-|>', mutation_scale=6, color=C2, lw=0.8))
    ax.plot([-a, -a], [h - c * a * a, -0.32], color=CB, lw=0.5, ls=':'); ax.plot([a, a], [h - c * a * a, -0.32], color=CB, lw=0.5, ls=':')
    ax.plot([-a, a], [-0.32, -0.32], color=CB, lw=0.8)
    ax.text(0, -0.43, 'contact set', ha='center', va='top', fontsize=7, color=CB)
    ax.text(0, 1.0, lab, ha='center', fontsize=7.5, color='0.25')
    if k == 0:
        ax.text(-1.0, 0.1, r'string $u$', fontsize=7.5, color=C1)
        ax.text(0.28, 0.22, r'obstacle $\psi$', fontsize=7.5, color=C2)
    if k == 2:
        ax.annotate('free boundary', (a, h - c * a * a), xytext=(0.55, 0.75), fontsize=7, color=CB,
                    arrowprops=dict(arrowstyle='-', lw=0.5, color=CB))
        ax.text(0.08, 0.38, r'$\lambda$', fontsize=8, color=C2)
    ax.set_xlim(-1.12, 1.12); ax.set_ylim(-0.6, 1.08)
fig.text(0.005, 0.95, 'A', fontweight='bold', fontsize=10)

# ---------- B: the stiff wire in an elastic lumen (computed) ----------
d = sio.loadmat(res + '/base.mat', squeeze_me=True, variable_names=['P', 'Pv', 'Pw', 'R', 'fn', 'S_il', 'contact_nodes'])
P, Pv, Pw, R, fn = d['P'], d['Pv'], d['Pw'], d['R'], d['fn']
hv = np.linalg.norm(P[1] - P[0]); S = float(d['S_il'])


def walls(Q, Rn):
    e = np.gradient(Q, axis=0); t = e / np.linalg.norm(e, axis=1, keepdims=True); n = np.stack([-t[:, 1], t[:, 0]], 1)
    return Q + Rn[:, None] * n, Q - Rn[:, None] * n


nv = int(S / hv) + 1                       # iliac segment only
bx = fig.add_axes([0.03, 0.0, 0.94, 0.54]); bx.set_aspect('equal'); bx.axis('off')
for Q, col, ls, lw in [(P, CG, (0, (3, 2)), 0.7), (Pv, C1, '-', 1.0)]:
    w1, w2 = walls(Q[:nv], R[:nv])
    bx.plot(*w1.T, color=col, ls=ls, lw=lw); bx.plot(*w2.T, color=col, ls=ls, lw=lw)
mw = np.nonzero(Pw[:, 0] <= P[nv - 1, 0])[0]
bx.plot(*Pw[mw].T, color='k', lw=1.4)
act = np.atleast_1d(d['contact_nodes']).astype(int)
act = act[Pw[act, 0] <= P[nv - 1, 0]]
# contact arcs, highlighted on the wire
groups = np.split(act, np.nonzero(np.diff(act) > 1)[0] + 1)
for gr in groups:
    if len(gr) == 0:
        continue
    bx.plot(*Pw[gr[0]:gr[-1] + 1].T, color=CY, lw=3.2, solid_capstyle='butt', zorder=4)
    # the force of the wire on the wall: outward from the vessel axis, total over the arc
    jm = gr[len(gr) // 2]
    dd = np.linalg.norm(Pv - Pw[jm], axis=1); iv = int(np.argmin(dd)); nh = (Pw[jm] - Pv[iv]) / dd[iv]
    F = fn[gr - 1].sum()
    L = 4.0 + 3.0 * F
    bx.add_patch(FancyArrowPatch(tuple(Pw[jm] + 0.5 * nh), tuple(Pw[jm] + L * nh), arrowstyle='-|>', mutation_scale=7,
                                 color=C2, lw=1.0, zorder=5))
    bx.text(*(Pw[jm] + (L + 2.2) * nh), f'{F:.1f} N', fontsize=7, color=C2, ha='center', va='center')
bx.text(P[0, 0] - 2, P[0, 1] - 9, 'sheath', fontsize=7, ha='center')
bx.text(P[nv - 1, 0] + 2, P[nv - 1, 1] - 2, 'to the aorta', fontsize=7)
bx.text(60, 34, r'lumen wall, State 1 (dashed) and State 2', fontsize=7.5, color=C1, ha='center')
jm2 = groups[-1][len(groups[-1]) // 2]
bx.annotate('contact arc', Pw[jm2], xytext=(150, -30), fontsize=7.5, color=CB, arrowprops=dict(arrowstyle='-', lw=0.5, color=CB))
jm1 = groups[0][len(groups[0]) // 2]
bx.annotate('contact arc', Pw[jm1], xytext=(5, 30), fontsize=7.5, color=CB, arrowprops=dict(arrowstyle='-', lw=0.5, color=CB))
iw = int(np.argmin(np.abs(Pw[:, 0] - 95)))
bx.annotate(r'stiff wire $p_w(s)$', Pw[iw], xytext=(70, -30), fontsize=7.5, arrowprops=dict(arrowstyle='-', lw=0.5, color='k'))
fig.text(0.005, 0.50, 'B', fontweight='bold', fontsize=10)
fig.savefig(out + 'fig_obstacle.pdf'); fig.savefig(out + 'fig_obstacle.png', dpi=300)
print('ok', [len(gr) for gr in groups])
