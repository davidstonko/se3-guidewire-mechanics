"""Figure 7 of the paper: the non-planar equilibrium (constant lumen radius).

Reads results/run3d_R6.mat (written by run_demos_3d(30, 6)) and writes
results/fig_nonplanar.png and .pdf.

Usage:  python3 fig_3d_paper.py [path/to/run3d_R6.mat]
Plotting only; every number is computed in MATLAB.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import loadmat

here = os.path.dirname(os.path.abspath(__file__))
matfile = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "results", "run3d_R6.mat")
outdir = os.path.dirname(os.path.abspath(matfile))
F = loadmat(matfile, squeeze_me=True, struct_as_record=False)["fig"]

plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
                     "font.family": "serif", "mathtext.fontset": "stix",
                     "xtick.direction": "out", "ytick.direction": "out",
                     "savefig.dpi": 400, "pdf.fonttype": 42})

P1 = np.atleast_2d(F.P1)
Pv = np.atleast_2d(F.Pv)
Pw = np.atleast_2d(F.Pw)
cn = np.atleast_1d(F.contact).astype(int)
fn = np.atleast_1d(F.fn)
s_w = np.atleast_1d(F.s_w)
h = float(F.h)
strong = cn[fn[cn - 1] > 1e-3]

fig = plt.figure(figsize=(7.0, 6.2))
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.25], hspace=0.55, wspace=0.42)

# A. Oblique view.
ax = fig.add_subplot(gs[0:2, 0], projection="3d")
ax.plot(*P1.T, color="0.55", ls="--", lw=0.9, label="vessel centerline, State 1")
ax.plot(*Pv.T, color="tab:red", lw=1.0, label="vessel centerline, State 2")
ax.plot(*Pw.T, color="k", lw=0.8, label="wire")
ax.scatter(*Pw[strong].T, s=5, color="tab:blue", depthshade=False, label="wire nodes in contact")
ax.set_xlabel("x (mm)", labelpad=-4)
ax.set_ylabel("y (mm)", labelpad=-4)
ax.set_zlabel("z (mm)", labelpad=-6)
ax.tick_params(pad=-2, labelsize=6)
ax.view_init(elev=24, azim=-62)
allp = np.vstack([P1, Pv, Pw])
ax.set_box_aspect(allp.max(0) - allp.min(0))
ax.set_yticks([0, 50, 100])
ax.set_zticks([-60, -30, 0, 30])
pos = ax.get_position()
ax.set_position([pos.x0 - 0.06, pos.y0, pos.width*0.95, pos.height])
ax.legend(loc="upper left", fontsize=6.5, frameon=False, bbox_to_anchor=(-0.05, 1.12))
ax.set_title("(A)", loc="left", fontsize=9, x=-0.05, y=1.12)


def proj(ax, j, lab, tag):
    ax.plot(P1[:, 0], P1[:, j], color="0.55", ls="--", lw=0.9)
    ax.plot(Pv[:, 0], Pv[:, j], color="tab:red", lw=1.0)
    ax.plot(Pw[:, 0], Pw[:, j], color="k", lw=0.8)
    ax.plot(Pw[strong, 0], Pw[strong, j], "o", ms=1.8, color="tab:blue")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel(lab)
    ax.set_title(tag, loc="left", fontsize=9)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


proj(fig.add_subplot(gs[0, 1]), 1, "y (mm)", "(B)")
proj(fig.add_subplot(gs[1, 1]), 2, "z (mm)", "(C)")

# D. Cumulative wall reaction.
ax = fig.add_subplot(gs[2, 0])
ax.step(s_w, F.cum_kkt, where="post", color="k", label="multipliers")
ax.plot(s_w, F.cum_line_vec, color="tab:red", ls="--", label=r"line load with $H=0$")
ax.axvline(float(F.S_il), color="0.8", lw=0.6)
ax.set_xlabel("arc length along the wire (mm)")
ax.set_ylabel("cumulative wall reaction (N)")
ax.legend(fontsize=7, frameon=False, loc="lower right")
ax.set_title("(D)", loc="left", fontsize=9)

# E. Material force.
He = np.atleast_1d(F.He)
hk = np.atleast_1d(F.halfBk2)
sm = (np.arange(len(He)) + 0.5) * h
ax2 = fig.add_subplot(gs[2, 1])
ax2.plot(sm, hk, color="0.55", ls="--", label=r"$\frac{1}{2}B\kappa^2$")
ax2.plot(sm[1:], He[1:], color="k", label=r"material force $H$")
ax2.axhline(0, color="0.8", lw=0.6)
ax2.set_xlabel("arc length along the wire (mm)")
ax2.set_ylabel("force (N)")
ax2.legend(fontsize=7, frameon=False)
ax2.set_title("(E)", loc="left", fontsize=9)
for a in (ax, ax2):
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(outdir, f"fig_nonplanar.{ext}"), bbox_inches="tight")
print("wrote fig_nonplanar (.png, .pdf) to", outdir)
