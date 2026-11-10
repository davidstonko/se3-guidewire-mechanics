"""Figures for the 3D wire-in-vessel computation (run_demos_3d.m).

Reads results/run3d.mat (written by run_demos_3d.m) and writes, in results/:
  fig3d_geometry.png/.pdf          State 1 and State 2 vessel centerlines and the wire,
                                   oblique and top views, contacts marked
  fig3d_reaction.png/.pdf          cumulative wall reaction along the wire:
                                   multipliers against the line load with H = 0
  fig3d_material_force.png/.pdf    material force H along the wire

Usage:  python3 fig_3d.py [path/to/run3d.mat]
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
matfile = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "results", "run3d.mat")
outdir = os.path.dirname(os.path.abspath(matfile))
F = loadmat(matfile, squeeze_me=True, struct_as_record=False)["fig"]

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
                     "xtick.direction": "out", "ytick.direction": "out",
                     "savefig.dpi": 300, "pdf.fonttype": 42})

P1 = np.atleast_2d(F.P1)          # State 1 vessel centerline (N x 3)
Pv = np.atleast_2d(F.Pv)          # State 2 vessel centerline
Pw = np.atleast_2d(F.Pw)          # wire nodes 0..M-1
cn = np.atleast_1d(F.contact).astype(int)   # wire node numbers j >= 1 in contact
fn = np.atleast_1d(F.fn)          # contact force at nodes 1..M-1
s_w = np.atleast_1d(F.s_w)
h = float(F.h)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)


# 1. Geometry: oblique and top views.
fig = plt.figure(figsize=(7.0, 3.4))
ax = fig.add_subplot(1, 2, 1, projection="3d")
ax.plot(*P1.T, color="0.6", ls="--", lw=0.9, label="vessel, State 1")
ax.plot(*Pv.T, color="tab:red", lw=1.0, label="vessel, State 2")
ax.plot(*Pw.T, color="k", lw=0.8, label="wire")
strong = cn[fn[cn - 1] > 1e-3]
ax.scatter(*Pw[strong].T, s=6, color="tab:blue", depthshade=False, label="contact")
ax.set_xlabel("x (mm)", labelpad=-2)
ax.set_ylabel("y (mm)", labelpad=-2)
ax.set_zlabel("z (mm)", labelpad=-2)
ax.tick_params(pad=-1, labelsize=7)
ax.view_init(elev=22, azim=-60)
# equal axes
mins = np.vstack([P1, Pv, Pw]).min(0)
maxs = np.vstack([P1, Pv, Pw]).max(0)
ax.set_box_aspect(maxs - mins)
ax.legend(loc="upper left", fontsize=7, frameon=False)

ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(P1[:, 0], P1[:, 1], color="0.6", ls="--", lw=0.9)
ax2.plot(Pv[:, 0], Pv[:, 1], color="tab:red", lw=1.0)
ax2.plot(Pw[:, 0], Pw[:, 1], color="k", lw=0.8)
ax2.plot(Pw[strong, 0], Pw[strong, 1], "o", ms=2.0, color="tab:blue")
ax2.set_aspect("equal")
ax2.set_xlabel("x (mm)")
ax2.set_ylabel("y (mm)")
for s in ("top", "right"):
    ax2.spines[s].set_visible(False)
fig.tight_layout()
save(fig, "fig3d_geometry")

# 2. Cumulative wall reaction along the wire.
fig, ax = plt.subplots(figsize=(4.6, 3.0))
ax.step(s_w, F.cum_kkt, where="post", color="k", label="multipliers, all contact nodes")
ax.plot(s_w, F.cum_line_vec, color="tab:red", ls="--", label=r"line load $B y'''' + \frac{3}{2}B(\kappa^2 y')'$")
ax.step(s_w, F.cum_kkt_frenet_nodes, where="post", color="0.55", lw=0.9,
        label=r"multipliers, nodes with $\kappa > 2\times10^{-3}$ mm$^{-1}$")
ax.plot(s_w, F.cum_line_frenet, color="tab:blue", ls=":", lw=1.2,
        label=r"$|B(\kappa'' + \kappa^3/2 - \kappa\tau^2)|$, same nodes")
ax.axvline(float(F.S_il), color="0.8", lw=0.6)
ax.set_xlabel("arc length along the wire (mm)")
ax.set_ylabel("cumulative wall reaction (N)")
ax.legend(fontsize=7, frameon=False, loc="upper left")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
save(fig, "fig3d_reaction")

# 3. Material force H on each edge.
He = np.atleast_1d(F.He)
hk = np.atleast_1d(F.halfBk2)
He0 = np.atleast_1d(F.He_planar)
sm = (np.arange(len(He)) + 0.5) * h
sm0 = (np.arange(len(He0)) + 0.5) * h
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8))
ax = axs[0]
ax.plot(sm, hk, color="0.6", ls="--", label=r"$B\kappa^2/2$")
ax.plot(sm[1:], He[1:], color="k", label=r"$H = 2h\nu$")
ax.set_xlabel("arc length along the wire (mm)")
ax.set_ylabel("force (N)")
ax.legend(fontsize=7, frameon=False)
ax = axs[1]
ax.plot(sm[1:], He[1:], color="k", label="non-planar vessel")
ax.plot(sm0[1:], He0[1:], color="tab:red", ls="--", label="planar vessel")
ax.axhline(0, color="0.8", lw=0.6)
ax.set_xlabel("arc length along the wire (mm)")
ax.set_ylabel(r"$H$ (N), first edge omitted")
ax.legend(fontsize=7, frameon=False)
for a in axs:
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
fig.tight_layout()
save(fig, "fig3d_material_force")
print("wrote fig3d_geometry, fig3d_reaction, fig3d_material_force (.png, .pdf) to", outdir)
