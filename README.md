# se3-guidewire-mechanics

MATLAB code for the computations in

> Stonko DP. *Solving the Wire Stiffening Problem: Guidewire-vessel mechanics on SE(3).* Manuscript submitted to Proceedings of the Royal Society A, 2026.

A stiff guidewire is clamped at the sheath in the groin, its tip floats free in the descending thoracic aorta, and it must stay inside the lumen of an elastic vessel held by surrounding tissue. The code computes the equilibrium of wire and vessel on a synthetic planar test vessel and checks it against the analytical results of the paper. It reproduces every number in Section 4 and the data in Figures 4 to 6.

All computation is in MATLAB. Python is used only to draw figures from the saved MATLAB results.

## Requirements

- MATLAB R2025b with the Optimization Toolbox (`fmincon`). Earlier releases with the interior-point `fmincon` and `HessianFcn` should work but were not tested.
- For the figures only: Python 3 with NumPy, SciPy and Matplotlib.

## Reproduce the paper

```matlab
cd matlab
run_demos
```

This takes under a minute on an Apple M1 laptop. It

1. checks every analytic derivative against central finite differences,
2. solves the base case (tip 120 mm into the aorta) by continuation from a nearly rigid vessel,
3. solves the same problem with the tip at 90, 60, 30, 0, 150 and 180 mm into the aorta,
4. writes `results/*.mat` and `results/summary.json`,
5. runs `audit_solution` on every saved equilibrium.

Then, for the figures:

```bash
cd figures/python
python fig_obstacle.py      # Figure 4 (panel B from results/base.mat)
python fig_from_matlab.py   # Figures 5 and 6
```

The output of the run used for the paper is in [`matlab/results/run_demos.log`](matlab/results/run_demos.log).

## Where each number in the paper comes from

| Paper | Value | Source |
|---|---|---|
| Section 2.4, Appendix B: derivative check | relative error below 1e-8 | `check_derivatives`, first lines of the log |
| Section 4: energy of the base case | 114.1455 N mm | `summary.json`, `base.E` |
| Section 4: peak and RMS displacement | 5.73 mm, 1.86 mm | `base.peak`, `base.rms` |
| Section 4: contact arcs and their forces | 2.01, 3.08, 1.12 N | `base.arcs` |
| Section 4: wall reaction, multipliers against line load with H = 0 | 6.21 N against 6.22 N | `base.F_kkt`, `base.F_line` |
| Section 4: material force after the clamped edge | below 0.02 N | `base.H_absmax_after_first_edge` |
| Section 4: compression at the sheath exit | 0.94 N | `base.halfBk2_exit` |
| Section 4: rigid vessel, total wall reaction | 11.3 N | `rigid.F_kkt`, `results/rigid.mat` |
| Section 4: rigid vessel, contacts | 4.32 and 4.57 N at the inner walls, 2.13 and 0.25 N in the aorta | `rigid.arcs`, rigid-vessel lines of the log |
| Section 4, Figure 6: other wire lengths | E = 114.1455 (90, 150, 180 mm), 114.1432 and 0.20 N tip force (60 mm) | `length` entries |
| Section 4: tip on the taper | not converged (exit flag 2) | `length` entries for 30 and 0 mm |

## Files

| File | Purpose |
|---|---|
| `matlab/run_demos.m` | Every computation reported in Section 4 |
| `matlab/evar_geometry.m` | State 1 test vessel: sine iliac segment and straight aorta, lumen radius 6 to 10 mm |
| `matlab/build_problem.m` | Problem data, unknowns, and the column of every node coordinate |
| `matlab/pack_nodes.m`, `unpack_nodes.m` | Node positions to and from the vector of unknowns |
| `matlab/bend_energy.m` | Turning-angle bending energy of a polygon, with exact gradient and Hessian |
| `matlab/objective.m` | Total energy: wire, vessel and foundation |
| `matlab/contact_terms.m` | Lumen constraint at the wire nodes, with exact Jacobian and weighted Hessian |
| `matlab/constraints.m` | All constraints in the form `fmincon` expects |
| `matlab/lagrangian_hessian.m` | Exact Hessian of the Lagrangian |
| `matlab/solve_equilibrium.m` | One interior-point solve |
| `matlab/continuation_solve.m` | Nearly rigid vessel, then raise the wire stiffness, then relax the foundation |
| `matlab/postprocess.m` | Contact forces, line load with H = 0, material force H from the multipliers, contact arcs |
| `matlab/check_derivatives.m` | Finite-difference check of every derivative |
| `matlab/audit_solution.m` | Independent checks of a saved equilibrium: feasibility, inactive box bounds, the signs required by Corollary 3.4.1 on each contact arc, and the force balance |
| `matlab/results/` | Saved equilibria, `summary.json` and the run log |
| `figures/python/` | Drawing scripts for Figures 4 to 6 |
| `figures/*.pdf` | Figures 1 to 6 as they appear in the paper. Figures 1 to 3 are schematics. |

## The discrete problem

    minimize   B/2 sum phi_j^2/l_j  +  B_v/2 sum (phi_i - phi_i^ref)^2/l_i  +  1/2 sum k_i |x_i - P_i|^2
    subject to dist(y_j, vessel polygon)^2 - delta_j^2 <= 0      (wire inside the lumen)
               |y_j - y_{j-1}|^2 = h^2                           (wire does not stretch)
               (y_1 - y_0) x t_0 = 0                              (direction of the sheath)

`phi` is the signed turning angle at a node and `l` the Voronoi length. `x_i` are the vessel nodes, `y_j` the wire nodes and `P_i` State 1. By Proposition 2.1.3 of the paper the turning-angle energy is the quadratic strain energy of the discrete Bishop frame.

Parameters: wire bending stiffness B = 158.4 GPa times the second moment of area of a 0.035 inch section (4856.6 N mm^2, Lunderquist Extra Stiff, Harrison et al. 2011), vessel bending stiffness 200 N mm^2, foundation 0.0158 N/mm^2, node spacing 1.598 mm (the spacing closest to 1.6 mm that divides the vessel evenly).

## Notes from the code audit (31 October 2026)

- `run_demos` was re-run from a clean copy. Every reported number matched the 16 October run to the printed precision.
- `audit_solution` confirms, for each saved equilibrium, that the lumen constraint and edge lengths hold to 2e-7 mm^2 and 1e-13 mm, that the 25 mm box bounds are inactive (largest node move 21.0 mm), and that on both inner-wall arcs the curvature toward the wall is positive and concave and the line load is positive, as Corollary 3.4.1 requires.
- In `continuation_solve`, an intermediate step that stops on step size (exit flag 2) is accepted and the continuation goes on. Only solves with exit flag 1 (first-order optimality met) are reported.
- When `run_demos` is driven from outside the MATLAB desktop (for example by an automation tool), `diary` may write an empty log. The log in `results/` was captured from the command window.

## License

MIT. See [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff). Please cite the paper once it is published.
