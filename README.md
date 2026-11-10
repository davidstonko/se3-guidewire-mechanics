# se3-guidewire-mechanics

MATLAB code for the computations in

> Stonko DP. *Solving the Wire Stiffening Problem: Guidewire-vessel mechanics on SE(3).* Manuscript submitted to Proceedings of the Royal Society A, 2026.

A stiff guidewire is clamped at the sheath in the groin, its tip floats free in the descending thoracic aorta, and it must stay inside the lumen of an elastic vessel held by surrounding tissue. The code computes the equilibrium of wire and vessel on a synthetic planar test vessel and checks it against the analytical results of the paper. It reproduces every number in Section 7 and the data in Figures 4 to 7. Sections 7.1 to 7.4 use the planar code in `matlab/`. Sections 7.5 to 7.7 (constant lumen radius, mesh refinement, a non-planar vessel and the clock position of a wall target) use `matlab/run_checks_v17.m` and the spatial code in `matlab/threed/`.

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
| Section 7.1, Appendix B: derivative check | relative error below 1e-8 | `check_derivatives`, first lines of the log |
| Section 7: energy of the base case | 114.1455 N mm | `summary.json`, `base.E` |
| Section 7: peak and RMS displacement | 5.73 mm, 1.86 mm | `base.peak`, `base.rms` |
| Section 7: contact arcs and their forces | 2.01, 3.08, 1.12 N | `base.arcs` |
| Section 7: wall reaction, multipliers against line load with H = 0 | 6.21 N against 6.22 N | `base.F_kkt`, `base.F_line` |
| Section 7: material force after the clamped edge | below 0.02 N | `base.H_absmax_after_first_edge` |
| Section 7: compression at the sheath exit | 0.94 N | `base.halfBk2_exit` |
| Section 7: rigid vessel, total wall reaction | 11.3 N | `rigid.F_kkt`, `results/rigid.mat` |
| Section 7: rigid vessel, contacts | 4.32 and 4.57 N at the inner walls, 2.13 and 0.25 N in the aorta | `rigid.arcs`, rigid-vessel lines of the log |
| Section 7: constraint qualification (LICQ) at each converged solution | smallest singular value of the row-scaled active Jacobian 0.012 to 0.017 | `kkt_diagnostics`, log |
| Section 7: second-order sufficient condition | reduced Hessian positive definite, smallest eigenvalue 3.8e-3 N/mm (base case), 3.8e-4 to 5.0e-2 N/mm (all cases) | `kkt_diagnostics`, log |
| Section 7: discrete thickness of the deformed vessel | min local radius 18.1 mm, min half-chord 14.6 mm, against delta at most 9.56 mm | `kkt_diagnostics`, log |
| Section 7: conditions of Proposition 5.3 | contacts project inside the vessel, no contacts on opposite walls of one cross-section | `kkt_diagnostics`, log |
| Section 7.4, Figure 6: other wire lengths | E = 114.1455 (90, 150, 180 mm), 114.1432 and 0.20 N tip force (60 mm) | `length` entries |
| Section 7: tip on the taper | not converged (exit flag 2) | `length` entries for 30 and 0 mm |

## Sections 7.5 to 7.7

```matlab
cd matlab
run_checks_v17('kkt')              % strict complementarity, nearest segment, second order in H^1
run_checks_v17('const')            % the planar problem with a constant lumen radius of 6 mm
run_checks_v17('nested', [], 6)    % constant radius at h = 1.6, 0.8 and 0.4 mm (Table 5), about 4 minutes
cd threed
run_demos_3d(30, 6)                % non-planar vessel, constant radius (Section 7.6, Table 6, Table 7)
run_demos_3d()                     % the same with the widening lumen
soft_mode3d(30, 6)                 % the softest admissible mode of the non-planar equilibrium
```

Then `python3 fig_3d_paper.py` in `matlab/threed/` draws Figure 7 from `results/run3d_R6.mat`.

| Paper | Value | Source |
|---|---|---|
| Table 5: wall reaction at h = 1.6, 0.8, 0.4 mm | 6.486, 6.538, 6.562 N | `results/checks_nested_R6.log` |
| Table 5: line load with H = 0 | 6.478, 6.535, 6.562 N | same |
| Table 5: largest material force after the first edge | 0.0191, 0.0066, 0.0019 N | same |
| Section 7.5: smallest H^1 eigenvalue of the reduced Hessian | 1.07e-3, 1.22e-3, 4.98e-4 N/mm | same |
| Section 7.5: widening lumen at h = 0.8 mm | no first-order certificate | `results/checks_nested.log`, `results/batch_2026-11-10_console.txt` |
| Section 7.6: spatial code on the planar geometry | wall reaction 6.472 N, peak 5.871 mm, rigid 10.485 N | `threed/results/run3d_R6.log`, part 3 |
| Section 7.6, Table 6: non-planar vessel | wall reaction 5.879 N (rigid 8.283 N), peak 4.901 mm | same, part 4 |
| Section 7.6: material force | below 0.0155 N, against 0.967 N for B kappa^2/2 | same |
| Section 7.6: line load, vector form, against multipliers | 5.8765 against 5.8792 N | same |
| Section 7.6: softest mode | 4.65e-4 N/mm, 99.7% in the aorta, 88% around the lumen axis | `threed/results/soft_mode3d_A30_R6.log` |
| Section 7.7, Table 7: clock shift | +2.81, -0.48, -6.86, -1.91, -0.24, -0.89 degrees | `threed/results/run3d_R6.log`, part 5 |

Notes on these runs (10 November 2026):

- The planar code uses a geometrically exact vessel (turning-angle energy). The spatial code uses the geometrically linear beam, `B_v/2 sum |u_{i+1} - 2u_i + u_{i-1}|^2/h^3`, and the curvature-binormal wire energy of Bergou et al. The comparison on the planar geometry therefore compares two vessel models as well as two codes.
- `run_checks_v17('nested')` starts each finer mesh from the interpolated solution of the coarser one and accepts only solves that end with the first-order certificate (exit flag 1). With a constant radius this works at 0.8 and 0.4 mm. With the widening lumen the solve at 0.8 mm stops on the step tolerance, by continuation (`continuation_solve_strict`, `mesh_strict`) and from the coarser solution, and that problem is reported at 1.6 mm only.
- Calls that take longer than about a minute should be run with `matlab -batch`.

## Files

| File | Purpose |
|---|---|
| `matlab/run_demos.m` | Every computation reported in Section 7 |
| `matlab/kkt_diagnostics.m` | LICQ, second-order check, discrete thickness and the conditions of Proposition 5.3 at a solution |
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
| `matlab/audit_solution.m` | Independent checks of a saved equilibrium: feasibility, inactive box bounds, the signs required by Corollary 5.6 on each contact arc, and the force balance |
| `matlab/run_checks_v17.m` | Further checks of the planar computation: strict complementarity, the second-order condition in H^1, constant radius, mesh refinement |
| `matlab/continuation_solve_strict.m`, `mesh_strict.m` | Continuation that accepts only exit flag 1, used to document the stall at h = 0.8 mm with the widening lumen |
| `matlab/run_logged.m` | Run a function and keep the command-window output in a file |
| `matlab/threed/run_demos_3d.m` | Every spatial computation of Sections 7.6 and 7.7 |
| `matlab/threed/geometry3d.m`, `build_problem3d.m`, `objective3d.m`, `wire_bend3d.m`, `contact_terms3d.m`, `constraints3d.m`, `lagrangian_*3d.m`, `solve_equilibrium3d.m`, `continuation_solve3d.m` | The spatial counterparts of the planar files |
| `matlab/threed/postprocess3d.m`, `kkt_diagnostics3d.m`, `clock_shift3d.m`, `soft_mode3d.m`, `check_derivatives3d.m` | Forces, line load in vector and Frenet form, KKT diagnostics, clock shift (Proposition 2.4), softest mode, derivative check |
| `matlab/threed/fig_3d_paper.py`, `fig_3d.py` | Figure 7, and working figures |
| `matlab/results/` | Saved equilibria, `summary.json` and the run log |
| `figures/python/` | Drawing scripts for Figures 4 to 6 |
| `figures/*.pdf` | Figures 1 to 7 as they appear in the paper. Figures 1 to 3 are schematics. |

## The discrete problem

    minimize   B/2 sum phi_j^2/l_j  +  B_v/2 sum (phi_i - phi_i^ref)^2/l_i  +  1/2 sum k_i |x_i - P_i|^2
    subject to dist(y_j, vessel polygon)^2 - delta_j^2 <= 0      (wire inside the lumen)
               |y_j - y_{j-1}|^2 = h^2                           (wire does not stretch)
               (y_1 - y_0) x t_0 = 0                              (direction of the sheath)

`phi` is the signed turning angle at a node and `l` the Voronoi length. `x_i` are the vessel nodes, `y_j` the wire nodes and `P_i` State 1. By Lemma 7.1 of the paper the turning-angle energy is the quadratic strain energy of the discrete Bishop frame.

Parameters: wire bending stiffness B = 158.4 GPa times the second moment of area of a 0.035 inch section (4856.6 N mm^2, Lunderquist Extra Stiff, Harrison et al. 2011), vessel bending stiffness 200 N mm^2, foundation 0.0158 N/mm^2, node spacing 1.598 mm (the spacing closest to 1.6 mm that divides the vessel evenly).

## Notes from the code audit (31 October 2026)

- `run_demos` was re-run from a clean copy. Every reported number matched the 16 October run to the printed precision.
- `audit_solution` confirms, for each saved equilibrium, that the lumen constraint and edge lengths hold to 2e-7 mm^2 and 1e-13 mm, that the 25 mm box bounds are inactive (largest node move 21.0 mm), and that on both inner-wall arcs the curvature toward the wall is positive and concave and the line load is positive, as Corollary 5.6 requires.
- In `continuation_solve`, an intermediate step that stops on step size (exit flag 2) is accepted and the continuation goes on. Only solves with exit flag 1 (first-order optimality met) are reported.
- When `run_demos` is driven from outside the MATLAB desktop (for example by an automation tool), `diary` may write an empty log. The log in `results/` was captured from the command window.

## License

MIT. See [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff). Please cite the paper once it is published.
