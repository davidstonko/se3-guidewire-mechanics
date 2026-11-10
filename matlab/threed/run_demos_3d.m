function S = run_demos_3d(A, Rconst)
%RUN_DEMOS_3D  3D version of the planar wire-in-vessel computation.
%
%   S = RUN_DEMOS_3D() runs, with the parameters of the planar RUN_DEMOS
%   (node spacing about 1.6 mm, B = 158.4 GPa * pi d^4/64 with d = 0.889 mm,
%   B_v = 200 N mm^2, foundation 0.0158 N/mm^2, tip 120 mm into the aorta):
%
%     1. the non-planar State 1 geometry of GEOMETRY3D (amplitude A, default
%        30 mm) and its curve checks, with a scan over A;
%     2. a derivative check of every analytic gradient and Hessian;
%     3. validation on the planar geometry (A = 0): rigid vessel and
%        deformable vessel, compared with the planar results in
%        ../results/summary.json (read only);
%     4. the non-planar base case: energy, wall reaction, contact arcs,
%        displacement, material force, line-load check, KKT diagnostics;
%     5. the clock shift of twist-free wall targets, for both geometries;
%     6. timings.
%
%   The command window output is also written to results/run3d.log, the
%   results to results/run3d.mat and results/run3d_summary.json. The function
%   changes no global state of the MATLAB session (no cd, path, diary,
%   random seed or timer): output is captured with evalc.
%
%   Figures: python3 fig_3d.py (reads results/run3d.mat).

if nargin < 1 || isempty(A), A = 30; end
if nargin < 2, Rconst = []; end
tag = 'run3d';  if ~isempty(Rconst), tag = sprintf('run3d_R%g', Rconst); end
here = fileparts(mfilename('fullpath'));
resdir = fullfile(here, 'results');
if ~exist(resdir, 'dir'), mkdir(resdir); end
S = struct();
txt = evalc('S = run_core(here, resdir, A, Rconst, tag);');
fprintf('%s', txt);
fid = fopen(fullfile(resdir, [tag '.log']), 'w');  fprintf(fid, '%s', txt);  fclose(fid);
if nargout == 0, clear S; end            % do not leave ans in the caller's workspace
end

function S = run_core(here, resdir, A, Rconst, tag)
S = struct();
try
    S = run_all(here, resdir, A, Rconst, tag);
catch ME
    fprintf('\nERROR: %s\n', ME.message);
    for q = 1:numel(ME.stack)
        fprintf('  in %s at line %d\n', ME.stack(q).name, ME.stack(q).line);
    end
end
end

function S = run_all(here, resdir, A, Rconst, tag)
t_all = tic;
fprintf('%s   MATLAB %s\n', char(datetime('now')), version);
if ~isempty(Rconst)
    fprintf('CONSTANT LUMEN RADIUS %g mm on the whole vessel. The planar-code numbers quoted for comparison are from ../results/summary.json, which has the widening lumen.\n', Rconst);
end
h = 1.6;  kbar = 0.0158;
B = 158.4e3*pi*0.889^4/64;               % N mm^2
a_ant = [0 0 1];                         % "anterior" vector for the clock frame
planar = jsondecode(fileread(fullfile(here, '..', 'results', 'summary.json')));

%% 1. Geometry.
fprintf('\n1. Geometry: iliac (x, 25 sin(2 pi x/150), A sin(pi x/150)), 0 <= x <= 150 mm, then a straight 200 mm aorta\n');
fprintf('   A (mm)  S_il (mm)  min R_c (mm)  max|tau| (1/mm)  RMS tau  polygon: min local radius, min dist of nodes > 10 pi apart (mm)\n');
for Ai = [0 10 20 30 40 50]
    Gi = geometry3d(h, Ai, Rconst);  c = Gi.check;
    fprintf('   %5.0f   %8.2f   %10.2f   %14.4f   %7.4f   %6.2f   %6.2f\n', Ai, Gi.S_il, c.min_Rc, c.max_abs_tau, c.rms_tau, ...
        c.poly_min_local_radius, c.poly_min_far_dist);
end
G = geometry3d(h, A, Rconst);  c = G.check;
fprintf('   chosen A = %g mm: iliac arc length %.2f mm, vessel length %.2f mm, %d nodes, spacing %.4f mm\n', ...
    A, G.S_il, G.S_il + 200, size(G.P,1), G.hv);
fprintf('   min radius of curvature %.2f mm (required >= 18), max |torsion| %.4f 1/mm (where kappa > 1e-3), RMS torsion %.4f 1/mm\n', ...
    c.min_Rc, c.max_abs_tau, c.rms_tau);
fprintf('   tube of radius 10 mm: min local radius of the node polygon %.2f mm (> 10), min distance of nodes more than 10 pi mm apart %.2f mm (> 20): %s\n', ...
    c.poly_min_local_radius, c.poly_min_far_dist, tf(c.tube_ok));
Tao = (G.P(end,:) - G.P(end-1,:))/norm(G.P(end,:) - G.P(end-1,:));
fprintf('   out-of-plane extent max z = %.2f mm; aortic axis (end tangent of the iliac segment) (%.4f, %.4f, %.4f)\n', ...
    max(G.P(:,3)), Tao);

%% 2. Derivative check.
fprintf('\n2. Derivative check (coarse mesh h = 6 mm, A = %g mm, random perturbation)\n', A);
S.deriv = check_derivatives3d(A);

%% 3. Planar validation.
fprintf('\n3. Validation on the planar geometry (A = 0)\n');
G0 = geometry3d(h, 0, Rconst);
t3 = tic;
[x0, out0, prob0, steps0, rig0] = continuation_solve3d(G0, G0.S_il + 120, B, kbar);
time_planar = toc(t3);
P0 = postprocess3d(prob0, x0, out0);
R0 = postprocess3d(rig0.prob, rig0.x, rig0.out);
fprintf('  wire nodes %d, vessel nodes %d, unknowns %d (rigid %d)\n', prob0.M, prob0.N, prob0.n, rig0.prob.n);
fprintf('  planarity of the 3D solution: max |z| wire %.1e mm, vessel %.1e mm\n', max(abs(P0.Pw(:,3))), max(abs(P0.Pv(:,3))));
fprintf('  rigid vessel (B_v -> infinity): energy %.6f (planar %.6f), wall reaction %.3f N (planar %.3f N, difference %+.2f%%)\n', ...
    R0.E, planar.rigid.E, R0.F_kkt, planar.rigid.F_kkt, 100*(R0.F_kkt/planar.rigid.F_kkt - 1));
print_arcs_compare(R0.arcs, planar.rigid.arcs);
fprintf('  deformable vessel: energy %.6f (planar %.6f), wall reaction %.3f N (planar %.3f N, difference %+.2f%%)\n', ...
    P0.E, planar.base.E, P0.F_kkt, planar.base.F_kkt, 100*(P0.F_kkt/planar.base.F_kkt - 1));
fprintf('    energy parts: wire %.4f (planar %.4f), vessel bending %.4f (planar %.4f), foundation %.4f (planar %.4f)\n', ...
    P0.E_wire, planar.base.E_wire, P0.E_vessel, planar.base.E_vessel, P0.E_foundation, planar.base.E_foundation);
print_arcs_compare(P0.arcs, planar.base.arcs);
fprintf('    peak displacement %.4f mm at s_v = %.1f mm (planar %.4f mm at %.1f mm, difference %+.2f%%), first bend %.2f (planar %.2f), aorta %.2f (planar %.2f), RMS %.2f (planar %.2f) mm\n', ...
    P0.peak, P0.peak_at_s, planar.base.peak, planar.base.peak_at_s, 100*(P0.peak/planar.base.peak - 1), ...
    P0.peak_first_bend, planar.base.peak_first_bend, P0.aorta_shift, planar.base.aorta_shift, P0.rms, planar.base.rms);
fprintf('    material force: max |H| after the first edge %.4f N (planar %.4f N), H on the first edge %.4f N (planar %.4f N)\n', ...
    P0.H_absmax_after_first_edge, planar.base.H_absmax_after_first_edge, P0.H_first_edge, planar.base.H_first_edge);
fprintf('    line load (vector form) %.4f N vs multipliers %.4f N (%.2f%%); planar code %.4f vs %.4f N\n', ...
    P0.F_line_vec, P0.F_kkt, 100*(P0.F_line_vec/P0.F_kkt - 1), planar.base.F_line, planar.base.F_kkt);
fprintf('  KKT diagnostics, planar geometry, rigid vessel\n');
S.planar_kkt_rigid = kkt_diagnostics3d(rig0.prob, rig0.x, rig0.out);
fprintf('  KKT diagnostics, planar geometry, deformable vessel\n');
S.planar_kkt = kkt_diagnostics3d(prob0, x0, out0);
fprintf('  clock shift sanity check on the planar geometry (State 1 -> State 2):\n');
tg0 = [50 100 150 200 250 300];
C00 = clock_shift3d(G0.P, P0.Pv, [0 0 1], tg0, G0.hv);
C01 = clock_shift3d(G0.P, P0.Pv, [0 1 0], tg0, G0.hv);
fprintf('    a = (0,0,1): max |shift| %.2e deg, max |integral shift| %.2e deg\n', max(abs(C00.shift)), max(abs(C00.shift_int)));
fprintf('    a = (0,1,0): max |shift| %.2e deg, max |integral shift| %.2e deg, max |a.t| %.3f / %.3f\n', ...
    max(abs(C01.shift)), max(abs(C01.shift_int)), C01.max_at1, C01.max_at2);

%% 4. Non-planar base case.
fprintf('\n4. Non-planar base case, A = %g mm, tip 120 mm into the aorta\n', A);
t4 = tic;
[x1, out1, prob1, steps1, rig1] = continuation_solve3d(G, G.S_il + 120, B, kbar);
time_np = toc(t4);
P1 = postprocess3d(prob1, x1, out1);
R1 = postprocess3d(rig1.prob, rig1.x, rig1.out);
fprintf('  wire nodes %d (length %.2f mm), vessel nodes %d, unknowns %d\n', prob1.M, prob1.ell, prob1.N, prob1.n);
fprintf('  final solve: flag %d, first-order optimality %.1e, constraint violation %.1e, %d iterations\n', ...
    P1.flag, P1.firstorderopt, P1.constrviolation, P1.iterations);
fprintf('  rigid vessel: energy %.6f, wall reaction %.3f N, contacts:\n', R1.E, R1.F_kkt);
print_arcs(R1.arcs);
fprintf('  deformable vessel: energy %.6f N mm (wire %.4f, vessel bending %.4f, foundation %.4f)\n', ...
    P1.E, P1.E_wire, P1.E_vessel, P1.E_foundation);
fprintf('  wall reaction from the multipliers %.3f N, contacts:\n', P1.F_kkt);
print_arcs(P1.arcs);
fprintf('  peak vessel displacement %.4f mm at s_v = %.1f mm (%s, node %d), vector (%.2f, %.2f, %.2f) mm\n', ...
    P1.peak, P1.peak_at_s, ifelse(P1.peak_in_iliac, 'iliac', 'aorta'), P1.peak_node, P1.peak_vector);
fprintf('  largest displacement in the first half of the iliac %.2f mm, in the aorta %.2f mm, RMS %.2f mm\n', ...
    P1.peak_first_bend, P1.aorta_shift, P1.rms);
fprintf('  curvature at the sheath exit %.4f 1/mm, B kappa^2/2 there %.3f N\n', P1.kappa_exit, P1.halfBk2_exit);

fprintf('  (c) material force H = 2 h nu on each edge: max |H| after the first edge %.4f N, H on the first edge %.4f N,\n', ...
    P1.H_absmax_after_first_edge, P1.H_first_edge);
fprintf('      max of B kappa^2/2 %.3f N, ratio %.4f\n', P1.halfBk2_max, P1.H_absmax_after_first_edge/P1.halfBk2_max);
fprintf('      (planar geometry, same code: max |H| %.4f N, max B kappa^2/2 %.3f N)\n', P0.H_absmax_after_first_edge, P0.halfBk2_max);

fprintf('  (d) line load with H = 0 against the multipliers\n');
print_lineload(P1, '      non-planar');
print_lineload(P0, '      planar    ');
fprintf('      per arc (non-planar, vector form):\n');
for q = 1:numel(P1.arcs)
    fprintf('        s = %6.1f to %6.1f mm   multipliers %.3f N   line load %.3f N\n', P1.arcs(q).s0, P1.arcs(q).s1, P1.arcs(q).force, P1.arcs(q).force_line);
end

fprintf('  (e) KKT diagnostics, non-planar, rigid vessel\n');
S.kkt_rigid = kkt_diagnostics3d(rig1.prob, rig1.x, rig1.out);
fprintf('  (e) KKT diagnostics, non-planar, deformable vessel\n');
S.kkt = kkt_diagnostics3d(prob1, x1, out1);

%% 5. Clock shift.
fprintf('\n5. Clock shift of a twist-free wall target, a = (%g,%g,%g)\n', a_ant);
tg = [50 100 150 200 250 300];
C1 = clock_shift3d(G.P, P1.Pv, a_ant, tg, G.hv);
fprintf('   max |a.t| on the centerline: State 1 %.3f, State 2 %.3f (a is never parallel to the centerline)\n', C1.max_at1, C1.max_at2);
fprintf('   convention: 12 o''clock = projection of a, clockwise looking downstream (sheath toward aorta);\n');
fprintf('   the target is at 12 o''clock at the first vessel node and carried by the Bishop frame (no twist)\n');
fprintf('   s (mm)  segment  State 1          State 2          shift (deg)  shift (clock min)  shift from rate integrals (deg)  difference (deg)\n');
for q = 1:numel(tg)
    fprintf('   %6.0f  %-7s  %6.2f deg %5s  %6.2f deg %5s  %+10.3f  %+16.2f  %+30.3f  %+12.1e\n', tg(q), ...
        ifelse(tg(q) <= G.S_il, 'iliac', 'aorta'), C1.theta1(q), C1.clock1{q}, C1.theta2(q), C1.clock2{q}, ...
        C1.shift(q), C1.shift_min(q), C1.shift_int(q), C1.agree(q));
end
fprintf('   State 1 clock angle from the rate integral vs the transported frame: max difference %.2e deg\n', max(abs(C1.int1 - C1.theta1)));
fprintf('   (g) planar geometry: max |shift| %.2e deg (a = (0,0,1)), %.2e deg (a = (0,1,0))\n', max(abs(C00.shift)), max(abs(C01.shift)));

%% 6. Timing.
fprintf('\n6. Timing (wall clock, s)\n');
print_steps('planar geometry', steps0, time_planar);
print_steps('non-planar', steps1, time_np);
fprintf('   total run %.1f s\n', toc(t_all));

%% Save.
S.A = A;  S.geometry = rmfield(G, 'P');  S.geometry.P = G.P;
S.planar = strip(P0);  S.planar_rigid = strip(R0);
S.base = strip(P1);    S.rigid = strip(R1);
S.steps_planar = steps0;  S.steps = steps1;
S.clock = C1;  S.clock_planar_z = C00;  S.clock_planar_y = C01;
S.reference_planar = struct('F_rigid', planar.rigid.F_kkt, 'F_deformable', planar.base.F_kkt, ...
    'peak', planar.base.peak, 'E', planar.base.E, 'E_rigid', planar.rigid.E);
% Arrays for the figures (plain numeric, readable by scipy.io.loadmat).
fig = struct('P1', G.P, 'Pv', P1.Pv, 'Pw', P1.Pw, 'contact', P1.contact_nodes, 'fn', P1.fn, ...
    's_w', P1.s_w, 'cum_kkt', P1.cum_kkt, 'cum_line_vec', P1.cum_line_vec, 'cum_line_frenet', P1.cum_line_frenet, ...
    'cum_kkt_frenet_nodes', P1.cum_kkt_frenet_nodes, 'He', P1.He, 'halfBk2', P1.halfBk2, 'S_il', G.S_il, ...
    'P1_planar', G0.P, 'Pv_planar', P0.Pv, 'Pw_planar', P0.Pw, 'He_planar', P0.He, 'halfBk2_planar', P0.halfBk2, ...
    'h', P1.h, 'R', G.R);
save(fullfile(resdir, [tag '.mat']), 'S', 'fig', '-v7');
J = S;
J.planar = rmfield(J.planar, big_fields());  J.planar_rigid = rmfield(J.planar_rigid, big_fields());
J.base = rmfield(J.base, big_fields());      J.rigid = rmfield(J.rigid, big_fields());
J.geometry = rmfield(J.geometry, {'P', 'R'});
fid = fopen(fullfile(resdir, [tag '_summary.json']), 'w');
fprintf(fid, '%s', jsonencode(J, 'PrettyPrint', true));  fclose(fid);
fprintf('\nsaved results/%s.mat, results/%s_summary.json, results/%s.log\ndone\n', tag, tag, tag);
end

function f = big_fields()
f = {'Pv', 'Pw', 'P', 'R', 'U', 'Uv', 's_v', 's_w', 'fn', 'fvec', 'contact_nodes', 'sigma', 'kap', 'tau', 'kpp', ...
    'lam_frenet', 'fb_frenet', 'fline', 'fline_mag', 'valid_frenet', 'angle_force_vs_lineload', 'cum_kkt', 'cum_line_vec', ...
    'cum_line_frenet', 'cum_kkt_frenet_nodes', 'He', 'N_edge', 'halfBk2'};
end

function T = strip(P)
% Keep everything; the struct is saved to the .mat file as is.
T = P;
end

function print_arcs(arcs)
for q = 1:numel(arcs)
    fprintf('    s = %6.1f to %6.1f mm (%2d nodes)   vessel sigma %6.1f to %6.1f mm   %-5s  angle to curvature vector %5.1f deg   %.3f N\n', ...
        arcs(q).s0, arcs(q).s1, arcs(q).nodes, arcs(q).sig0, arcs(q).sig1, arcs(q).side, arcs(q).angle, arcs(q).force);
end
end

function print_arcs_compare(arcs, ref)
fprintf('    3D code:\n');
print_arcs(arcs);
fprintf('    planar code:\n');
for q = 1:numel(ref)
    fprintf('    s = %6.1f to %6.1f mm   %-5s  %.3f N\n', ref(q).s0, ref(q).s1, ref(q).side, ref(q).force);
end
end

function print_lineload(P, label)
fprintf('%s: vector form f = B y'''''''' + 3/2 B (kappa^2 y'')'' on %d of %d contact nodes: %.4f N vs multipliers %.4f N (%+.2f%%)\n', ...
    label, P.n_vec_contact, P.n_contact, P.F_line_vec, P.F_kkt, 100*(P.F_line_vec/P.F_kkt - 1));
fprintf('%s  angle between contact force and line-load vector on contact nodes: median %.1f deg, max %.1f deg\n', ...
    blanks(numel(label)), P.median_angle_force_vs_lineload, P.max_angle_force_vs_lineload);
b = blanks(numel(label));
fprintf('%s  Frenet form on the %d contact nodes with kappa > %.0e 1/mm at j-2..j+2 (contact force there %.4f N by multipliers):\n', ...
    b, P.n_frenet_valid_contact, P.kappa_min, P.F_kkt_on_frenet_nodes);
fprintf('%s    normal part   -B(kappa'''' + kappa^3/2 - kappa tau^2): %.4f N  vs normal part of the contact forces %.4f N (%+.2f%%)\n', ...
    b, P.F_line_frenet, P.F_kkt_normal, 100*(P.F_line_frenet/max(P.F_kkt_normal, eps) - 1));
fprintf('%s    binormal part  B(2 kappa'' tau + kappa tau''):        %.4f N  vs binormal part of the contact forces %.4f N\n', ...
    b, P.F_frenet_binormal, P.F_kkt_binormal);
fprintf('%s    magnitude of (normal, binormal):                   %.4f N  vs contact forces %.4f N (%+.2f%%)\n', ...
    b, P.F_frenet_full, P.F_kkt_on_frenet_nodes, 100*(P.F_frenet_full/max(P.F_kkt_on_frenet_nodes, eps) - 1));
fprintf('%s    normal part alone vs full contact force: %+.2f%%; largest binormal share |f.b|/|f| of the line load on these nodes %.3f\n', ...
    b, 100*(P.F_line_frenet/max(P.F_kkt_on_frenet_nodes, eps) - 1), P.max_binormal_share);
fprintf('%s  contact force on nodes where the Frenet form is not evaluated (kappa too small): %.4f N of %.4f N\n', ...
    b, P.F_kkt - P.F_kkt_on_frenet_nodes, P.F_kkt);
end

function print_steps(label, steps, total)
fprintf('   %s:\n', label);
for q = 1:numel(steps)
    if isinf(steps(q).k), ks = 'rigid'; else, ks = sprintf('k %.4g N/mm', steps(q).k); end
    fprintf('     %-5s  B %7.1f  %-16s flag %d  %4d it  %6.2f s\n', steps(q).stage, steps(q).B, ks, steps(q).flag, steps(q).iterations, steps(q).time);
end
fprintf('     continuation total %.1f s\n', total);
end

function s = tf(b)
if b, s = 'no self-intersection'; else, s = 'CHECK FAILED'; end
end

function out = ifelse(c, a, b)
if c, out = a; else, out = b; end
end
