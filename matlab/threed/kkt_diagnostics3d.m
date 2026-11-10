function D = kkt_diagnostics3d(prob, x, out)
%KKT_DIAGNOSTICS3D  Constraint qualification, second-order check and discrete
%   thickness at a computed 3D equilibrium (the 3D counterpart of the planar
%   KKT_DIAGNOSTICS).
%
%   D = KKT_DIAGNOSTICS3D(prob, x, out) reports:
%     1. LICQ: the smallest singular value of the active constraint Jacobian,
%        each row scaled to unit length. The rows are the edge-length and the
%        two sheath-direction equalities and the active lumen constraints.
%     2. Second order: the smallest eigenvalue of the Hessian of the
%        Lagrangian on the null space of the equalities and the strongly
%        active lumen constraints (contact force above 1e-6 N), and on the
%        smaller null space of all active constraints. The Hessian is built
%        by central finite differences (step 1e-6 mm) of the exact gradient
%        of the Lagrangian at the solution, then symmetrized. Its relative
%        difference from the analytic Hessian of LAGRANGIAN_HESSIAN3D is
%        also reported, as a check of both.
%     3. Box bounds: the smallest distance of any unknown to its 25 mm bound.
%     4. Discrete thickness of the deformed vessel polygon: the smallest local
%        radius of curvature (Voronoi length over turning angle) and the
%        smallest half-distance between vessel nodes more than pi*delta_max
%        apart in arc length.
%     5. As in the planar code: where the contacts project on the vessel, and
%        the number of contact pairs within one node spacing of the same
%        vessel parameter whose wall normals point to opposite sides.

[Pv, Pw] = unpack_nodes3d(prob, x);
[~, ~, Gc, Geq] = constraints3d(x, prob);
[~, ~, ~, info] = contact_terms3d(Pw(2:end,:), Pv, prob.R, prob.rw, prob.wcol(2:end,:), prob.vcol, prob.n);
lam = out.lambda.ineqnonlin(:);
fn = 2*info.delta.*abs(lam);                         % nodal contact force (N)
gap = info.delta - info.d;                           % mm, >= 0 inside the lumen
act = gap <= 1e-3*info.delta;
strong = act & fn > 1e-6;

% 1. LICQ.
A = [Geq'; Gc(:, act)'];
A = A ./ sqrt(sum(A.^2, 2));
sv = svd(full(A));
D.n_unknowns = prob.n;
D.n_active_rows = size(A, 1);
D.n_contacts = nnz(act);  D.n_strong = nnz(strong);
D.licq_smin = sv(end);  D.licq_cond = sv(1)/sv(end);

% 2. Second order, with the Hessian from finite differences of the exact gradient.
n = prob.n;  hs = 1e-6;
Hfd = zeros(n);
for i = 1:n
    xp = x;  xm = x;  xp(i) = xp(i) + hs;  xm(i) = xm(i) - hs;
    Hfd(:,i) = (lagrangian_gradient3d(xp, out.lambda, prob) - lagrangian_gradient3d(xm, out.lambda, prob))/(2*hs);
end
H = (Hfd + Hfd')/2;
Han = full(lagrangian_hessian3d(x, out.lambda, prob));
D.hessian_fd_vs_analytic = norm(H - (Han + Han')/2, 'fro')/norm(H, 'fro');
D.hessian_fd_asymmetry = norm(Hfd - Hfd', 'fro')/norm(Hfd, 'fro');
Zs = null(full([Geq'; Gc(:, strong)']));
Za = null(full([Geq'; Gc(:, act)']));
D.dim_null_strong = size(Zs, 2);  D.dim_null_active = size(Za, 2);
[Vs, es] = eig((Zs'*H*Zs + (Zs'*H*Zs)')/2, 'vector');
ea = eig((Za'*H*Za + (Za'*H*Za)')/2);
[D.min_eig_strong, imin] = min(es);  D.max_eig_strong = max(es);
D.min_eig_active = min(ea);
% Where the softest admissible mode lives: share of its norm on the wire,
% and the wire node (arc length) or vessel node with the largest amplitude.
v = Zs*Vs(:, imin);
vw = reshape(v(prob.nv+1:end), 3, [])';  aw = sqrt(sum(vw.^2, 2));
D.soft_mode_wire_share = sum(aw.^2)/sum(v.^2);
[~, jw] = max(aw);  D.soft_mode_wire_peak_s = jw*prob.h0;          % wire node y_jw
t_loc = Pw(min(jw+2, prob.M),:) - Pw(jw,:);  t_loc = t_loc/norm(t_loc);
D.soft_mode_wire_peak_axial = abs(dot(vw(jw,:)/aw(jw), t_loc));  % 1 = along the wire, 0 = across it
if ~prob.rigid
    vv = reshape(v(1:prob.nv), 3, [])';  av = sqrt(sum(vv.^2, 2));
    [~, iv] = max(av);  D.soft_mode_vessel_peak_s = iv*prob.hv;   % vessel node iv+1, State 1 arc length
else
    D.soft_mode_vessel_peak_s = NaN;
end

% 3. Box bounds.
D.box_margin = min([x - (prob.x_ref - 25); (prob.x_ref + 25) - x]);

% 4. Discrete thickness of the deformed vessel polygon.
e = diff(Pv);  L = sqrt(sum(e.^2, 2));
cr = cross(e(1:end-1,:), e(2:end,:), 2);
phi = atan2(sqrt(sum(cr.^2, 2)), sum(e(1:end-1,:).*e(2:end,:), 2));
D.min_local_radius = min(0.5*(L(1:end-1) + L(2:end))./max(phi, eps));
s = [0; cumsum(L)];
dmax = max(prob.R) - prob.rw;
far = abs(s - s') > pi*dmax;
DX = sqrt((Pv(:,1) - Pv(:,1)').^2 + (Pv(:,2) - Pv(:,2)').^2 + (Pv(:,3) - Pv(:,3)').^2);
D.min_half_chord = 0.5*min(DX(far));
D.delta_max = dmax;

% 5. Where the contacts project, and opposite-wall pairs.
J = find(act);
sig = s(info.seg(J)) + info.t(J).*L(info.seg(J));     % vessel parameter of each foot point
D.contact_sigma_range = [min(sig), max(sig)];
D.vessel_length = s(end);
nrm = zeros(numel(J), 3);
for q = 1:numel(J)
    f = Pv(info.seg(J(q)),:) + info.t(J(q))*e(info.seg(J(q)),:);
    y = Pw(J(q)+1,:);
    nrm(q,:) = (y - f)/norm(y - f);
end
opp = 0;
for p = 1:numel(J)
    for q = p+1:numel(J)
        if abs(sig(p) - sig(q)) <= prob.h0 && dot(nrm(p,:), nrm(q,:)) < 0, opp = opp + 1; end
    end
end
D.opposite_wall_pairs = opp;

fprintf('    LICQ: %d active rows, %d unknowns, smallest singular value %.3e (condition %.1e)\n', ...
    D.n_active_rows, D.n_unknowns, D.licq_smin, D.licq_cond);
fprintf('    Hessian of the Lagrangian by finite differences of the exact gradient: asymmetry %.1e, difference from the analytic Hessian %.1e\n', ...
    D.hessian_fd_asymmetry, D.hessian_fd_vs_analytic);
fprintf('    second order: min eigenvalue of the reduced Hessian %.3e N/mm (strongly active, dim %d), %.3e N/mm (all active, dim %d)\n', ...
    D.min_eig_strong, D.dim_null_strong, D.min_eig_active, D.dim_null_active);
fprintf('    largest eigenvalue of the reduced Hessian %.2e N/mm\n', D.max_eig_strong);
fprintf('    softest mode: %.0f%% of its norm on the wire, largest wire amplitude at s = %.1f mm (direction %.2f along the wire), largest vessel amplitude at s_v = %.1f mm\n', ...
    100*D.soft_mode_wire_share, D.soft_mode_wire_peak_s, D.soft_mode_wire_peak_axial, D.soft_mode_vessel_peak_s);
fprintf('    box bounds: smallest margin %.2f mm\n', D.box_margin);
fprintf('    vessel polygon: min local radius %.1f mm, min half-chord %.1f mm, delta_max %.2f mm\n', ...
    D.min_local_radius, D.min_half_chord, D.delta_max);
fprintf('    contacts project to sigma in [%.1f, %.1f] of [0, %.1f] mm, opposite-wall pairs %d\n', ...
    D.contact_sigma_range(1), D.contact_sigma_range(2), D.vessel_length, D.opposite_wall_pairs);
end
