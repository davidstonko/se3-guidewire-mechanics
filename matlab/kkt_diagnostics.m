function D = kkt_diagnostics(prob, x, out)
%KKT_DIAGNOSTICS  Constraint qualification, second-order check and discrete
%   thickness at a computed equilibrium.
%
%   D = KKT_DIAGNOSTICS(prob, x, out) uses the solution x and the solver
%   output out (with its Lagrange multipliers) and reports:
%     1. LICQ: the smallest singular value of the active constraint Jacobian,
%        each row scaled to unit length. The rows are the edge-length and
%        sheath-direction equalities and the active lumen constraints. A
%        positive value means the multipliers, and so the contact forces, are
%        unique.
%     2. Second order: the smallest eigenvalue of the Hessian of the
%        Lagrangian on the null space of the equalities and the strongly
%        active lumen constraints (multiplier above 1e-6 N), and on the
%        smaller null space of all active constraints. A positive value on the
%        first space is the strong second-order sufficient condition, so the
%        point is a strict local minimizer of the discrete problem.
%     3. Box bounds: the smallest distance of any unknown to its 25 mm bound.
%     4. Discrete thickness of the deformed vessel polygon: the smallest local
%        radius of curvature (Voronoi length over turning angle) and the
%        smallest half-distance between vessel nodes more than pi*delta_max
%        apart in arc length. If both exceed delta_max, the tube of radius
%        delta_max about the polygon does not overlap itself.
%     5. The conditions of Proposition 3.2.3: every contact projects inside
%        the vessel, and no two contacts within one node spacing of the same
%        vessel parameter press on opposite walls.

[Pv, Pw] = unpack_nodes(prob, x);
[c, ~, Gc, Geq] = constraints(x, prob);
[~, ~, ~, info] = contact_terms(Pw(2:end,:), Pv, prob.R, prob.rw, prob.wcol(2:end,:), prob.vcol, prob.n);
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

% 2. Second order.
H = lagrangian_hessian(x, out.lambda, prob);
H = full((H + H')/2);
Zs = null(full([Geq'; Gc(:, strong)']));
Za = null(full([Geq'; Gc(:, act)']));
D.dim_null_strong = size(Zs, 2);  D.dim_null_active = size(Za, 2);
D.min_eig_strong = min(eig(Zs'*H*Zs));
D.min_eig_active = min(eig(Za'*H*Za));
D.max_eig_strong = max(eig(Zs'*H*Zs));

% 3. Box bounds.
D.box_margin = min([x - (prob.x_ref - 25); (prob.x_ref + 25) - x]);

% 4. Discrete thickness of the deformed vessel polygon.
e = diff(Pv);  L = sqrt(sum(e.^2, 2));
phi = atan2(e(1:end-1,1).*e(2:end,2) - e(1:end-1,2).*e(2:end,1), sum(e(1:end-1,:).*e(2:end,:), 2));
D.min_local_radius = min(0.5*(L(1:end-1) + L(2:end))./max(abs(phi), eps));
s = [0; cumsum(L)];
dmax = max(prob.R) - prob.rw;
far = abs(s - s') > pi*dmax;
DX = sqrt((Pv(:,1) - Pv(:,1)').^2 + (Pv(:,2) - Pv(:,2)').^2);
D.min_half_chord = 0.5*min(DX(far));
D.delta_max = dmax;

% 5. Conditions of Proposition 3.2.3.
J = find(act);
sig = s(info.seg(J)) + info.t(J).*L(info.seg(J));     % vessel parameter of each foot point
D.contact_sigma_range = [min(sig), max(sig)];
D.vessel_length = s(end);
nrm = zeros(numel(J), 2);
for q = 1:numel(J)
    a = Pv(info.seg(J(q)),:);  ee = e(info.seg(J(q)),:);
    f = a + info.t(J(q))*ee;  y = Pw(J(q)+1,:);
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
fprintf('    second order: min eigenvalue of the reduced Hessian %.3e N/mm (strongly active, dim %d), %.3e N/mm (all active, dim %d)\n', ...
    D.min_eig_strong, D.dim_null_strong, D.min_eig_active, D.dim_null_active);
fprintf('    box bounds: smallest margin %.2f mm\n', D.box_margin);
fprintf('    vessel polygon: min local radius %.1f mm, min half-chord %.1f mm, delta_max %.2f mm\n', ...
    D.min_local_radius, D.min_half_chord, D.delta_max);
fprintf('    Proposition 3.2.3: contacts project to sigma in [%.1f, %.1f] of [0, %.1f] mm, opposite-wall pairs %d\n', ...
    D.contact_sigma_range(1), D.contact_sigma_range(2), D.vessel_length, D.opposite_wall_pairs);
end
