function err = check_derivatives3d(A)
%CHECK_DERIVATIVES3D  Compare the analytic gradient, constraint Jacobians and
%   Lagrangian Hessian with central finite differences on a coarse mesh of
%   the non-planar vessel, at a perturbed point where some foot points lie
%   inside segments and the lumen radius varies. Also checks the wire
%   bending energy alone at a strongly bent random polygon (turning angles
%   up to about 1 rad), where the tan(phi/2) form differs from phi.
%
%   err = CHECK_DERIVATIVES3D(A) returns the relative errors. A private random
%   stream is used, so the global random state is not touched.

if nargin < 1, A = 40; end
rs = RandStream('mt19937ar', 'Seed', 1);
G = geometry3d(6.0, A);
prob = build_problem3d(G, G.S_il + 120, 4.86e3, 0.0158*G.hv);
x = prob.x_ref + 0.8*randn(rs, prob.n, 1);
[~, g] = objective3d(x, prob);
[c, ceq, Gc, Geq] = constraints3d(x, prob);
lam.ineqnonlin = rand(rs, numel(c), 1);  lam.eqnonlin = randn(rs, numel(ceq), 1);
H = lagrangian_hessian3d(x, lam, prob);
hstep = 1e-6;  n = prob.n;
gfd = zeros(n,1);  Jfd = zeros(n, numel(c));  Jeqfd = zeros(n, numel(ceq));  Hfd = zeros(n);
for i = 1:n
    xp = x;  xm = x;  xp(i) = xp(i) + hstep;  xm(i) = xm(i) - hstep;
    gfd(i) = (objective3d(xp, prob) - objective3d(xm, prob))/(2*hstep);
    [cp, eqp] = constraints3d(xp, prob);  [cm, eqm] = constraints3d(xm, prob);
    Jfd(i,:) = (cp - cm)'/(2*hstep);  Jeqfd(i,:) = (eqp - eqm)'/(2*hstep);
    Hfd(:,i) = (lagrangian_gradient3d(xp, lam, prob) - lagrangian_gradient3d(xm, lam, prob))/(2*hstep);
end
err.gradient = norm(g - gfd)/norm(g);
err.jacobian = norm(full(Gc) - Jfd, 'fro')/norm(Jfd, 'fro');
err.jacobian_eq = norm(full(Geq) - Jeqfd, 'fro')/norm(Jeqfd, 'fro');
err.hessian = norm(full(H) - Hfd, 'fro')/norm(Hfd, 'fro');

% Wire energy alone on a strongly bent polygon.
Q = cumsum([0 0 0; randn(rs, 12, 3)], 1);
cols = reshape(1:3*13, 3, [])';  m = 3*13;
q = reshape(Q', [], 1);
fw = @(v) wire_bend3d(reshape(v, 3, [])', 1.0, cols, m);
[~, gw, Hw] = wire_bend3d(Q, 1.0, cols, m);
gwfd = zeros(m,1);  Hwfd = zeros(m);
for i = 1:m
    qp = q;  qm = q;  qp(i) = qp(i) + hstep;  qm(i) = qm(i) - hstep;
    gwfd(i) = (fw(qp) - fw(qm))/(2*hstep);
    [~, gp] = wire_bend3d(reshape(qp, 3, [])', 1.0, cols, m);
    [~, gm] = wire_bend3d(reshape(qm, 3, [])', 1.0, cols, m);
    Hwfd(:,i) = (gp - gm)/(2*hstep);
end
err.wire_gradient = norm(gw - gwfd)/norm(gwfd);
err.wire_hessian = norm(full(Hw) - Hwfd, 'fro')/norm(Hwfd, 'fro');
fprintf('relative error  gradient %.2e   Jacobian %.2e   Jacobian(eq) %.2e   Hessian %.2e\n', ...
    err.gradient, err.jacobian, err.jacobian_eq, err.hessian);
fprintf('  wire energy alone, random bent polygon: gradient %.2e   Hessian %.2e\n', ...
    err.wire_gradient, err.wire_hessian);
end
