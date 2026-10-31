function check_derivatives()
%CHECK_DERIVATIVES  Compare the analytic gradient, constraint Jacobian and
%   Lagrangian Hessian with central finite differences on a coarse mesh,
%   at a perturbed point where some foot points lie inside segments and the
%   lumen radius varies.

G = evar_geometry(6.0);
prob = build_problem(G, G.S_il + 120, 4.86e3, 0.0158*G.hv);
rng(1);
x = prob.x_ref + 0.8*randn(prob.n, 1);
[~, g] = objective(x, prob);
[c, ceq, Gc, Geq] = constraints(x, prob);
lam.ineqnonlin = rand(numel(c), 1);  lam.eqnonlin = randn(numel(ceq), 1);
H = lagrangian_hessian(x, lam, prob);
hstep = 1e-6;  n = prob.n;
gfd = zeros(n,1);  Jfd = zeros(n, numel(c));  Jeqfd = zeros(n, numel(ceq));  Hfd = zeros(n);
lagr_grad = @(z) lagrangian_gradient(z, lam, prob);
for i = 1:n
    xp = x;  xm = x;  xp(i) = xp(i) + hstep;  xm(i) = xm(i) - hstep;
    gfd(i) = (objective(xp, prob) - objective(xm, prob))/(2*hstep);
    [cp, eqp] = constraints(xp, prob);  [cm, eqm] = constraints(xm, prob);
    Jfd(i,:) = (cp - cm)'/(2*hstep);  Jeqfd(i,:) = (eqp - eqm)'/(2*hstep);
    Hfd(:,i) = (lagr_grad(xp) - lagr_grad(xm))/(2*hstep);
end
fprintf('relative error  gradient %.2e   Jacobian %.2e   Jacobian(eq) %.2e   Hessian %.2e\n', ...
    norm(g - gfd)/norm(g), norm(full(Gc) - Jfd, 'fro')/norm(Jfd, 'fro'), ...
    norm(full(Geq) - Jeqfd, 'fro')/norm(Jeqfd, 'fro'), norm(full(H) - Hfd, 'fro')/norm(Hfd, 'fro'));
end

function gl = lagrangian_gradient(x, lam, prob)
[~, g] = objective(x, prob);
[~, ~, Gc, Geq] = constraints(x, prob);
gl = g + Gc*lam.ineqnonlin + Geq*lam.eqnonlin;
end
