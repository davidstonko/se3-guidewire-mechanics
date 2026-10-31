function H = lagrangian_hessian(x, lambda, prob)
%LAGRANGIAN_HESSIAN  Exact Hessian of the Lagrangian for fmincon,
%
%   L(x) = E(x) + sum_j lambda_ineq_j c_j(x) + sum_j lambda_eq_j ceq_j(x).
%
%   The edge-length constraint |y_j - y_{j-1}|^2 - h^2 has Hessian
%   [2I -2I; -2I 2I] in (y_j, y_{j-1}); the sheath constraint is linear.

[~, ~, H] = objective(x, prob);
[Pv, Pw] = unpack_nodes(prob, x);
n = prob.n;  M = prob.M;
[~, ~, Hc] = contact_terms(Pw(2:end,:), Pv, prob.R, prob.rw, prob.wcol(2:end,:), prob.vcol, n, lambda.ineqnonlin);
nu = lambda.eqnonlin(1:M-1);
cj = prob.wcol(2:end,:);  cm = prob.wcol(1:end-1,:);
I = []; J = []; V = [];
for d = 1:2
    a = cj(:,d);  b = cm(:,d);
    I = [I; a; b; a; b];  J = [J; a; b; b; a];          %#ok<AGROW>
    V = [V; 2*nu; 2*nu; -2*nu; -2*nu];                  %#ok<AGROW>
end
keep = I > 0 & J > 0;
He = sparse(I(keep), J(keep), V(keep), n, n);
H = H + Hc + He;
end
