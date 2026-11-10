function gl = lagrangian_gradient3d(x, lam, prob)
%LAGRANGIAN_GRADIENT3D  Exact gradient of the Lagrangian,
%   grad E + sum_j lambda_ineq_j grad c_j + sum_j lambda_eq_j grad ceq_j,
%   with lam in the fmincon layout (fields ineqnonlin, eqnonlin). Used for
%   the derivative check and for the finite-difference Hessian in
%   KKT_DIAGNOSTICS3D.
[~, g] = objective3d(x, prob);
[~, ~, Gc, Geq] = constraints3d(x, prob);
gl = g + Gc*lam.ineqnonlin + Geq*lam.eqnonlin;
end
