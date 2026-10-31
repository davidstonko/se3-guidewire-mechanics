function [x, out] = solve_equilibrium(prob, x0, warm, maxit)
%SOLVE_EQUILIBRIUM  One interior-point solve of the discrete problem.
%
%   [x, out] = SOLVE_EQUILIBRIUM(prob, x0, warm, maxit) minimizes the energy
%   of OBJECTIVE subject to CONSTRAINTS with fmincon (interior-point), using
%   the exact gradient, constraint Jacobian and Hessian of the Lagrangian.
%   Node positions are bounded to within 25 mm of the reference point, a
%   bound that is inactive at every reported solution.
%
%   warm = true restarts from a nearby equilibrium with a small initial
%   barrier parameter. out holds the exit flag, the first-order optimality
%   measure, the Lagrange multipliers and the wall-clock time.

if nargin < 3 || isempty(warm), warm = false; end
if nargin < 4 || isempty(maxit), maxit = 3000; end
box = 25;
lb = prob.x_ref - box;  ub = prob.x_ref + box;
mu0 = 0.1;  if warm, mu0 = 1e-4; end
opts = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
    'HessianFcn', @(x, lam) lagrangian_hessian(x, lam, prob), ...
    'SubproblemAlgorithm', 'factorization', ...
    'InitBarrierParam', mu0, ...
    'OptimalityTolerance', 1e-9, 'ConstraintTolerance', 1e-9, 'StepTolerance', 1e-15, ...
    'MaxIterations', maxit, 'MaxFunctionEvaluations', 1e7, 'Display', 'off');
tic;
[x, fval, flag, op, lam] = fmincon(@(x) objective(x, prob), x0, [], [], [], [], lb, ub, ...
    @(x) constraints(x, prob), opts);
out = struct('flag', flag, 'E', fval, 'firstorderopt', op.firstorderopt, ...
    'constrviolation', op.constrviolation, 'iterations', op.iterations, ...
    'lambda', lam, 'time', toc);
end
