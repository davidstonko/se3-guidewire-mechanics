function [x, out] = solve_equilibrium3d(prob, x0, warm, maxit)
%SOLVE_EQUILIBRIUM3D  One interior-point solve of the discrete 3D problem.
%
%   [x, out] = SOLVE_EQUILIBRIUM3D(prob, x0, warm, maxit) minimizes the energy
%   of OBJECTIVE3D subject to CONSTRAINTS3D with fmincon (interior-point),
%   using the exact gradient, constraint Jacobian and Hessian of the
%   Lagrangian. Settings are those of the planar SOLVE_EQUILIBRIUM: node
%   positions bounded to within 25 mm of the reference point, tolerances 1e-9.
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
    'HessianFcn', @(x, lam) lagrangian_hessian3d(x, lam, prob), ...
    'SubproblemAlgorithm', 'factorization', ...
    'InitBarrierParam', mu0, ...
    'OptimalityTolerance', 1e-9, 'ConstraintTolerance', 1e-9, 'StepTolerance', 1e-15, ...
    'MaxIterations', maxit, 'MaxFunctionEvaluations', 1e7, 'Display', 'off');
t0 = tic;                               % private timer: the shared session timer is not touched
[x, fval, flag, op, lam] = fmincon(@(x) objective3d(x, prob), x0, [], [], [], [], lb, ub, ...
    @(x) constraints3d(x, prob), opts);
out = struct('flag', flag, 'E', fval, 'firstorderopt', op.firstorderopt, ...
    'constrviolation', op.constrviolation, 'iterations', op.iterations, ...
    'lambda', lam, 'time', toc(t0));
end
