function [x, out, prob, steps, rigid] = continuation_solve(G, ell, B, kbar)
%CONTINUATION_SOLVE  Equilibrium by continuation from a nearly rigid vessel.
%
%   The problem is not convex, so the starting point matters. The wire starts
%   on the State 1 centerline and the vessel is held almost rigid by per-node
%   springs of 1e3*h/1.6 N/mm. The wire stiffness is raised to B in three
%   steps, then the foundation is relaxed to kbar*h per node in steps of at
%   most half a decade, halving the step whenever a solve fails.
%
%   Returns the final solution x, its solver output, the problem struct and a
%   log of every step (stiffnesses, exit flag, energy, time). rigid holds the
%   solution with the full wire stiffness in the nearly rigid vessel.

hv = G.hv;
k_rigid = 1e3*hv/1.6;
steps = struct('B', {}, 'k', {}, 'flag', {}, 'E', {}, 'time', {});
x = [];
for Bi = [291, 1849, B]
    prob = build_problem(G, ell, Bi, k_rigid);
    if isempty(x), x = prob.x_ref; warm = false; else, warm = true; end
    [x, out] = solve_equilibrium(prob, x, warm);
    steps(end+1) = struct('B', Bi, 'k', k_rigid, 'flag', out.flag, 'E', out.E, 'time', out.time); %#ok<AGROW>
    fprintf('  rigid vessel  B = %7.1f   flag %d   E = %.6f   %.1f s\n', Bi, out.flag, out.E, out.time);
end
rigid = struct('x', x, 'out', out, 'prob', prob);
j = 2.0;  step = 0.5;                              % log10(k / (kbar*h))
while j > 0
    jn = max(0, j - step);
    prob = build_problem(G, ell, B, kbar*hv*10^jn);
    [xt, ot] = solve_equilibrium(prob, x, true, 400);
    steps(end+1) = struct('B', B, 'k', kbar*hv*10^jn, 'flag', ot.flag, 'E', ot.E, 'time', ot.time); %#ok<AGROW>
    fprintf('  relax  log10(k/kbar h) = %.3f   flag %d   E = %.6f   %.1f s\n', jn, ot.flag, ot.E, ot.time);
    if ot.flag > 0
        x = xt;  out = ot;  j = jn;  step = min(0.5, 1.5*step);
    else
        step = step/2;
        if step < 1e-3, error('continuation failed'); end
    end
end
end
