function [x, out, prob, steps, rigid] = continuation_solve3d(G, ell, B, kbar)
%CONTINUATION_SOLVE3D  Equilibrium by continuation from a rigid vessel.
%
%   Mirrors the planar CONTINUATION_SOLVE. The wire starts on the State 1
%   centerline in a rigid vessel (held fixed, B_v -> infinity) and its
%   stiffness is raised to B in three steps. The vessel is then released and
%   the foundation relaxed from 10^2 kbar*h per node to kbar*h in steps of
%   at most half a decade (first solve at 10^1.5), halving the step whenever
%   a solve fails.
%
%   (The planar code holds the vessel "almost rigid" with springs of
%   1e3*h/1.6 N/mm per node; here the rigid stage holds it exactly, which
%   is the B_v -> infinity limit of the linear beam with held ends.)
%
%   Returns the final solution x, its solver output, the problem struct, a
%   log of every step (stiffnesses, exit flag, energy, time), and rigid,
%   the solution with the full wire stiffness in the rigid vessel.

hv = G.hv;
steps = struct('stage', {}, 'B', {}, 'k', {}, 'flag', {}, 'E', {}, 'firstorderopt', {}, 'iterations', {}, 'time', {});
x = [];
for Bi = [291, 1849, B]
    prob = build_problem3d(G, ell, Bi, kbar*hv, [], true);
    if isempty(x), x = prob.x_ref; warm = false; else, warm = true; end
    [x, out] = solve_equilibrium3d(prob, x, warm);
    steps(end+1) = struct('stage', 'rigid', 'B', Bi, 'k', Inf, 'flag', out.flag, 'E', out.E, ...
        'firstorderopt', out.firstorderopt, 'iterations', out.iterations, 'time', out.time); %#ok<AGROW>
    fprintf('  rigid vessel  B = %7.1f   flag %d   E = %.6f   optimality %.1e   %d it   %.1f s\n', ...
        Bi, out.flag, out.E, out.firstorderopt, out.iterations, out.time);
end
rigid = struct('x', x, 'out', out, 'prob', prob);
[~, Pw] = unpack_nodes3d(prob, x);

prob = build_problem3d(G, ell, B, kbar*hv, [], false);
x = pack_nodes3d(prob, G.P, Pw);                  % vessel at State 1, wire from the rigid solve
j = 2.0;  step = 0.5;                              % log10(k / (kbar*h)), as in the planar code
while j > 0
    jn = max(0, j - step);
    prob = build_problem3d(G, ell, B, kbar*hv*10^jn, [], false);
    [xt, ot] = solve_equilibrium3d(prob, x, true, 400);
    steps(end+1) = struct('stage', 'relax', 'B', B, 'k', kbar*hv*10^jn, 'flag', ot.flag, 'E', ot.E, ...
        'firstorderopt', ot.firstorderopt, 'iterations', ot.iterations, 'time', ot.time); %#ok<AGROW>
    fprintf('  relax  log10(k/kbar h) = %.3f   flag %d   E = %.6f   optimality %.1e   %d it   %.1f s\n', ...
        jn, ot.flag, ot.E, ot.firstorderopt, ot.iterations, ot.time);
    if ot.flag > 0
        x = xt;  out = ot;  j = jn;  step = min(0.5, 1.5*step);
    else
        step = step/2;
        if step < 1e-3, error('continuation failed'); end
    end
end
end
