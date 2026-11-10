function [x, out, prob, steps, rigid] = continuation_solve_strict(G, ell, B, kbar)
%CONTINUATION_SOLVE_STRICT  CONTINUATION_SOLVE that accepts a step only with
%   the first-order certificate (fmincon exit flag 1).
%
%   CONTINUATION_SOLVE accepts any positive exit flag. At h = 0.8 mm that
%   lets through solves that stop on the step tolerance (flag 2) far from an
%   equilibrium. Here a solve that ends with another flag is restarted from
%   the point it reached, up to three times, and the continuation step is
%   halved only if that also fails. Same inputs and outputs.

hv = G.hv;
k_rigid = 1e3*hv/1.6;
steps = struct('B', {}, 'k', {}, 'flag', {}, 'E', {}, 'opt', {}, 'time', {});
x = [];
for Bi = [291, 1849, B]
    prob = build_problem(G, ell, Bi, k_rigid);
    if isempty(x), x = prob.x_ref; warm = false; else, warm = true; end
    [x, out] = polish(prob, x, warm, 3000);
    fprintf('  rigid vessel  B = %7.1f   flag %d   E = %.6f   opt %.1e   %.1f s\n', Bi, out.flag, out.E, out.firstorderopt, out.time);
    if out.flag ~= 1, error('rigid stage did not reach the first-order certificate'); end
end
rigid = struct('x', x, 'out', out, 'prob', prob);
j = 2.0;  step = 0.5;                              % log10(k / (kbar*h))
while j > 0
    jn = max(0, j - step);
    prob = build_problem(G, ell, B, kbar*hv*10^jn);
    [xt, ot] = polish(prob, x, true, 400);
    fprintf('  relax  log10(k/kbar h) = %.3f   flag %d   E = %.6f   opt %.1e   %.1f s\n', jn, ot.flag, ot.E, ot.firstorderopt, ot.time);
    if ot.flag == 1
        x = xt;  out = ot;  j = jn;  step = min(0.5, 1.5*step);
    else
        step = step/2;
        if step < 1e-2, error('continuation failed'); end
    end
end

    function [xs, os] = polish(prob, x0, warm, maxit)
        tt = 0;
        [xs, os] = solve_equilibrium(prob, x0, warm, maxit);  tt = tt + os.time;
        for r = 1:3
            if os.flag == 1, break, end
            [xs, os] = solve_equilibrium(prob, xs, true, maxit);  tt = tt + os.time;
        end
        os.time = tt;
        steps(end+1) = struct('B', prob.B, 'k', prob.k(1), 'flag', os.flag, 'E', os.E, 'opt', os.firstorderopt, 'time', tt);
    end
end
