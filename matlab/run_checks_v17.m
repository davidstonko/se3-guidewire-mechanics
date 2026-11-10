function R = run_checks_v17(part, hs, Rc)
%RUN_CHECKS_V17  Further checks of the planar computation of Section 7.
%
%   R = RUN_CHECKS_V17(part), with part one of:
%     'kkt'    base case at h = 1.6 mm, rigid and deformable vessel:
%              strict complementarity, the nearest-segment choice in
%              CONTACT_TERMS, and the second-order condition measured in a
%              discrete H^1 metric as well as in the Euclidean one
%     'const'  the same vessel with a constant lumen radius of 6 mm, so that
%              no contact can sit where the lumen widens (Assumption 5.1(c))
%     'mesh'   the base case at h = 1.6, 0.8 and 0.4 mm, or at the spacings
%              hs given as a second argument (the log is then
%              results/checks_mesh_h<hs>.log)
%     'nested' the base case at h = 1.6, 0.8 and 0.4 mm (or hs), with each
%              finer mesh started from the interpolated solution of the
%              coarser one at the target stiffnesses, not by continuation
%              in the foundation. Only solves that end with the first-order
%              certificate (exit flag 1) are accepted. A third argument Rc
%              sets a constant lumen radius (mm).
%
%   Output goes to the command window, to results/checks_<part>.log and to
%   results/checks_<part>.json. No other file is changed.
%
%   Second order in H^1. KKT_DIAGNOSTICS reports the smallest eigenvalue of
%   Z'*H*Z with Z an orthonormal basis (MATLAB null) of the null space of the
%   active constraint gradients. That number depends on the coordinates and
%   scales with h. Here we also report the smallest eigenvalue of the pencil
%   (Z'*H*Z, Z'*G*Z), where G is the Gram matrix of the discrete H^1 norm
%   h*sum|y_j|^2 + sum|y_j - y_{j-1}|^2/h of the node positions (length
%   unit 1 mm), summed over wire and vessel. This is the discrete form of
%   the second variation divided by the squared H^1 norm, and should
%   converge as h -> 0.

here = fileparts(mfilename('fullpath'));
old = cd(here);  restore = onCleanup(@() cd(old));
B = 158.4e3*pi*0.889^4/64;  kbar = 0.0158;
tag = part;
if nargin > 1 && ~isempty(hs), tag = [part '_h' strrep(strjoin(compose('%g', hs), '_'), '.', 'p')]; end
if nargin > 2 && ~isempty(Rc), tag = sprintf('%s_R%g', tag, Rc); end
logf = fopen(fullfile(here, 'results', ['checks_' tag '.log']), 'w');
closer = onCleanup(@() fclose(logf));
say('%s   MATLAB %s   part %s\n', char(datetime('now')), version, part);

switch part
    case 'kkt'
        G = evar_geometry(1.6);
        [x, out, prob, ~, rigid] = continuation_solve(G, G.S_il + 120, B, kbar);
        say('Rigid vessel\n');      R.rigid = extra(rigid.prob, rigid.x, rigid.out);
        say('Deformable vessel\n'); R.base  = extra(prob, x, out);

    case 'const'
        G = evar_geometry(1.6);  G.R(:) = 6;
        [x, out, prob, ~, rigid] = continuation_solve(G, G.S_il + 120, B, kbar);
        for k = 1:2
            if k == 1, p = rigid.prob; xx = rigid.x; oo = rigid.out; nm = 'rigid';
            else,      p = prob;       xx = x;       oo = out;       nm = 'deformable'; end
            S = postprocess(p, xx, oo);
            say('Constant radius 6 mm, %s vessel: flag %d, E = %.4f, wall reaction %.3f N (line load with H = 0: %.3f N)\n', ...
                nm, S.flag, S.E, S.F_kkt, S.F_line);
            for a = 1:numel(S.arcs)
                say('    s = %6.1f to %6.1f mm   %-6s %.3f N\n', S.arcs(a).s0, S.arcs(a).s1, S.arcs(a).side, S.arcs(a).force);
            end
            say('    peak displacement %.4f mm at s_v = %.1f mm, |H| after the first edge < %.4f N, tip force %.4f N\n', ...
                S.peak, S.peak_at_s, S.H_absmax_after_first_edge, S.tip_force);
            say('    last contact s = %.1f mm, straight continuation stays in the lumen up to s = %.1f mm\n', S.s_last_contact, S.s_free_until);
            D = extra(p, xx, oo);
            R.(nm) = struct('flag', S.flag, 'E', S.E, 'F_kkt', S.F_kkt, 'F_line', S.F_line, 'arcs', S.arcs, ...
                'peak', S.peak, 'peak_at_s', S.peak_at_s, 'H_absmax', S.H_absmax_after_first_edge, ...
                'tip_force', S.tip_force, 's_last_contact', S.s_last_contact, 'kkt', D);
        end

    case 'mesh'
        if nargin < 2 || isempty(hs), hs = [1.6 0.8 0.4]; end
        for i = 1:numel(hs)
            G = evar_geometry(hs(i));
            say('h = %.2f mm (actual %.4f), %d vessel nodes\n', hs(i), G.hv, size(G.P,1));
            [x, out, prob, ~, rigid] = continuation_solve(G, G.S_il + 120, B, kbar);
            S  = postprocess(prob, x, out);
            Sr = postprocess(rigid.prob, rigid.x, rigid.out);
            say('  rigid vessel: wall reaction %.3f N, contacts:\n', Sr.F_kkt);
            for a = 1:numel(Sr.arcs)
                say('    s = %6.1f to %6.1f mm   %-6s %.3f N\n', Sr.arcs(a).s0, Sr.arcs(a).s1, Sr.arcs(a).side, Sr.arcs(a).force);
            end
            say('  deformable vessel: flag %d, E = %.4f, wall reaction %.3f N (line load with H = 0: %.3f N), contacts:\n', ...
                S.flag, S.E, S.F_kkt, S.F_line);
            for a = 1:numel(S.arcs)
                say('    s = %6.1f to %6.1f mm   %-6s %.3f N\n', S.arcs(a).s0, S.arcs(a).s1, S.arcs(a).side, S.arcs(a).force);
            end
            say('  peak %.4f mm at s_v = %.1f, first bend %.3f, aorta %.3f, RMS %.3f mm\n', ...
                S.peak, S.peak_at_s, S.peak_first_bend, S.aorta_shift, S.rms);
            say('  |H| after the first edge < %.4f N (first edge %.4f N), kappa at exit %.4f 1/mm, max B kappa^2/2 %.3f N\n', ...
                S.H_absmax_after_first_edge, S.H_first_edge, S.kappa_exit, S.halfBk2_max);
            say('  last contact s = %.1f mm, straight continuation to s = %.1f mm\n', S.s_last_contact, S.s_free_until);
            D = extra(prob, x, out);
            R.mesh(i) = struct('h', G.hv, 'n', prob.n, 'flag', S.flag, 'E', S.E, ...
                'F_kkt', S.F_kkt, 'F_line', S.F_line, 'F_rigid', Sr.F_kkt, 'arcs', S.arcs, 'arcs_rigid', Sr.arcs, ...
                'peak', S.peak, 'peak_at_s', S.peak_at_s, 'peak_first_bend', S.peak_first_bend, ...
                'aorta_shift', S.aorta_shift, 'rms', S.rms, 'H_absmax', S.H_absmax_after_first_edge, ...
                'kappa_exit', S.kappa_exit, 'halfBk2_max', S.halfBk2_max, ...
                's_last_contact', S.s_last_contact, 's_free_until', S.s_free_until, ...
                'gen_min_eig', D.gen_min_eig, 'min_eig_strong', D.min_eig_strong, ...
                'licq_smin', D.licq_smin, 'sc_min_force_active', D.sc_min_force_active);
        end
    case 'nested'
        if nargin < 2 || isempty(hs), hs = [1.6 0.8 0.4]; end
        if nargin < 3, Rc = []; end
        if ~isempty(Rc), say('Constant lumen radius %g mm\n', Rc); end
        for i = 1:numel(hs)
            G = evar_geometry(hs(i));
            if ~isempty(Rc), G.R(:) = Rc; end
            say('h = %.2f mm (actual %.4f), %d vessel nodes\n', hs(i), G.hv, size(G.P,1));
            tl = tic;
            if i == 1
                [x, out, prob, ~, rigid] = continuation_solve(G, G.S_il + 120, B, kbar);
            else
                rigid = nested(G, rigid.prob, rigid.x, 1e3*G.hv/1.6);
                say('  rigid stage from the coarser mesh: flag %d, first-order optimality %.1e, %d solves\n', rigid.out.flag, rigid.out.firstorderopt, rigid.nsolve);
                N = nested(G, prob, x, kbar*G.hv);
                x = N.x;  out = N.out;  prob = N.prob;
                say('  deformable from the coarser mesh: flag %d, first-order optimality %.1e, %d solves\n', out.flag, out.firstorderopt, N.nsolve);
            end
            say('  wall time %.0f s\n', toc(tl));
            if out.flag ~= 1 || rigid.out.flag ~= 1
                say('  NOT CONVERGED to the first-order certificate at this mesh. Stopping.\n');
                break
            end
            S  = postprocess(prob, x, out);
            Sr = postprocess(rigid.prob, rigid.x, rigid.out);
            say('  rigid vessel: wall reaction %.3f N, contacts:\n', Sr.F_kkt);
            for a = 1:numel(Sr.arcs)
                say('    s = %6.1f to %6.1f mm   %-6s %.3f N\n', Sr.arcs(a).s0, Sr.arcs(a).s1, Sr.arcs(a).side, Sr.arcs(a).force);
            end
            say('  deformable vessel: flag %d, E = %.4f, wall reaction %.3f N (line load with H = 0: %.3f N), contacts:\n', ...
                S.flag, S.E, S.F_kkt, S.F_line);
            for a = 1:numel(S.arcs)
                say('    s = %6.1f to %6.1f mm   %-6s %.3f N\n', S.arcs(a).s0, S.arcs(a).s1, S.arcs(a).side, S.arcs(a).force);
            end
            say('  peak %.4f mm at s_v = %.1f, first bend %.3f, aorta %.3f, RMS %.3f mm\n', ...
                S.peak, S.peak_at_s, S.peak_first_bend, S.aorta_shift, S.rms);
            say('  |H| after the first edge < %.4f N (first edge %.4f N), kappa at exit %.4f 1/mm, max B kappa^2/2 %.3f N\n', ...
                S.H_absmax_after_first_edge, S.H_first_edge, S.kappa_exit, S.halfBk2_max);
            say('  last contact s = %.1f mm, straight continuation to s = %.1f mm\n', S.s_last_contact, S.s_free_until);
            D = extra(prob, x, out);
            R.mesh(i) = struct('h', G.hv, 'n', prob.n, 'flag', S.flag, 'E', S.E, ...
                'F_kkt', S.F_kkt, 'F_line', S.F_line, 'F_rigid', Sr.F_kkt, 'arcs', S.arcs, 'arcs_rigid', Sr.arcs, ...
                'peak', S.peak, 'peak_at_s', S.peak_at_s, 'peak_first_bend', S.peak_first_bend, ...
                'aorta_shift', S.aorta_shift, 'rms', S.rms, 'H_absmax', S.H_absmax_after_first_edge, ...
                'kappa_exit', S.kappa_exit, 'halfBk2_max', S.halfBk2_max, ...
                's_last_contact', S.s_last_contact, 's_free_until', S.s_free_until, ...
                'gen_min_eig', D.gen_min_eig, 'min_eig_strong', D.min_eig_strong, ...
                'licq_smin', D.licq_smin, 'sc_min_force_active', D.sc_min_force_active);
        end
        if ~exist('R', 'var'), R = struct(); end
    otherwise
        error('part must be ''kkt'', ''const'', ''mesh'' or ''nested''');
end
fid = fopen(fullfile(here, 'results', ['checks_' tag '.json']), 'w');
fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true));  fclose(fid);
say('done\n');

    function say(varargin)
        fprintf(varargin{:});  fprintf(logf, varargin{:});
    end

    function N = nested(G, pc, xc, knode)
        % Start the mesh G from the solution xc of the coarser problem pc,
        % with the wire stiffness B and the nodal foundation knode.
        [Pvc, Pwc] = unpack_nodes(pc, xc);
        L  = (pc.N - 1)*pc.h0;
        sc = linspace(0, L, pc.N)';  sf = linspace(0, L, size(G.P,1))';
        Pvf = G.P + interp1(sc, Pvc - pc.P, sf, 'pchip');
        p  = build_problem(G, G.S_il + 120, B, knode);
        swc = (0:pc.M-1)'*pc.h0;  q = (0:p.M-1)'*p.h0;
        te = Pwc(end,:) - Pwc(end-1,:);  te = te/norm(te);
        Pwf = interp1(swc, Pwc, min(q, swc(end)), 'spline') + max(q - swc(end), 0).*te;
        Pwf(1,:) = p.a;
        x0 = pack_nodes(p, Pvf, Pwf);
        [xs, os] = solve_equilibrium(p, x0, true, 3000);  ns = 1;
        while os.flag ~= 1 && ns < 4
            [xs, os] = solve_equilibrium(p, xs, true, 3000);  ns = ns + 1;
        end
        if os.flag ~= 1
            [xs, os] = solve_equilibrium(p, xs, false, 3000);  ns = ns + 1;
        end
        N = struct('x', xs, 'out', os, 'prob', p, 'nsolve', ns);
    end

    function D = extra(prob, x, out)
        D = kkt_diagnostics(prob, x, out);
        [Pv, Pw] = unpack_nodes(prob, x);
        [~, ~, Gc, Geq] = constraints(x, prob);
        [~, ~, ~, info] = contact_terms(Pw(2:end,:), Pv, prob.R, prob.rw, prob.wcol(2:end,:), prob.vcol, prob.n);
        fn  = 2*info.delta.*abs(out.lambda.ineqnonlin(:));
        gap = info.delta - info.d;
        act = gap <= 1e-3*info.delta;

        % Strict complementarity.
        D.sc_min_force_active  = min(fn(act));
        D.sc_weak_active       = nnz(act & fn <= 1e-6);
        D.sc_max_gap_active    = max(gap(act));
        D.sc_min_gap_inactive  = min(gap(~act));
        D.sc_max_force_inactive = max([fn(~act); 0]);
        D.sc_n_below = [nnz(act & fn < 1e-4), nnz(act & fn < 1e-3), nnz(act & fn < 1e-2)];
        say('    strict complementarity: %d active nodes, smallest active force %.3e N, %d active with force <= 1e-6 N\n', ...
            nnz(act), D.sc_min_force_active, D.sc_weak_active);
        say('      active nodes with force below 1e-4, 1e-3, 1e-2 N: %d, %d, %d\n', D.sc_n_below);
        say('      largest active gap %.2e mm, smallest inactive gap %.2e mm, largest inactive force %.1e N\n', ...
            D.sc_max_gap_active, D.sc_min_gap_inactive, D.sc_max_force_inactive);

        % Nearest-segment choice: recompute by brute force over all segments.
        dmax = max(prob.R) - prob.rw;
        Y = Pw(2:end,:);  A0 = Pv(1:end-1,:);  E = diff(Pv);  q = sum(E.^2, 2);
        m = size(Y,1);  d1 = zeros(m,1);  d2 = d1;  tv = d1;  seg = d1;
        for j = 1:m
            w = Y(j,:) - A0;  t = min(max(sum(w.*E, 2)./q, 0), 1);
            d = sqrt(sum((w - t.*E).^2, 2));
            [d1(j), ii] = min(d);
            d2(j) = min(d(abs((1:numel(d))' - ii)*prob.h0 > pi*dmax));
            tv(j) = t(ii);  seg(j) = ii;
        end
        D.seg_agree = all(seg(:) == info.seg(:));
        D.seg_margin_active = min(d2(act) - d1(act));
        D.seg_margin_all = min(d2 - d1);
        D.foot_at_vertex_active = nnz(act & (tv < 1e-9 | tv > 1 - 1e-9));
        D.foot_min_dist_to_vertex = min(min(tv(act), 1 - tv(act)).*sqrt(q(seg(act))));
        say('    nearest segment: brute force agrees %d, margin to the nearest segment more than pi*delta_max away in arc length %.2f mm (contacts), %.2f mm (all nodes)\n', ...
            D.seg_agree, D.seg_margin_active, D.seg_margin_all);
        say('      contact foot points at a vertex %d, smallest distance of a contact foot point to a vertex %.3f mm\n', ...
            D.foot_at_vertex_active, D.foot_min_dist_to_vertex);

        % Second order in the discrete H^1 metric.
        H = lagrangian_hessian(x, out.lambda, prob);  H = full((H + H')/2);
        % The critical cone is largest when the fewest contacts count as
        % strongly active, so the condition is checked for several force
        % thresholds; a larger threshold gives a larger subspace.
        Gm = h1gram(prob);
        thr = [1e-6 1e-4 1e-3 1e-2];
        D.so_threshold = thr;  D.so_n_strong = zeros(size(thr));  D.so_dim = D.so_n_strong;
        D.so_gen_min_eig = D.so_n_strong;
        for k = 1:numel(thr)
            strong = act & fn > thr(k);
            Z = null(full([Geq'; Gc(:, strong)']));
            Hz = Z'*H*Z;  Gz = full(Z'*Gm*Z);
            ev = eig((Hz + Hz')/2, (Gz + Gz')/2, 'chol');
            D.so_n_strong(k) = nnz(strong);  D.so_dim(k) = size(Z,2);  D.so_gen_min_eig(k) = min(ev);
            if k == 1
                D.orth_err = norm(Z'*Z - eye(size(Z,2)), 'fro');
                D.gen_min_eig = min(ev);  D.gen_max_eig = max(ev);
                say('    second order in H^1: smallest eigenvalue %.4e N/mm, largest %.3e N/mm (||Z''Z - I|| = %.1e)\n', ...
                    D.gen_min_eig, D.gen_max_eig, D.orth_err);
            end
            say('      force threshold %.0e N: %d strongly active contacts, subspace dimension %d, smallest H^1 eigenvalue %.4e N/mm\n', ...
                thr(k), D.so_n_strong(k), D.so_dim(k), D.so_gen_min_eig(k));
        end
    end
end

function Gm = h1gram(prob)
% Gram matrix of h*sum|y|^2 + sum|dy|^2/h over the free nodes of vessel and wire.
h = prob.h0;  n = prob.n;
Gm = sparse(n, n);
for part = 1:2
    if part == 1, cols = prob.vcol; else, cols = prob.wcol; end
    K = size(cols, 1);
    Dd = spdiags([-ones(K-1,1), ones(K-1,1)], [0 1], K-1, K);
    Mloc = h*speye(K) + (Dd'*Dd)/h;
    for c = 1:2
        free = cols(:,c) > 0;  idx = cols(free, c);
        Gm(idx, idx) = Gm(idx, idx) + Mloc(free, free); %#ok<SPRIX>
    end
end
end
