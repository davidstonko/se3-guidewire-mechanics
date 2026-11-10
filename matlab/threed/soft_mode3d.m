function R = soft_mode3d(A, Rconst)
%SOFT_MODE3D  Locate the softest admissible modes of a computed 3D equilibrium.
%
%   R = SOFT_MODE3D(A) solves the base case of RUN_DEMOS_3D for the iliac
%   lift A (mm, default 30; A = 0 is the planar geometry), in the rigid and
%   in the deformable vessel, and for the five smallest eigenvalues of the
%   Hessian of the Lagrangian on the null space of the equalities and the
%   strongly active lumen constraints it reports:
%     - the eigenvalue in the Euclidean metric and the smallest eigenvalue
%       in the discrete H^1 metric h*sum|y|^2 + sum|dy|^2/h (N/mm),
%     - the share of the mode on the wire, and of the wire part on nodes in
%       contact, on the free tail beyond the last contact, and in the aorta,
%     - the split of the wire part into axial (along the wire), radial
%       (along the line to the nearest centerline point) and circumferential
%       (their cross product) components,
%     - the Rayleigh quotient split into the elastic energy and the
%       curvature of the active constraints, sum_j lambda_j Hess c_j.
%   A mode that slides the wire around the wall has a circumferential
%   share near 1 and a stiffness that comes from the constraint term.
%
%   SOFT_MODE3D(A, Rconst) uses a constant lumen radius Rconst (mm).
%
%   Output goes to the command window and to results/soft_mode3d_A<A>.log
%   (soft_mode3d_A<A>_R<Rconst>.log with a constant radius).

if nargin < 1 || isempty(A), A = 30; end
if nargin < 2, Rconst = []; end
tag = sprintf('soft_mode3d_A%g', A);  if ~isempty(Rconst), tag = sprintf('%s_R%g', tag, Rconst); end
here = fileparts(mfilename('fullpath'));
old = cd(here);  restore = onCleanup(@() cd(old));
logf = fopen(fullfile(here, 'results', [tag '.log']), 'w');
closer = onCleanup(@() fclose(logf));
h = 1.6;  kbar = 0.0158;  B = 158.4e3*pi*0.889^4/64;
say('%s   MATLAB %s   soft modes, A = %g mm\n', char(datetime('now')), version, A);
G = geometry3d(h, A, Rconst);
[x, out, prob, ~, rigid] = continuation_solve3d(G, G.S_il + 120, B, kbar);
say('Rigid vessel\n');       R.rigid = modes(rigid.prob, rigid.x, rigid.out);
say('Deformable vessel\n');  R.deformable = modes(prob, x, out);
say('done\n');

    function say(varargin)
        fprintf(varargin{:});  fprintf(logf, varargin{:});
    end

    function D = modes(prob, x, out)
        [Pv, Pw] = unpack_nodes3d(prob, x);
        [~, ~, Gc, Geq] = constraints3d(x, prob);
        n = prob.n;  M = prob.M;
        [~, ~, Hc, info] = contact_terms3d(Pw(2:end,:), Pv, prob.R, prob.rw, ...
            prob.wcol(2:end,:), prob.vcol, n, out.lambda.ineqnonlin);
        fn  = 2*info.delta.*abs(out.lambda.ineqnonlin(:));
        gap = info.delta - info.d;
        act = gap <= 1e-3*info.delta;  strong = act & fn > 1e-6;
        H  = full(lagrangian_hessian3d(x, out.lambda, prob));  H = (H + H')/2;
        Hc = full(Hc);  Hc = (Hc + Hc')/2;
        [~, ~, He] = objective3d(x, prob);  He = full(He);  He = (He + He')/2;
        Z  = null(full([Geq'; Gc(:, strong)']));
        Hz = Z'*H*Z;  Hz = (Hz + Hz')/2;
        [V, ev] = eig(Hz, 'vector');  [ev, order] = sort(ev);  V = V(:, order);
        Gm = h1gram(prob);  Gz = full(Z'*Gm*Z);
        eh = eig(Hz, (Gz + Gz')/2, 'chol');
        D.min_eig = ev(1);  D.min_eig_H1 = min(eh);

        % Local frame at the wire nodes 2..M: tangent, radial, circumferential.
        Y = Pw(2:end,:);  m = M - 1;
        T = zeros(m, 3);
        T(1:m-1,:) = Pw(3:end,:) - Pw(1:end-2,:);  T(m,:) = Pw(end,:) - Pw(end-1,:);
        T = T./sqrt(sum(T.^2, 2));
        e = diff(Pv);
        F = Pv(info.seg,:) + info.t.*e(info.seg,:);
        Nr = Y - F;  dn = sqrt(sum(Nr.^2, 2));
        Nr = Nr - sum(Nr.*T, 2).*T;  Nr = Nr./max(sqrt(sum(Nr.^2, 2)), eps);
        Cc = cross(T, Nr, 2);
        s  = (1:m)'*prob.h0;
        last = max(s(act));
        s_ao = prob.S_il*prob.ell/(prob.S_il + 120);   % wire arc length at the end of the iliac, approximately
        say('  %d unknowns, null space dimension %d, %d strongly active contacts, last contact s = %.1f mm, wire length %.1f mm\n', ...
            n, size(Z,2), nnz(strong), last, prob.ell);
        say('  smallest eigenvalue: Euclidean %.3e N/mm, H^1 %.3e N/mm\n', D.min_eig, D.min_eig_H1);
        say('  mode  eigenvalue   wire   | of the wire part: contact  tail   aorta | axial  radial  circ. | near-centerline | elastic    constraint\n');
        for q = 1:5
            v  = Z*V(:, q);
            vw = reshape(v(prob.nv+1:end), 3, [])';  aw2 = sum(vw.^2, 2);  W = sum(aw2);
            share_wire = W/sum(v.^2);
            ax = sum(vw.*T, 2).^2;  ra = sum(vw.*Nr, 2).^2;  ci = sum(vw.*Cc, 2).^2;
            ok = dn > 0.2;      % the radial direction is undefined on the centerline
            D.mode(q) = struct('eig', ev(q), 'wire', share_wire, ...
                'contact', sum(aw2(act))/W, 'tail', sum(aw2(s > last))/W, 'aorta', sum(aw2(s > s_ao))/W, ...
                'axial', sum(ax)/W, 'radial', sum(ra(ok))/W, 'circ', sum(ci(ok))/W, 'undefined', sum(aw2(~ok))/W, ...
                'elastic', v'*He*v, 'constraint', v'*Hc*v);
            c = D.mode(q);
            say('  %d   %.3e   %5.3f  | %5.3f   %5.3f  %5.3f | %5.3f  %5.3f  %5.3f | %5.3f | %+.3e %+.3e\n', ...
                q, c.eig, c.wire, c.contact, c.tail, c.aorta, c.axial, c.radial, c.circ, c.undefined, c.elastic, c.constraint);
        end
        % The softest mode along the wire, every 10th node.
        v = Z*V(:, 1);  vw = reshape(v(prob.nv+1:end), 3, [])';
        say('  softest mode along the wire (s, |v|, circumferential part, contact force):\n');
        for j = 10:10:m
            say('    s = %6.1f mm   |v| %.4f   circ %+.4f   radial %+.4f   axial %+.4f   gap %.3f mm   force %.4f N\n', ...
                s(j), norm(vw(j,:)), dot(vw(j,:), Cc(j,:)), dot(vw(j,:), Nr(j,:)), dot(vw(j,:), T(j,:)), gap(j), fn(j));
        end
        D.aorta_contact_force = sum(fn(act & s > s_ao));
        D.mean_force_over_delta = mean(fn(act & s > s_ao)./info.delta(act & s > s_ao));
        say('  aortic contact force %.3f N on %d nodes, mean nodal force/delta there %.3e N/mm\n', ...
            D.aorta_contact_force, nnz(act & s > s_ao), D.mean_force_over_delta);
    end
end

function Gm = h1gram(prob)
% Gram matrix of h*sum|y|^2 + sum|dy|^2/h over the free nodes of vessel and wire.
h = prob.h0;  n = prob.n;
Gm = sparse(n, n);
for part = 1:2
    if part == 1
        if prob.rigid, continue, end
        cols = prob.vcol;
    else
        cols = prob.wcol;
    end
    K = size(cols, 1);
    Dd = spdiags([-ones(K-1,1), ones(K-1,1)], [0 1], K-1, K);
    Mloc = h*speye(K) + (Dd'*Dd)/h;
    for c = 1:3
        free = cols(:,c) > 0;  idx = cols(free, c);
        Gm(idx, idx) = Gm(idx, idx) + Mloc(free, free); %#ok<SPRIX>
    end
end
end
