function S = postprocess(prob, x, out)
%POSTPROCESS  Forces, material force and geometry of a computed equilibrium.
%
%   S = POSTPROCESS(prob, x, out) evaluates the quantities reported in
%   Section 4 of the paper from the solution x and the fmincon output.
%
%   Contact forces. The multiplier lambda_j of c_j = d_j^2 - delta_j^2 <= 0
%   is a force per unit of c_j. Since |grad c_j| = 2 d_j = 2 delta_j on
%   contact, the nodal contact force is f_j = 2 delta_j lambda_j (N).
%
%   Line load. With H = 0 (Theorem 3.3.1) the wall load is
%   lambda = -B (kappa'' + kappa^3/2) (Corollary 3.4.1), evaluated with the
%   discrete curvature kappa_j = phi_j/h and central differences.
%
%   Material force. The wire energy depends on the edge lengths through the
%   Voronoi lengths, which carries the part -B kappa^2/2 of the axial force.
%   The remainder, 2 h nu_e with nu_e the multiplier of the length constraint
%   on edge e, is the discrete material force H on that edge.

[Pv, Pw] = unpack_nodes(prob, x);
[~, ~, ~, parts] = objective(x, prob);
h = prob.h0;  B = prob.B;  M = prob.M;
s_v = (0:prob.N-1)'*h;                             % arc length along the State 1 vessel
s_w = (0:M-1)'*h;

% Vessel displacement.
U = sqrt(sum((Pv - prob.P).^2, 2));
[peak, ip] = max(U);
il = s_v <= prob.S_il;
first = s_v < prob.S_il/2;

% Contact.
[c, ~, ~, info] = contact_terms(Pw(2:end,:), Pv, prob.R, prob.rw, prob.wcol(2:end,:), prob.vcol, prob.n);
gapv = -c;                                          % >= 0 inside the lumen
lam_in = out.lambda.ineqnonlin;
fn = 2*info.delta.*abs(lam_in);                     % nodal contact force, nodes 1..M-1
act = gapv <= 1e-3*info.delta;                      % nodes in contact
actn = find(act);                                   % wire node numbers j (y_j, j >= 1)

% Discrete curvature at interior wire nodes, padded with zeros at the sheath and tip.
e = diff(Pw);
phi = atan2(e(1:end-1,1).*e(2:end,2) - e(1:end-1,2).*e(2:end,1), sum(e(1:end-1,:).*e(2:end,:), 2));
kap = [0; phi/h; 0];                                % nodes 0..M-1
kpp = zeros(M,1);  kpp(3:end-2) = (kap(4:end-1) - 2*kap(3:end-2) + kap(2:end-3))/h^2;
lineload = -B*(kpp + 0.5*kap.^3);                   % nodes 0..M-1
inc = false(M,1);  inc(actn + 1) = true;  inc(end) = false;
F_line = sum(abs(lineload(inc)))*h;

% Material force on each edge, and the axial force N = H - B kappa^2/2.
nu = out.lambda.eqnonlin(1:M-1);
He = 2*h*nu;
kmid = 0.5*(kap(1:end-1) + kap(2:end));
halfBk2 = 0.5*B*kmid.^2;
N_edge = He - halfBk2;

% Side of each contact: inner wall of a bend, outer wall, or straight aorta.
[kv, ~] = curvature_vectors(Pv);
side = strings(numel(actn), 1);
for q = 1:numel(actn)
    yj = Pw(actn(q)+1,:);
    [~, i] = min(sum((Pv - yj).^2, 2));
    nh = (yj - Pv(i,:))/norm(yj - Pv(i,:));
    if norm(kv(i,:)) < 2e-3, side(q) = "aorta";
    elseif dot(nh, kv(i,:)) > 0, side(q) = "inner";
    else, side(q) = "outer";
    end
end
% Contact arcs: runs of consecutive contact nodes.
actn = actn(:);
brk = [0; find(diff(actn) > 1); numel(actn)];
arcs = struct('s0', {}, 's1', {}, 'side', {}, 'force', {});
for q = 1:numel(brk)-1
    r = (brk(q)+1:brk(q+1))';
    if isempty(r), continue; end
    arcs(end+1) = struct('s0', actn(r(1))*h, 's1', actn(r(end))*h, ...
        'side', mode_string(side(r)), 'force', sum(fn(actn(r)))); %#ok<AGROW>
end

% Straight continuation of the wire beyond the last contact.
lc = max(actn);
t = Pw(end,:) - Pw(lc+1,:);  t = t/norm(t);
ss = (0:0.05:200)';
pts = Pw(lc+1,:) + ss.*t;
cc = contact_terms(pts, Pv, prob.R, prob.rw, zeros(numel(ss),2), prob.vcol, prob.n);
bad = find(cc > 0, 1);
if isempty(bad), s_free = Inf; else, s_free = lc*h + ss(bad); end

S = struct();
S.flag = out.flag;  S.firstorderopt = out.firstorderopt;  S.time = out.time;
S.iterations = out.iterations;
S.h = h;  S.ell = prob.ell;  S.M = M;  S.N = prob.N;  S.B = B;  S.k = prob.k(1);
S.E = sum(parts);  S.E_wire = parts(1);  S.E_vessel = parts(2);  S.E_foundation = parts(3);
S.peak = peak;  S.peak_at_s = s_v(ip);  S.peak_first_bend = max(U(first));
S.aorta_shift = max(U(~il));  S.rms = sqrt(mean(U.^2));
S.F_kkt = sum(fn);  S.F_line = F_line;  S.tip_force = fn(end);
S.H_absmax_after_first_edge = max(abs(He(2:end)));  S.H_first_edge = He(1);
S.halfBk2_max = max(halfBk2);  S.kappa_exit = kap(2);  S.halfBk2_exit = 0.5*B*kap(2)^2;
S.s_last_contact = lc*h;  S.s_free_until = s_free;
S.arcs = arcs;
S.Pv = Pv;  S.Pw = Pw;  S.P = prob.P;  S.R = prob.R;  S.S_il = prob.S_il;  S.U = U;
S.s_v = s_v;  S.s_w = s_w;  S.fn = fn;  S.contact_nodes = actn;  S.kap = kap;
S.cum_kkt = cumsum([0; fn]);  S.cum_line = cumsum(abs(lineload).*inc*h);
S.He = He;  S.N_edge = N_edge;  S.halfBk2 = halfBk2;
end

function [k, t] = curvature_vectors(Q)
% Discrete curvature vector at each node: change of unit tangent per Voronoi length.
e = diff(Q);  L = sqrt(sum(e.^2, 2));  t = e./L;
k = zeros(size(Q));
k(2:end-1,:) = (t(2:end,:) - t(1:end-1,:))./(0.5*(L(2:end) + L(1:end-1)));
end

function m = mode_string(s)
[u, ~, j] = unique(s);  m = u(mode(j));
end
