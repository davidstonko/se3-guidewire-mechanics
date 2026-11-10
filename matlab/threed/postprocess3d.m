function S = postprocess3d(prob, x, out, kappa_min)
%POSTPROCESS3D  Forces, material force, line load and geometry of a computed
%   3D equilibrium (the 3D counterpart of the planar POSTPROCESS).
%
%   S = POSTPROCESS3D(prob, x, out) evaluates, from the solution x and the
%   fmincon output out:
%
%   Contact forces. The multiplier lambda_j of c_j = d_j^2 - delta_j^2 <= 0 is
%   a force per unit of c_j; since |grad c_j| = 2 d_j = 2 delta_j on contact,
%   the nodal contact force is f_j = 2 delta_j |lambda_j| (N). The force
%   vector on the wire is -lambda_j grad_{y_j} c_j.
%
%   Material force. The wire energy depends on the edge lengths only through
%   the Voronoi lengths, which carries the part -B kappa^2/2 of the axial
%   force. The remainder H_e = 2 h nu_e, with nu_e the multiplier of the
%   length constraint on edge e, is the discrete material force on that edge.
%
%   Line load with H = 0. For a twist-free isotropic rod the internal force is
%   F = -B r''' + (H - 3/2 B kappa^2) r', so with H = 0 the load the wall
%   applies per unit length is
%       f = B r'''' + 3/2 B (kappa^2 r')',
%   whose component along the Frenet normal n is
%       f.n = B (kappa'' + kappa^3/2 - kappa tau^2) = -lambda,
%   and along the binormal b is B (2 kappa' tau + kappa tau').
%   Two discrete versions are compared with the multipliers on the contact nodes:
%     (i)  Frenet form: lambda = -B (kappa'' + kappa^3/2 - kappa tau^2) with
%          kappa_j = phi_j/h, kappa'' by central differences, and tau from the
%          rotation of the discrete binormal between neighbouring vertices;
%          also the binormal part B (2 kappa' tau + kappa tau'). It is used
%          only where kappa > kappa_min (default 2e-3 1/mm) at the five nodes
%          j-2..j+2, since the Frenet frame is undefined at kappa = 0. The
%          normal and binormal parts are compared with the components of the
%          contact force vectors along the discrete n and b at the same nodes;
%     (ii) vector form: f_j = B D4 y_j/h^4 + 3/2 B D(kappa^2 D y)_j/h^2, with
%          D4 the fourth difference, which needs no Frenet frame and is
%          evaluated at every interior node j = 2..n_w-3.
%
%   Side of a contact: the wall normal at the foot point (from the vessel
%   axis toward the wire) is compared with the curvature vector of the
%   deformed vessel there: inner wall if they point the same way, outer if
%   opposite, 'aorta' where the vessel curvature is below 2e-3 1/mm (as in
%   the planar code). The angle between the two is also reported.

if nargin < 4 || isempty(kappa_min), kappa_min = 2e-3; end
[Pv, Pw] = unpack_nodes3d(prob, x);
[~, ~, ~, parts] = objective3d(x, prob);
h = prob.h0;  B = prob.B;  M = prob.M;  N = prob.N;
s_v = (0:N-1)'*prob.hv;                            % arc length along the State 1 vessel
s_w = (0:M-1)'*h;

% Vessel displacement.
Uv = Pv - prob.P;
U = sqrt(sum(Uv.^2, 2));
[peak, ip] = max(U);
il = s_v <= prob.S_il;
first = s_v < prob.S_il/2;

% Contact.
[c, Gc, ~, info] = contact_terms3d(Pw(2:end,:), Pv, prob.R, prob.rw, prob.wcol(2:end,:), prob.vcol, prob.n);
lam_in = out.lambda.ineqnonlin(:);
fn = 2*info.delta.*abs(lam_in);                     % nodal contact force, nodes 1..M-1
act = -c <= 1e-3*info.delta;                        % nodes in contact (planar criterion)
actn = find(act);                                   % wire node numbers j (y_j, j >= 1)
fvec = zeros(M-1, 3);                               % contact force vector on the wire
for j = 1:M-1
    fvec(j,:) = -lam_in(j)*full(Gc(prob.wcol(j+1,:), j))';
end
foot = Pv(info.seg,:) + info.t.*(Pv(info.seg+1,:) - Pv(info.seg,:));
sigma = (info.seg - 1 + info.t)*prob.hv;            % State 1 vessel arc length of each foot point

% Discrete curvature, Frenet normal and binormal at the wire nodes (0..M-1).
e  = diff(Pw);  T = e./sqrt(sum(e.^2, 2));
cr = cross(T(1:end-1,:), T(2:end,:), 2);
phi = atan2(sqrt(sum(cr.^2, 2)), sum(T(1:end-1,:).*T(2:end,:), 2));
kap = [0; phi/h; 0];
bn = zeros(M,3);  nn = zeros(M,3);
bn(2:M-1,:) = cr./max(sqrt(sum(cr.^2, 2)), eps);
dT = T(2:end,:) - T(1:end-1,:);
nn(2:M-1,:) = dT./max(sqrt(sum(dT.^2, 2)), eps);   % Frenet normal (toward the centre of curvature)
% Torsion: rotation of the binormal about the edge tangent between neighbouring vertices.
taue = nan(M-1,1);                                  % edge e joins nodes e-1 and e (1-based e)
for q = 2:M-2                                       % edges between two interior nodes
    b1 = bn(q,:);  b2 = bn(q+1,:);
    taue(q) = atan2(dot(cross(b1, b2), T(q,:)), dot(b1, b2))/h;
end
tau = nan(M,1);
tau(3:M-2) = 0.5*(taue(2:M-3) + taue(3:M-2));
kpp = zeros(M,1);  kpp(3:end-2) = (kap(4:end-1) - 2*kap(3:end-2) + kap(2:end-3))/h^2;
kp = nan(M,1);   kp(2:M-1) = (kap(3:M) - kap(1:M-2))/(2*h);     % kappa'
taup = nan(M,1); taup(4:M-3) = (tau(5:M-2) - tau(3:M-4))/(2*h); % tau'
valid = false(M,1);
for q = 4:M-3
    valid(q) = all(kap(q-2:q+2) > kappa_min) && ~isnan(taup(q));
end
lam_frenet = nan(M,1);  fb_frenet = nan(M,1);
lam_frenet(valid) = -B*(kpp(valid) + 0.5*kap(valid).^3 - kap(valid).*tau(valid).^2);   % -(f.n)
fb_frenet(valid)  =  B*(2*kp(valid).*tau(valid) + kap(valid).*taup(valid));          %   f.b

% Vector form of the line load (no Frenet frame needed).
fline = nan(M,3);
k2e = 0.5*(kap(1:end-1).^2 + kap(2:end).^2);        % kappa^2 on each edge
for q = 3:M-2                                       % node q is y_{q-1}
    D4 = Pw(q+2,:) - 4*Pw(q+1,:) + 6*Pw(q,:) - 4*Pw(q-1,:) + Pw(q-2,:);
    fline(q,:) = B*D4/h^4 + 1.5*B*(k2e(q)*e(q,:) - k2e(q-1)*e(q-1,:))/h^2;
end
fline_mag = sqrt(sum(fline.^2, 2));

% Force balance on the contact nodes.
inc = false(M,1);  inc(actn + 1) = true;  inc(end) = false;
incv = inc & ~isnan(fline_mag);
F_vec = sum(fline_mag(incv))*h;
incf = inc & valid;
F_frenet = sum(abs(lam_frenet(incf)))*h;
F_frenet_full = sum(sqrt(lam_frenet(incf).^2 + fb_frenet(incf).^2))*h;
fn0 = [0; fn];                                      % padded to nodes 0..M-1
fv0 = [0 0 0; fvec];
F_kkt_frenet_nodes = sum(fn0(incf));
F_kkt_normal = sum(abs(sum(fv0(incf,:).*nn(incf,:), 2)));    % normal part of the contact forces
F_kkt_binormal = sum(abs(sum(fv0(incf,:).*bn(incf,:), 2)));  % binormal part
F_frenet_binormal = sum(abs(fb_frenet(incf)))*h;
% Direction check: angle between the contact force and the line-load vector on contact nodes.
ang = nan(M,1);
for q = find(incv & fn0 > 1e-6)'
    ang(q) = acosd(max(-1, min(1, dot(fv0(q,:), fline(q,:))/(norm(fv0(q,:))*fline_mag(q)))));
end
% Binormal share of the line load on contact nodes with a defined Frenet frame.
fb_share = nan(M,1);
for q = find(incf)'
    fb_share(q) = abs(dot(fline(q,:), bn(q,:)))/fline_mag(q);
end

% Material force on each edge, and the axial force N = H - B kappa^2/2.
nu = out.lambda.eqnonlin(1:M-1);
He = 2*h*nu;
kmid = 0.5*(kap(1:end-1) + kap(2:end));
halfBk2 = 0.5*B*kmid.^2;
N_edge = He - halfBk2;

% Side of each contact against the curvature vector of the deformed vessel.
kv = curvature_vectors(Pv);
side = cell(numel(actn), 1);  angk = nan(numel(actn), 1);
for q = 1:numel(actn)
    j = actn(q);  s = info.seg(j);  t = info.t(j);
    kf = (1 - t)*kv(s,:) + t*kv(s+1,:);            % vessel curvature vector at the foot point
    w = Pw(j+1,:) - foot(j,:);  nh = w/norm(w);    % wall normal, axis -> wire
    if norm(kf) < 2e-3
        side{q} = 'aorta';
    elseif dot(nh, kf) > 0
        side{q} = 'inner';
    else
        side{q} = 'outer';
    end
    if norm(kf) > 0, angk(q) = acosd(max(-1, min(1, dot(nh, kf)/norm(kf)))); end
end
% Contact arcs: runs of consecutive contact nodes.
actn = actn(:);
brk = [0; find(diff(actn) > 1); numel(actn)];
arcs = struct('s0', {}, 's1', {}, 'sig0', {}, 'sig1', {}, 'side', {}, 'angle', {}, 'force', {}, 'nodes', {});
for q = 1:numel(brk)-1
    r = (brk(q)+1:brk(q+1))';
    if isempty(r), continue; end
    jr = actn(r);
    wts = fn(jr);
    if sum(wts) > 0, am = sum(wts.*angk(r), 'omitnan')/sum(wts(~isnan(angk(r)))); else, am = mean(angk(r), 'omitnan'); end
    arcs(end+1) = struct('s0', jr(1)*h, 's1', jr(end)*h, 'sig0', min(sigma(jr)), 'sig1', max(sigma(jr)), ...
        'side', mode_char(side(r)), 'angle', am, 'force', sum(fn(jr)), 'nodes', numel(jr)); %#ok<AGROW>
end
% Force on each arc from the line load (vector form).
for q = 1:numel(arcs)
    in = false(M,1);  in(round(arcs(q).s0/h)+1 : round(arcs(q).s1/h)+1) = true;
    arcs(q).force_line = sum(fline_mag(in & incv))*h;
end

S = struct();
S.flag = out.flag;  S.firstorderopt = out.firstorderopt;  S.time = out.time;
S.iterations = out.iterations;  S.constrviolation = out.constrviolation;
S.h = h;  S.ell = prob.ell;  S.M = M;  S.N = N;  S.B = B;  S.k = prob.k(1);  S.rigid = prob.rigid;
S.E = sum(parts);  S.E_wire = parts(1);  S.E_vessel = parts(2);  S.E_foundation = parts(3);
S.peak = peak;  S.peak_at_s = s_v(ip);  S.peak_node = ip;  S.peak_vector = Uv(ip,:);
S.peak_in_iliac = il(ip);
S.peak_first_bend = max(U(first));  S.aorta_shift = max(U(~il));  S.rms = sqrt(mean(U.^2));
S.F_kkt = sum(fn);  S.F_line_vec = F_vec;  S.F_line_frenet = F_frenet;
S.F_kkt_on_frenet_nodes = F_kkt_frenet_nodes;
S.F_frenet_full = F_frenet_full;  S.F_frenet_binormal = F_frenet_binormal;
S.F_kkt_normal = F_kkt_normal;    S.F_kkt_binormal = F_kkt_binormal;
S.n_contact = numel(actn);  S.n_frenet_valid_contact = nnz(incf);  S.n_vec_contact = nnz(incv);
S.kappa_min = kappa_min;
S.max_angle_force_vs_lineload = max(ang, [], 'omitnan');
S.median_angle_force_vs_lineload = median(ang, 'omitnan');
S.max_binormal_share = max(fb_share, [], 'omitnan');
S.tip_force = fn(end);
S.H_absmax_after_first_edge = max(abs(He(2:end)));  S.H_first_edge = He(1);
S.halfBk2_max = max(halfBk2);  S.kappa_exit = kap(2);  S.halfBk2_exit = 0.5*B*kap(2)^2;
S.s_last_contact = max(actn)*h;
S.arcs = arcs;
S.Pv = Pv;  S.Pw = Pw;  S.P = prob.P;  S.R = prob.R;  S.S_il = prob.S_il;  S.U = U;  S.Uv = Uv;
S.s_v = s_v;  S.s_w = s_w;  S.fn = fn;  S.fvec = fvec;  S.contact_nodes = actn;  S.sigma = sigma;
S.kap = kap;  S.tau = tau;  S.kpp = kpp;  S.lam_frenet = lam_frenet;  S.fb_frenet = fb_frenet;  S.fline = fline;
S.fline_mag = fline_mag;  S.valid_frenet = valid;  S.angle_force_vs_lineload = ang;
S.cum_kkt = cumsum(fn0);
S.cum_line_vec = cumsum(fline_mag.*incv*h, 'omitnan');
lf = abs(lam_frenet);  lf(~incf) = 0;
S.cum_line_frenet = cumsum(lf*h);
S.cum_kkt_frenet_nodes = cumsum(fn0.*incf);
S.He = He;  S.N_edge = N_edge;  S.halfBk2 = halfBk2;
end

function k = curvature_vectors(Q)
% Discrete curvature vector at each node: change of unit tangent per Voronoi length.
e = diff(Q);  L = sqrt(sum(e.^2, 2));  t = e./L;
k = zeros(size(Q));
k(2:end-1,:) = (t(2:end,:) - t(1:end-1,:))./(0.5*(L(2:end) + L(1:end-1)));
end

function m = mode_char(s)
[u, ~, j] = unique(s);  m = u{mode(j)};
end
