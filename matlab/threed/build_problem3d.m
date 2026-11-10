function prob = build_problem3d(G, ell, B, k, Pw_init, rigid)
%BUILD_PROBLEM3D  Discrete 3D wire-in-vessel problem (free tip, fixed length).
%
%   prob = BUILD_PROBLEM3D(G, ell, B, k, Pw_init, rigid) sets up
%
%     minimize   E(x) = E_wire + E_vessel + E_foundation
%     subject to c_j(x)  = dist(y_j, Gamma_h)^2 - delta_j^2 <= 0,  j = 1..n_w-1
%                ceq_j(x) = |y_j - y_{j-1}|^2 - h_w^2       = 0,  j = 1..n_w-1
%                (y_1 - y_0).f1 = 0,  (y_1 - y_0).f2 = 0     (sheath direction)
%
%   for the geometry G from GEOMETRY3D. f1, f2 is an orthonormal basis of the
%   plane normal to the sheath direction t0 (the first State 1 vessel edge),
%   so the two equalities say that the first wire edge is parallel to t0.
%
%   Energies:
%     E_wire = (B/2) sum_i |kb_i|^2 / l_i     (discrete elastic rod, twist free)
%     E_vessel     = (B_v/2) sum_i |u_{i+1} - 2 u_i + u_{i-1}|^2 / h^3
%     E_foundation = (1/2)   sum_i k_i |u_i|^2
%   with u_i = x_i - P_i the vessel displacement (geometrically linear beam on a
%   Winkler foundation) and the two vessel end nodes held at u = 0.
%
%   Unknowns x: the interior vessel nodes x_2..x_{N-1} and the wire nodes
%   y_1..y_{n_w-1}, three coordinates each, stored as [x1 y1 z1 x2 ...].
%   With rigid = true the vessel is held at State 1 (B_v -> infinity) and the
%   unknowns are the wire nodes only.
%
%   Inputs
%     G       geometry struct (P, R, S_il, hv)
%     ell     length of wire in the vessel (mm), rounded to a whole number of edges
%     B       wire bending stiffness (N mm^2)
%     k       foundation stiffness per node (N/mm), scalar or N-by-1 (k_bar * h)
%     Pw_init (optional) initial wire nodes, n_w-by-3 (default: State 1 centerline)
%     rigid   (optional) hold the vessel fixed, default false

if nargin < 6 || isempty(rigid), rigid = false; end
P = G.P;  N = size(P,1);  hv = G.hv;
M = round(ell/hv) + 1;                      % number of wire nodes
prob.P    = P;
prob.R    = G.R;
prob.S_il = G.S_il;
prob.hv   = hv;
prob.N    = N;
prob.M    = M;
prob.ell  = (M - 1)*hv;
prob.h0   = prob.ell/(M - 1);               % wire edge length
prob.rw   = 0.4445;                         % wire radius, 0.035 inch (mm)
prob.B    = B;
prob.Bv   = 200;                            % vessel bending stiffness (N mm^2)
prob.k    = k(:).*ones(N,1);                % foundation stiffness per node (N/mm)
prob.rigid = rigid;
prob.a    = P(1,:);                         % sheath position
t0 = P(2,:) - P(1,:);  t0 = t0/norm(t0);
prob.t0 = t0;                               % sheath direction
% Orthonormal basis of the plane normal to t0.
ref = [0 0 1];  if abs(dot(ref, t0)) > 0.9, ref = [1 0 0]; end
f1 = cross(t0, ref);  f1 = f1/norm(f1);
f2 = cross(t0, f1);   f2 = f2/norm(f2);
prob.f1 = f1;  prob.f2 = f2;

% Column of each node coordinate in x (0 = held node).
if rigid, prob.nv = 0; else, prob.nv = 3*(N - 2); end
prob.nw = 3*(M - 1);
prob.n  = prob.nv + prob.nw;
vcol = zeros(N,3);
if ~rigid, vcol(2:N-1,:) = reshape(1:prob.nv, 3, [])'; end
wcol = zeros(M,3);  wcol(2:M,:) = prob.nv + reshape(1:prob.nw, 3, [])';
prob.vcol = vcol;  prob.wcol = wcol;

% Vessel energy is quadratic in the free displacements: 1/2 u'(Kb + Kk) u.
if ~rigid
    D2 = spdiags(ones(N-2,1)*[1 -2 1], 0:2, N-2, N);        % second difference, rows i = 2..N-1
    Kn_b = (prob.Bv/hv^3)*(D2'*D2);                          % node-level bending matrix
    Kn_b = Kn_b(2:N-1, 2:N-1);                               % ends held at u = 0
    Kn_k = spdiags(prob.k(2:N-1), 0, N-2, N-2);
    prob.Kb = kron(Kn_b, speye(3));
    prob.Kk = kron(Kn_k, speye(3));
    prob.xP = reshape(P(2:N-1,:)', [], 1);                   % State 1 interior vessel nodes
    prob.Hv = blkdiag(prob.Kb + prob.Kk, sparse(prob.nw, prob.nw));
end

% Initial wire: equal arc length along the State 1 centerline polygon.
if nargin < 5 || isempty(Pw_init)
    Pw_init = resample_polyline(P, prob.ell, M);
end
prob.Pw1 = Pw_init;
prob.x_ref = pack_nodes3d(prob, P, resample_polyline(P, prob.ell, M));   % centre of the box bounds
end

function Q = resample_polyline(V, ell, M)
s = [0; cumsum(sqrt(sum(diff(V).^2, 2)))];
q = linspace(0, ell, M)';
Q = interp1(s, V, q, 'linear', 'extrap');
end
