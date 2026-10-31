function prob = build_problem(G, ell, B, k, Pw_init)
%BUILD_PROBLEM  Discrete wire-in-vessel problem (free tip, fixed length).
%
%   prob = BUILD_PROBLEM(G, ell, B, k) sets up the nonlinear program of
%   Section 2.4 of the paper for the geometry G from EVAR_GEOMETRY:
%
%     minimize   E(x) = E_wire + E_vessel + E_foundation
%     subject to c_j(x)  = dist(y_j, Gamma_h)^2 - delta_j^2 <= 0,  j = 1..n_w-1
%                ceq_j(x) = |y_j - y_{j-1}|^2 - h_w^2       = 0,  j = 1..n_w-1
%                ceq_0(x) = (y_1 - y_0) x t_0                = 0  (sheath direction)
%
%   Unknowns x: the interior vessel nodes x_2..x_{N-1} (the vessel ends are
%   held) and the wire nodes y_1..y_{n_w-1} (y_0 is held at the sheath).
%   Each node contributes two coordinates, stored as [x1 y1 x2 y2 ...].
%
%   Inputs
%     G    geometry struct (P, R, S_il, hv)
%     ell  length of wire in the vessel (mm), rounded to a whole number of edges
%     B    wire bending stiffness (N mm^2)
%     k    foundation stiffness per node (N/mm), scalar or N-by-1
%     Pw_init (optional) initial wire nodes, n_w-by-2
%
%   The wire and vessel use the same node spacing h = G.hv.

P = G.P;  N = size(P,1);  hv = G.hv;
M = round(ell/hv) + 1;                      % number of wire nodes
prob.P   = P;
prob.R   = G.R;
prob.S_il = G.S_il;
prob.N   = N;
prob.M   = M;
prob.ell = (M - 1)*hv;
prob.h0  = prob.ell/(M - 1);                % wire edge length
prob.rw  = 0.4445;                          % wire radius, 0.035 inch (mm)
prob.B   = B;
prob.Bv  = 200;                             % vessel bending stiffness (N mm^2)
prob.k   = k(:).*ones(N,1);                 % foundation stiffness per node (N/mm)
prob.phi_ref = turning_angles(P);           % State 1 turning angles of the vessel
prob.a   = P(1,:);                          % sheath position
t0 = P(2,:) - P(1,:);  prob.t0 = t0/norm(t0);   % sheath direction

% Column of each node coordinate in x (0 = held node).
prob.nv = 2*(N - 2);                        % number of vessel unknowns
prob.nw = 2*(M - 1);                        % number of wire unknowns
prob.n  = prob.nv + prob.nw;
vcol = zeros(N,2);  vcol(2:N-1,:) = reshape(1:prob.nv, 2, [])';
wcol = zeros(M,2);  wcol(2:M,:)   = prob.nv + reshape(1:prob.nw, 2, [])';
prob.vcol = vcol;  prob.wcol = wcol;

% Initial wire: equal arc length along the State 1 centerline polygon.
if nargin < 5 || isempty(Pw_init)
    Pw_init = resample_polyline(P, prob.ell, M);
end
prob.Pw1 = Pw_init;
prob.x_ref = pack_nodes(prob, P, Pw_init);        % reference point for the box bounds
end

function phi = turning_angles(Q)
e = diff(Q);
phi = atan2(e(1:end-1,1).*e(2:end,2) - e(1:end-1,2).*e(2:end,1), sum(e(1:end-1,:).*e(2:end,:), 2));
end

function Q = resample_polyline(V, ell, M)
s = [0; cumsum(sqrt(sum(diff(V).^2, 2)))];
q = linspace(0, ell, M)';
Q = [interp1(s, V(:,1), q, 'linear', 'extrap'), interp1(s, V(:,2), q, 'linear', 'extrap')];
end
