function A = audit_solution(file)
%AUDIT_SOLUTION  Independent checks of a saved equilibrium (results/*.mat).
%
%   A = AUDIT_SOLUTION('results/base.mat') rebuilds the problem from the
%   saved data and checks, without re-solving:
%     1. feasibility: lumen constraint and edge lengths at the solution,
%     2. that the 25 mm box bounds on the node positions are inactive,
%     3. the signs of Corollary 3.4.1 on every contact arc: on inner-wall
%        contact the curvature toward the wall is positive and concave, on
%        outer-wall contact it is negative, and the line load with H = 0 is
%        nonnegative at nodes inside an arc,
%     4. the force balance: wall reaction from the multipliers against the
%        line load with H = 0.
%   It prints a short report and returns the numbers in a struct.

S = load(file);
G = evar_geometry(1.6);
prob = build_problem(G, S.ell, S.B, S.k, S.Pw);
Pv = S.Pv;  Pw = S.Pw;  M = size(Pw,1);  h = S.h;

% 1. Feasibility.
c = contact_terms(Pw(2:end,:), Pv, prob.R, prob.rw, prob.wcol(2:end,:), prob.vcol, prob.n);
edge = sqrt(sum(diff(Pw).^2, 2));
A.max_lumen_violation = max(c);                          % mm^2, <= 0 is feasible
A.max_edge_error = max(abs(edge - h));                   % mm

% 2. Box bounds (25 mm about the starting point, which is State 1 for the vessel).
A.max_vessel_displacement = max(S.U);                    % mm, must be < 25
sP = [0; cumsum(sqrt(sum(diff(G.P).^2, 2)))];  q = (0:M-1)'*h;
W0 = [interp1(sP, G.P(:,1), q), interp1(sP, G.P(:,2), q)];   % wire start: State 1 centerline
A.max_wire_move_from_centerline = max(sqrt(sum((Pw - W0).^2, 2)));

% 3. Signs of Corollary 3.4.1 at the contact nodes.
J = S.contact_nodes(:);                                  % wire node numbers j >= 1
kap = S.kap(:);                                          % nodes 0..M-1
kpp = zeros(M,1);  kpp(3:end-2) = (kap(4:end-1) - 2*kap(3:end-2) + kap(2:end-3))/h^2;
sig = zeros(numel(J),1);
for q = 1:numel(J)
    j = J(q);  y = Pw(j+1,:);
    [foot] = nearest_point(y, Pv);
    nh = (y - foot)/norm(y - foot);                      % outward normal of the wall
    t = Pw(min(j+2,M),:) - Pw(j,:);  t = t/norm(t);
    sig(q) = sign(nh(1)*(-t(2)) + nh(2)*t(1));           % +1 if the wall is on the left
end
ktw  = sig.*kap(J+1);                                    % curvature toward the wall
kppw = sig.*kpp(J+1);
lam  = -S.B*(kppw + 0.5*ktw.^3);                         % line load with H = 0
inside = ismember(J-1, J) & ismember(J+1, J) & ismember(J-2, J) & ismember(J+2, J);
A.arcs = S.arcs;
for a = 1:numel(S.arcs)
    r = J*h >= S.arcs(a).s0 - 1e-9 & J*h <= S.arcs(a).s1 + 1e-9;
    A.arc(a).side = S.arcs(a).side;
    A.arc(a).curvature_toward_wall = [min(ktw(r)), max(ktw(r))];
    A.arc(a).kpp_toward_wall_inside = [min(kppw(r & inside)), max(kppw(r & inside))];
    A.arc(a).lineload_inside = [min(lam(r & inside)), max(lam(r & inside))];
end

% 4. Force balance.
A.F_kkt = S.F_kkt;  A.F_line = S.F_line;  A.rel_diff = abs(S.F_line - S.F_kkt)/S.F_kkt;
A.H_absmax_after_first_edge = max(abs(S.He(2:end)));

fprintf('%s\n', file);
fprintf('  max lumen violation %.2e mm^2, max edge-length error %.2e mm\n', A.max_lumen_violation, A.max_edge_error);
fprintf('  max vessel displacement %.3f mm, max wire node move from the State 1 centerline %.3f mm (box bound 25 mm)\n', ...
    A.max_vessel_displacement, A.max_wire_move_from_centerline);
for a = 1:numel(A.arc)
    fprintf('  arc %d (%s): kappa toward wall %s 1/mm, kappa'''' inside %s, line load inside %s N/mm\n', ...
        a, A.arc(a).side, span(A.arc(a).curvature_toward_wall, '%.4f'), ...
        span(A.arc(a).kpp_toward_wall_inside, '%.2e'), span(A.arc(a).lineload_inside, '%.3f'));
end
fprintf('  wall reaction %.4f N (multipliers) vs %.4f N (line load, H = 0), %.2f%%\n', A.F_kkt, A.F_line, 100*A.rel_diff);
fprintf('  |H| after the first edge < %.4f N\n', A.H_absmax_after_first_edge);
end

function str = span(v, fmt)
% Interval [min, max] as text, or a note when the arc has too few nodes.
if numel(v) < 2, str = '(arc too short)'; else, str = sprintf(['[' fmt ', ' fmt ']'], v(1), v(2)); end
end

function f = nearest_point(y, V)
best = Inf;  f = V(1,:);
for i = 1:size(V,1)-1
    e = V(i+1,:) - V(i,:);  t = min(max(dot(y - V(i,:), e)/dot(e,e), 0), 1);
    p = V(i,:) + t*e;  d = norm(y - p);
    if d < best, best = d; f = p; end
end
end
