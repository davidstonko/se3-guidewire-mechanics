function C = clock_shift3d(P1, P2, a, s_targets, hv)
%CLOCK_SHIFT3D  Clock position of a twist-free wall target before and after
%   the wire is inserted.
%
%   C = CLOCK_SHIFT3D(P1, P2, a, s_targets, hv) takes the vessel centerline in
%   State 1 (P1) and State 2 (P2), both N-by-3 with node i the same material
%   point, a fixed unit vector a ("anterior"), State 1 arc lengths of the
%   targets, and the State 1 node spacing hv.
%
%   Frames, one per edge of the polygon (tangent T_e):
%     anterior frame  u2_e = normalized projection of a on the plane normal
%                     to T_e (12 o'clock), u1_e = T_e x u2_e (3 o'clock);
%     Bishop frame    d_e, obtained by parallel transport from d_1 = u2_1 at
%                     the first (held) vessel edge: d_{e+1} is d_e rotated by
%                     the smallest rotation taking T_e to T_{e+1}.
%   psi_e is the angle from d_e to u2_e, right-handed about T_e. A wall
%   target fixed in the Bishop frame (no vessel twist), at 12 o'clock at the
%   first node, is then at clock angle theta = -psi measured from 12 o'clock
%   toward 3 o'clock, i.e. clockwise for a viewer looking downstream along
%   T (from the sheath toward the aorta). Clock hours = theta/30.
%
%   Cross-check. Along a smooth curve psi' = (a.t)(a.(t x t'))/(1 - (a.t)^2).
%   This rate is summed over the vertices (t the normalized mean of the two
%   edge tangents, t' ds the change of tangent across the vertex), from the
%   first edge to the target, independently of the transported frame.
%
%   Targets are located by material node index s/hv (the same material
%   point in both states) and the edge values are interpolated linearly
%   between edge midpoints.
%
%   Output C (one row per target in each field):
%     s, theta1, theta2 (deg), shift (deg) = theta2 - theta1, shift_min
%     (clock minutes, 1 hour = 30 deg), shift_int (deg, from the rate
%     integrals), agree (deg) = shift - shift_int, clock1/clock2 (text),
%     and max_at1, max_at2: max |a.t| over each polygon.

a = a(:)'/norm(a);
F1 = frames(P1, a);  F2 = frames(P2, a);
xi = s_targets(:)/hv;                                 % fractional node index (0-based)
em = (1:numel(F1.psi))' - 0.5;                        % edge midpoints in node index
psi1 = interp1(em, F1.psi, xi, 'linear', 'extrap');  psi2 = interp1(em, F2.psi, xi, 'linear', 'extrap');
I1   = interp1(em, F1.I,   xi, 'linear', 'extrap');  I2   = interp1(em, F2.I,   xi, 'linear', 'extrap');
C.s = s_targets(:);
C.theta1 = -rad2deg(psi1);  C.theta2 = -rad2deg(psi2);
C.shift = C.theta2 - C.theta1;
C.shift_min = 2*C.shift;                              % 30 deg per hour -> 0.5 deg per clock minute
C.shift_int = -rad2deg(I2 - I1);
C.agree = C.shift - C.shift_int;
C.int1 = -rad2deg(I1);  C.int2 = -rad2deg(I2);       % integral-based clock angles
C.clock1 = arrayfun(@clock_text, C.theta1, 'UniformOutput', false);
C.clock2 = arrayfun(@clock_text, C.theta2, 'UniformOutput', false);
C.max_at1 = F1.max_at;  C.max_at2 = F2.max_at;
C.a = a;
end

function F = frames(P, a)
e = diff(P);  L = sqrt(sum(e.^2, 2));  T = e./L;  Ne = size(T,1);
at = T*a';
u2 = a - at.*T;  u2 = u2./sqrt(sum(u2.^2, 2));       % anterior axis on each edge
d = zeros(Ne,3);  d(1,:) = u2(1,:);
I = zeros(Ne,1);
for k = 2:Ne
    t1 = T(k-1,:);  t2 = T(k,:);  v = d(k-1,:);
    % Parallel transport: smallest rotation taking t1 to t2, applied to v.
    w = v - (dot(t1 + t2, v)/(1 + dot(t1, t2)))*(t1 + t2) + 2*dot(t1, v)*t2;
    w = w - dot(w, t2)*t2;  d(k,:) = w/norm(w);
    % Rate integral across the vertex between edges k-1 and k.
    tm = (t1 + t2)/norm(t1 + t2);  dt = t2 - t1;  am = dot(a, tm);
    I(k) = I(k-1) + am*dot(a, cross(tm, dt))/(1 - am^2);
end
psi = atan2(sum(cross(d, u2, 2).*T, 2), sum(d.*u2, 2));
F.psi = unwrap(psi);  F.I = I;  F.max_at = max(abs(at));
end

function str = clock_text(theta)
% Clock position of an angle theta (deg) measured clockwise from 12 o'clock.
mins = round(mod(theta, 360)*2);                      % clock minutes from 12:00, 0..719
if mins >= 720, mins = mins - 720; end
hh = floor(mins/60);  mm = mins - 60*hh;
if hh == 0, hh = 12; end
str = sprintf('%d:%02d', hh, mm);
end
