function G = geometry3d(h, A, Rconst)
%GEOMETRY3D  Non-planar test vessel (State 1) for the 3D wire problem.
%
%   G = GEOMETRY3D(h, A) returns the State 1 centerline sampled at uniform
%   arc length with spacing close to h (mm), the lumen radius at every node,
%   and numerical checks of the curve.
%
%   The iliac segment is the planar sine of EVAR_GEOMETRY lifted out of its
%   plane,
%       C(x) = (x, 25 sin(2 pi x/150), A sin(pi x/150)),   0 <= x <= 150 mm,
%   and the aorta continues straight along the end tangent for 200 mm. The
%   lift vanishes with its second derivative at both ends, so the centerline
%   is C^2 at the iliac-aortic junction (curvature zero on both sides) and is
%   straight at the sheath, as in the planar vessel. A = 0 gives exactly the
%   planar vessel of EVAR_GEOMETRY, embedded in the plane z = 0.
%   The lumen radius is 6 mm in the iliac segment and rises to 10 mm over the
%   first 20 mm of the aorta with the same C^1 smoothstep.
%
%   Fields of G:
%     P      N-by-3 node positions of the State 1 centerline (mm)
%     R      N-by-1 lumen radius at each node (mm)
%     S_il   arc length of the iliac segment (mm)
%     hv     actual node spacing (mm)
%     A      out-of-plane amplitude (mm)
%     check  struct of curve checks:
%            min_Rc          smallest radius of curvature of the smooth curve (mm)
%            max_abs_tau     largest |torsion| where the curvature exceeds 1e-3 1/mm
%            rms_tau         RMS torsion over the iliac arc length (same points)
%            poly_min_local_radius  smallest Voronoi length / turning angle of the node polygon
%            poly_min_far_dist      smallest distance between nodes more than pi*10 mm apart in arc length
%            tube_ok         true if the 10 mm tube about the centerline cannot overlap itself

L_ao = 200; R_il = 6; R_ao = 10; L_tr = 20; Lx = 150;
k = 2*pi/Lx;  p = pi/Lx;

% Densely sampled iliac curve and its exact derivatives in the parameter x.
x  = linspace(0, Lx, 200001)';
C  = [x, 25*sin(k*x), A*sin(p*x)];
d1 = [ones(size(x)),  25*k*cos(k*x),   A*p*cos(p*x)];
d2 = [zeros(size(x)), -25*k^2*sin(k*x), -A*p^2*sin(p*x)];
d3 = [zeros(size(x)), -25*k^3*cos(k*x), -A*p^3*cos(p*x)];
cr = cross(d1, d2, 2);
ncr = sqrt(sum(cr.^2, 2));
kappa = ncr./sqrt(sum(d1.^2, 2)).^3;               % curvature
tau = sum(cr.*d3, 2)./ncr.^2;                       % torsion (undefined where kappa = 0)
sC = [0; cumsum(sqrt(sum(diff(C).^2, 2)))];
S_il = sC(end);
T  = C(end,:) - C(end-1,:);  T = T/norm(T);        % end tangent (last chord, as in the planar code)

L = S_il + L_ao;
N = round(L/h) + 1;
s = linspace(0, L, N)';
P = zeros(N, 3);
il = s <= S_il;
for d = 1:3
    P(il,d) = interp1(sC, C(:,d), s(il));
end
P(~il,:) = C(end,:) + (s(~il) - S_il).*T;          % straight aorta

z = min(max((s - S_il)/L_tr, 0), 1);
R = R_il + (R_ao - R_il)*z.^2.*(3 - 2*z);          % C^1 smoothstep taper
if nargin > 2 && ~isempty(Rconst), R(:) = Rconst; end   % constant lumen radius (Assumption 5.1(c))

% Checks of the smooth curve.
ok = kappa > 1e-3;
dsx = sqrt(sum(d1.^2, 2));                          % ds/dx
chk.min_Rc = 1/max(kappa);
if any(ok)
    chk.max_abs_tau = max(abs(tau(ok)));
    chk.rms_tau = sqrt(trapz(x(ok), tau(ok).^2.*dsx(ok))/trapz(x(ok), dsx(ok)));
else
    chk.max_abs_tau = 0;  chk.rms_tau = 0;
end
chk.total_torsion = trapz(x(ok), tau(ok).*dsx(ok));  % integral of tau ds (mostly cancels by symmetry)

% Checks of the node polygon: local radius and distant-node separation.
e = diff(P);  Le = sqrt(sum(e.^2, 2));
cc = sum(e(1:end-1,:).*e(2:end,:), 2);
xx = sqrt(sum(cross(e(1:end-1,:), e(2:end,:), 2).^2, 2));
phi = atan2(xx, cc);
chk.poly_min_local_radius = min(0.5*(Le(1:end-1) + Le(2:end))./max(phi, eps));
sp = [0; cumsum(Le)];
far = abs(sp - sp') > pi*R_ao;
DX = sqrt((P(:,1) - P(:,1)').^2 + (P(:,2) - P(:,2)').^2 + (P(:,3) - P(:,3)').^2);
chk.poly_min_far_dist = min(DX(far));
chk.tube_ok = chk.min_Rc > R_ao && chk.poly_min_local_radius > R_ao && chk.poly_min_far_dist > 2*R_ao;

G = struct('P', P, 'R', R, 'S_il', S_il, 'hv', s(2) - s(1), 'A', A, 'check', chk);
end
