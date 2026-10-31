function G = evar_geometry(h)
%EVAR_GEOMETRY  Planar test vessel for the EVAR problem (State 1).
%
%   G = EVAR_GEOMETRY(h) returns the State 1 centerline sampled at uniform
%   arc length with spacing close to h (mm), and the lumen radius at every node.
%
%   The centerline is one period of the sine (x, 25 sin(2 pi x/150)) mm,
%   0 <= x <= 150 (the iliac segment), continued straight along its end
%   tangent for 200 mm (the aorta). The lumen radius is 6 mm in the iliac
%   segment and rises smoothly to 10 mm over the first 20 mm of the aorta.
%
%   Fields of G:
%     P     N-by-2 node positions of the State 1 centerline (mm)
%     R     N-by-1 lumen radius at each node (mm)
%     S_il  arc length of the iliac segment (mm)
%     hv    actual node spacing (mm)

L_ao = 200; R_il = 6; R_ao = 10; L_tr = 20;

% Densely sampled iliac sine, used only to interpolate at uniform arc length.
x  = linspace(0, 150, 200001)';
C  = [x, 25*sin(2*pi*x/150)];
sC = [0; cumsum(sqrt(sum(diff(C).^2, 2)))];
S_il = sC(end);
T  = C(end,:) - C(end-1,:);  T = T/norm(T);      % end tangent of the iliac segment

L = S_il + L_ao;
N = round(L/h) + 1;
s = linspace(0, L, N)';
P = zeros(N, 2);
il = s <= S_il;
P(il,1) = interp1(sC, C(:,1), s(il));
P(il,2) = interp1(sC, C(:,2), s(il));
P(~il,:) = C(end,:) + (s(~il) - S_il).*T;          % straight aorta

z = min(max((s - S_il)/L_tr, 0), 1);
R = R_il + (R_ao - R_il)*z.^2.*(3 - 2*z);          % C^1 smoothstep taper

G = struct('P', P, 'R', R, 'S_il', S_il, 'hv', s(2) - s(1));
end
