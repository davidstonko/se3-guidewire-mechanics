function [c, Gc, Hc, info] = contact_terms3d(Y, Pv, R, rw, ycols, vcols, n, lam)
%CONTACT_TERMS3D  Lumen constraint at the wire nodes in 3D, with exact derivatives.
%
%   c_j = dist(y_j, Gamma_h)^2 - delta_j^2 <= 0,   delta_j = R(sigma*_j) - r_w,
%
%   where Gamma_h is the vessel polygon Pv, the distance is to its nearest
%   segment [A, B] (searched over the whole polygon), and R is interpolated
%   linearly along that segment at the foot point A + t (B - A). The nearest
%   segment is chosen without differentiation and then held fixed. This is
%   CONTACT_TERMS of the planar code with 3-vectors in place of 2-vectors;
%   the formulas below do not depend on the dimension.
%
%   Inputs
%     Y      m-by-3 constrained wire nodes
%     Pv     N-by-3 vessel nodes;  R  N-by-1 lumen radii;  rw  wire radius
%     ycols  m-by-3 columns of the wire coordinates in x (0 = held)
%     vcols  N-by-3 columns of the vessel coordinates in x (0 = held)
%     n      length of x
%     lam    (optional) m-by-1 multipliers; then Hc = sum_j lam_j * Hess c_j
%   Outputs
%     c      m-by-1 constraint values (<= 0 inside the lumen)
%     Gc     n-by-m sparse constraint gradients (columns)
%     Hc     n-by-n sparse weighted constraint Hessian
%     info   struct: segment index seg, foot parameter t, distance d, delta
%
%   Derivatives, in w = y - A and e = B - A, with p = w.e, q = e.e, t = p/q:
%     interior foot point (0 < t < 1):
%       d^2 = w.w - p^2/q,  grad_w d^2 = 2 r,  grad_e d^2 = -2 t r,  r = w - t e,
%       Hess_ww = 2 I - 2 e e'/q,
%       Hess_we = -2 (e w'/q + t I - 2 t e e'/q),
%       Hess_ee = -2 (w - 2 t e)(w - 2 t e)'/q + 2 t^2 I,
%       grad_w t = e/q,  grad_e t = (w - 2 t e)/q,
%       Hess_we t = I/q - 2 e e'/q^2,
%       Hess_ee t = -2 [e (w - 2 t e)' + (w - 2 t e) e']/q^2 - 2 t I/q,
%     foot point at a vertex (t clamped to 0 or 1): d^2 = |y - vertex|^2 with
%       the radius of that vertex.

m = size(Y,1);
A = Pv(1:end-1,:);  Bn = Pv(2:end,:);  E = Bn - A;  qE = sum(E.^2, 2);
% Nearest segment of every point (no derivatives).
seg = zeros(m,1);  tt = zeros(m,1);  d2 = zeros(m,1);
for j = 1:m
    W  = Y(j,:) - A;
    t  = min(max(sum(W.*E, 2)./qE, 0), 1);
    D2 = sum((W - t.*E).^2, 2);
    [d2(j), s] = min(D2);
    seg(j) = s;  tt(j) = t(s);
end
Rt = (1 - tt).*R(seg) + tt.*R(seg+1);
del = Rt - rw;
c = d2 - del.^2;
info = struct('seg', seg, 't', tt, 'd', sqrt(d2), 'delta', del);
if nargout < 2, return; end

wantH = nargout > 2 && nargin > 7;
I3 = eye(3);  Z3 = zeros(3);
Tm = [I3 -I3 Z3; Z3 -I3 I3];           % [w; e] = Tm [y; A; B]
GI = zeros(9*m,1); GJ = GI; GV = GI;  HI = zeros(81*m,1); HJ = HI; HV = HI;
for j = 1:m
    s = seg(j);  w = (Y(j,:) - A(s,:))';  e = E(s,:)';  q = qE(s);  t = tt(j);
    dR = R(s+1) - R(s);
    if t > 0 && t < 1
        rr = w - t*e;  u = w - 2*t*e;
        gd = [2*rr; -2*t*rr];
        Hd = [2*I3 - 2*(e*e')/q, -2*(e*w'/q + t*I3 - 2*t*(e*e')/q)];
        Hd = [Hd; Hd(:,4:6)', -2*(u*u')/q + 2*t^2*I3]; %#ok<AGROW>
        gt = [e/q; u/q];
        Hte = I3/q - 2*(e*e')/q^2;
        Htee = -2*(e*u' + u*e')/q^2 - 2*t*I3/q;
        Ht = [Z3, Hte; Hte', Htee];
        gw = gd - 2*del(j)*dR*gt;
        Hw = Hd - 2*dR^2*(gt*gt') - 2*del(j)*dR*Ht;
    elseif t == 0                       % foot point at A
        gw = [2*w; 0; 0; 0];
        Hw = blkdiag(2*I3, Z3);
    else                                % foot point at B
        gw = [2*(w - e); -2*(w - e)];
        Hw = [2*I3, -2*I3; -2*I3, 2*I3];
    end
    cl = [ycols(j,:), vcols(s,:), vcols(s+1,:)];
    gz = Tm'*gw;
    k9 = (j-1)*9 + (1:9);
    GI(k9) = cl;  GJ(k9) = j;  GV(k9) = gz;
    if wantH
        Hz = lam(j)*(Tm'*Hw*Tm);
        [Ii, Jj] = ndgrid(cl, cl);
        k81 = (j-1)*81 + (1:81);
        HI(k81) = Ii(:);  HJ(k81) = Jj(:);  HV(k81) = Hz(:);
    end
end
keep = GI > 0;
Gc = sparse(GI(keep), GJ(keep), GV(keep), n, m);
Hc = [];
if wantH
    keep = HI > 0 & HJ > 0;
    Hc = sparse(HI(keep), HJ(keep), HV(keep), n, n);
end
end
