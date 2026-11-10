function [f, g, H, parts] = objective3d(x, prob)
%OBJECTIVE3D  Total elastic energy E = E_wire + E_vessel + E_foundation,
%   with exact gradient and Hessian.
%
%   E_wire       = B/2   * sum |kb_i|^2 / l_i                    (WIRE_BEND3D)
%   E_vessel     = B_v/2 * sum |u_{i+1} - 2 u_i + u_{i-1}|^2 / h^3
%   E_foundation = 1/2   * sum k_i |u_i|^2,     u_i = x_i - P_i
%
%   The vessel terms are quadratic, 1/2 u'(Kb + Kk) u, with the matrices built
%   once in BUILD_PROBLEM3D. With prob.rigid the vessel terms are zero.
%   parts = [E_wire, E_vessel, E_foundation].

[~, Pw] = unpack_nodes3d(prob, x);
n = prob.n;  nv = prob.nv;
needH = nargout > 2;
if needH
    [fw, g, Hw] = wire_bend3d(Pw, prob.B, prob.wcol, n);
else
    [fw, g] = wire_bend3d(Pw, prob.B, prob.wcol, n);
end
fb = 0;  ff = 0;
if ~prob.rigid
    u  = x(1:nv) - prob.xP;
    gb = prob.Kb*u;  gk = prob.Kk*u;
    fb = 0.5*(u'*gb);  ff = 0.5*(u'*gk);
    g(1:nv) = g(1:nv) + gb + gk;
end
f = fw + fb + ff;
parts = [fw, fb, ff];
if needH
    H = Hw;
    if ~prob.rigid, H = H + prob.Hv; end
end
end
