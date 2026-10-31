function [f, g, H, parts] = objective(x, prob)
%OBJECTIVE  Total elastic energy E = E_wire + E_vessel + E_foundation,
%   with exact gradient and Hessian.
%
%   E_wire       = B/2   * sum phi_j^2 / l_j            (straight, untorqued wire)
%   E_vessel     = B_v/2 * sum (phi_i - phi_i^ref)^2 / l_i
%   E_foundation = 1/2   * sum k_i |x_i - P_i|^2        (Winkler foundation)
%
%   parts = [E_wire, E_vessel, E_foundation].

[Pv, Pw] = unpack_nodes(prob, x);
n = prob.n;
needH = nargout > 2;
if needH
    [fw, gw, Hw] = bend_energy(Pw, zeros(prob.M-2,1), prob.B,  prob.wcol, n);
    [fv, gv, Hv] = bend_energy(Pv, prob.phi_ref,      prob.Bv, prob.vcol, n);
else
    [fw, gw] = bend_energy(Pw, zeros(prob.M-2,1), prob.B,  prob.wcol, n);
    [fv, gv] = bend_energy(Pv, prob.phi_ref,      prob.Bv, prob.vcol, n);
end
U  = Pv - prob.P;                                 % vessel displacement
ff = 0.5*sum(prob.k.*sum(U.^2, 2));
cols = prob.vcol(2:end-1,:);  kk = prob.k(2:end-1);
gf = zeros(n,1);  gf(cols(:)) = [kk; kk].*reshape(U(2:end-1,:), [], 1);
f = fw + fv + ff;
g = gw + gv + gf;
parts = [fw, fv, ff];
if needH
    Hf = sparse(cols(:), cols(:), [kk; kk], n, n);
    H = Hw + Hv + Hf;
end
end
