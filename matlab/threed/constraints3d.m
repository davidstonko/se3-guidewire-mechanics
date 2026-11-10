function [c, ceq, Gc, Geq] = constraints3d(x, prob)
%CONSTRAINTS3D  Nonlinear constraints in the form fmincon expects.
%
%   c    (n_w-1)-by-1  lumen constraint at wire nodes 1..n_w-1   (<= 0)
%   ceq  (n_w+1)-by-1  edge lengths |y_j - y_{j-1}|^2 - h_w^2 (= 0), j = 1..n_w-1,
%                      followed by the sheath direction (y_1 - y_0).f1 = 0 and
%                      (y_1 - y_0).f2 = 0, which pin the first edge to t0
%   Gc, Geq            constraint gradients, one column per constraint

[Pv, Pw] = unpack_nodes3d(prob, x);
n = prob.n;  M = prob.M;
if nargout > 2
    [c, Gc] = contact_terms3d(Pw(2:end,:), Pv, prob.R, prob.rw, prob.wcol(2:end,:), prob.vcol, n);
else
    c = contact_terms3d(Pw(2:end,:), Pv, prob.R, prob.rw, prob.wcol(2:end,:), prob.vcol, n);
end
E = diff(Pw);                                     % wire edges
ceq = [sum(E.^2, 2) - prob.h0^2; E(1,:)*prob.f1'; E(1,:)*prob.f2'];
if nargout > 2
    % d/dy_j |y_j - y_{j-1}|^2 = 2 e_j,  d/dy_{j-1} = -2 e_j
    jj = (1:M-1)';
    cj = prob.wcol(2:end,:);  cm = prob.wcol(1:end-1,:);
    I = [cj(:); cm(:)];  J = repmat(jj, 6, 1);  V = [2*E(:); -2*E(:)];
    keep = I > 0;
    % sheath direction: linear in y_1 (y_0 is held)
    I = [I(keep); prob.wcol(2,:)'; prob.wcol(2,:)'];
    J = [J(keep); M; M; M; M+1; M+1; M+1];
    V = [V(keep); prob.f1'; prob.f2'];
    Geq = sparse(I, J, V, n, M+1);
end
end
