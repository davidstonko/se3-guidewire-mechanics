function [c, ceq, Gc, Geq] = constraints(x, prob)
%CONSTRAINTS  Nonlinear constraints in the form fmincon expects.
%
%   c   (n_w-1)-by-1  lumen constraint at wire nodes 1..n_w-1   (<= 0)
%   ceq  n_w-by-1     edge lengths  |y_j - y_{j-1}|^2 - h_w^2   (= 0), j = 1..n_w-1,
%                     followed by the sheath direction (y_1 - y_0) x t_0 = 0
%   Gc, Geq           constraint gradients, one column per constraint

[Pv, Pw] = unpack_nodes(prob, x);
n = prob.n;  M = prob.M;
if nargout > 2
    [c, Gc] = contact_terms(Pw(2:end,:), Pv, prob.R, prob.rw, prob.wcol(2:end,:), prob.vcol, n);
else
    c = contact_terms(Pw(2:end,:), Pv, prob.R, prob.rw, prob.wcol(2:end,:), prob.vcol, n);
end
E = diff(Pw);                                     % wire edges
cs = prob.t0;
ceq = [sum(E.^2, 2) - prob.h0^2; E(1,1)*cs(2) - E(1,2)*cs(1)];
if nargout > 2
    % d/dy_j |y_j - y_{j-1}|^2 = 2 e_j,  d/dy_{j-1} = -2 e_j
    jj = (1:M-1)';
    cj = prob.wcol(2:end,:);  cm = prob.wcol(1:end-1,:);
    I = [cj(:); cm(:)];  J = [jj; jj; jj; jj];  V = [2*E(:); -2*E(:)];
    keep = I > 0;
    % sheath direction: linear in y_1
    I = [I(keep); prob.wcol(2,1); prob.wcol(2,2)];
    J = [J(keep); M; M];
    V = [V(keep); cs(2); -cs(1)];
    Geq = sparse(I, J, V, n, M);
end
end
