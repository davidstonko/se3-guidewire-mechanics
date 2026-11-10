function [f, g, H] = wire_bend3d(Q, c, cols, n)
%WIRE_BEND3D  Bending energy of a twist-free isotropic discrete rod in 3D,
%   with exact gradient and Hessian.
%
%   f = (c/2) sum_i |kb_i|^2 / l_i,
%
%   kb_i = 2 a x b / (|a||b| + a.b)   (curvature binormal, Bergou et al. 2008),
%   l_i  = (|a| + |b|)/2              (Voronoi length),
%
%   with a = q_i - q_{i-1} and b = q_{i+1} - q_i the edges at the interior
%   vertex i. Since |a x b|^2 = (|a||b| - a.b)(|a||b| + a.b), with the cosine
%   cs = a.b/(|a||b|) of the turning angle phi,
%       |kb|^2 = g(cs) = 4 (1 - cs)/(1 + cs) = 4 tan^2(phi/2),
%   so at small angles f reduces to (c/2) sum phi_i^2 / l_i, the planar energy.
%
%   Inputs
%     Q     n_nodes-by-3 node positions
%     c     bending stiffness (N mm^2)
%     cols  n_nodes-by-3 column of each coordinate in x (0 for a held node)
%     n     length of x
%   Outputs
%     f     energy (N mm);  g  n-by-1 gradient;  H  n-by-n sparse Hessian
%
%   Derivatives, per vertex, in the edge variables (a, b), with na = |a|,
%   nb = |b|:
%     grad_a cs = b/(na nb) - cs a/na^2,       grad_b cs = a/(na nb) - cs b/nb^2,
%     Hess_aa cs = -(b a' + a b')/(na^3 nb) - cs I/na^2 + 3 cs a a'/na^4,
%     Hess_ab cs = I/(na nb) - b b'/(na nb^3) - a a'/(na^3 nb) + cs a b'/(na^2 nb^2),
%     Hess_bb cs = (a <-> b in Hess_aa),
%     grad l = (a/na, b/nb)/2,  Hess l = blkdiag((I - a a'/na^2)/na, (I - b b'/nb^2)/nb)/2,
%     g' = -8/(1 + cs)^2,  g'' = 16/(1 + cs)^3,
%   and f_i = (c/2) g(cs)/l is differentiated by the chain rule. The vertex
%   derivatives are mapped to the three nodes by (a, b) = T (q_{i-1}, q_i, q_{i+1}).

a  = Q(2:end-1,:) - Q(1:end-2,:);
b  = Q(3:end,:)   - Q(2:end-1,:);
na = sqrt(sum(a.^2, 2));  nb = sqrt(sum(b.^2, 2));
cs = sum(a.*b, 2)./(na.*nb);
l  = 0.5*(na + nb);
gk = 4*(1 - cs)./(1 + cs);                 % |kb|^2
f  = 0.5*c*sum(gk./l);
if nargout < 2, return; end

nvx = size(a,1);
g1 = -8./(1 + cs).^2;                      % dg/dcs
ca = b./(na.*nb) - cs.*a./na.^2;           % grad_a cs
cb = a./(na.*nb) - cs.*b./nb.^2;           % grad_b cs
la = 0.5*a./na;  lb = 0.5*b./nb;           % grad l
fc = 0.5*c*g1./l;                          % df/dcs
fl = -0.5*c*gk./l.^2;                      % df/dl
gab = [fc.*ca + fl.*la, fc.*cb + fl.*lb];  % nvx-by-6

I3 = eye(3);  Z3 = zeros(3);
T = [-I3 I3 Z3; Z3 -I3 I3];                % (a, b) = T * (q_{i-1}, q_i, q_{i+1})
gq = gab*T;                                % nvx-by-9
C9 = [cols(1:end-2,:), cols(2:end-1,:), cols(3:end,:)];
keep = C9 > 0;
g = accumarray(C9(keep), gq(keep), [n 1]);
if nargout < 3, return; end

g2 = 16./(1 + cs).^3;                      % d2g/dcs2
Hq = zeros(9, 9, nvx);
for i = 1:nvx
    av = a(i,:)';  bv = b(i,:)';  A_ = na(i);  B_ = nb(i);  C_ = cs(i);  L_ = l(i);
    Haa = -(bv*av' + av*bv')/(A_^3*B_) - C_*I3/A_^2 + 3*C_*(av*av')/A_^4;
    Hbb = -(av*bv' + bv*av')/(B_^3*A_) - C_*I3/B_^2 + 3*C_*(bv*bv')/B_^4;
    Hab = I3/(A_*B_) - (bv*bv')/(A_*B_^3) - (av*av')/(A_^3*B_) + C_*(av*bv')/(A_^2*B_^2);
    Hcs = [Haa, Hab; Hab', Hbb];
    Hl  = 0.5*blkdiag((I3 - (av*av')/A_^2)/A_, (I3 - (bv*bv')/B_^2)/B_);
    dc = [ca(i,:), cb(i,:)]';  dl = [la(i,:), lb(i,:)]';
    fcc = 0.5*c*g2(i)/L_;  fcl = -0.5*c*g1(i)/L_^2;  fll = c*gk(i)/L_^3;
    H6 = fcc*(dc*dc') + fcl*(dc*dl' + dl*dc') + fll*(dl*dl') + fc(i)*Hcs + fl(i)*Hl;
    Hq(:,:,i) = T'*H6*T;
end
Ii = repmat(reshape(C9', 9, 1, nvx), 1, 9, 1);
Jj = repmat(reshape(C9', 1, 9, nvx), 9, 1, 1);
keep = Ii > 0 & Jj > 0;
H = sparse(Ii(keep), Jj(keep), Hq(keep), n, n);
end
