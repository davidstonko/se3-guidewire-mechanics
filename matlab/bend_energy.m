function [f, g, H] = bend_energy(Q, phi0, c, cols, n)
%BEND_ENERGY  Turning-angle bending energy of a planar polygon, with exact
%   gradient and Hessian.
%
%   f = sum_i  c/2 * (phi_i - phi0_i)^2 / l_i
%
%   where, at each interior vertex i, phi_i is the signed turning angle
%   between the incoming edge a = q_i - q_{i-1} and the outgoing edge
%   b = q_{i+1} - q_i, and l_i = (|a| + |b|)/2 is the Voronoi length.
%   By Lemma 7.1 of the paper this is the quadratic strain energy of
%   the discrete Bishop frame of the polygon.
%
%   Inputs
%     Q     n_nodes-by-2 node positions
%     phi0  (n_nodes-2)-by-1 reference turning angles (zeros for a straight wire)
%     c     bending stiffness (N mm^2)
%     cols  n_nodes-by-2 column of each coordinate in x (0 for a held node)
%     n     length of x
%   Outputs
%     f     energy (N mm)
%     g     n-by-1 gradient
%     H     n-by-n sparse Hessian (only if requested)
%
%   Derivatives. Locally phi = theta(b) - theta(a) with theta(v) = atan2(v_y, v_x),
%     grad theta(v) = (-v_y, v_x)/|v|^2,
%     Hess theta(v) = [2 v_x v_y, v_y^2 - v_x^2; v_y^2 - v_x^2, -2 v_x v_y]/|v|^4,
%   and grad |v| = v/|v|, Hess |v| = (I - v v'/|v|^2)/|v|. The energy of one
%   vertex is differentiated in the edge variables (a, b) and mapped to the
%   three nodes by a = q_i - q_{i-1}, b = q_{i+1} - q_i.

a  = Q(2:end-1,:) - Q(1:end-2,:);
b  = Q(3:end,:)   - Q(2:end-1,:);
phi = atan2(a(:,1).*b(:,2) - a(:,2).*b(:,1), sum(a.*b, 2));
na = sqrt(sum(a.^2, 2));  nb = sqrt(sum(b.^2, 2));
l  = 0.5*(na + nb);
r  = phi - phi0;
f  = 0.5*c*sum(r.^2./l);
if nargout < 2, return; end

nvx = size(a,1);
% First derivatives in the edge variables, [a_x a_y b_x b_y] per vertex.
dphi = [a(:,2)./na.^2, -a(:,1)./na.^2, -b(:,2)./nb.^2, b(:,1)./nb.^2];
dl   = 0.5*[a./na, b./nb];
fphi = c*r./l;                 % df/dphi
fl   = -0.5*c*r.^2./l.^2;      % df/dl
gab  = fphi.*dphi + fl.*dl;    % nvx-by-4

% Node map: [q_{i-1}; q_i; q_{i+1}] (6) -> [a; b] (4), (a, b) = T * q.
T = [-1 0 1 0 0 0; 0 -1 0 1 0 0; 0 0 -1 0 1 0; 0 0 0 -1 0 1];
gq = gab*T;                    % nvx-by-6, gradient per vertex in node coordinates
C6 = [cols(1:end-2,:), cols(2:end-1,:), cols(3:end,:)];   % nvx-by-6 columns
keep = C6 > 0;
g = accumarray(C6(keep), gq(keep), [n 1]);
if nargout < 3, return; end

% Second derivatives in (a, b), 4-by-4 per vertex.
Hab = zeros(4,4,nvx);
for i = 1:nvx
    ax = a(i,1); ay = a(i,2); bx = b(i,1); by = b(i,2);
    Hth_a = [2*ax*ay, ay^2-ax^2; ay^2-ax^2, -2*ax*ay]/na(i)^4;
    Hth_b = [2*bx*by, by^2-bx^2; by^2-bx^2, -2*bx*by]/nb(i)^4;
    Hphi = blkdiag(-Hth_a, Hth_b);
    ua = a(i,:)'/na(i);  ub = b(i,:)'/nb(i);
    Hl = 0.5*blkdiag((eye(2) - ua*ua')/na(i), (eye(2) - ub*ub')/nb(i));
    gp = dphi(i,:)';  gl = dl(i,:)';
    Hab(:,:,i) = c/l(i)*(gp*gp') + fphi(i)*Hphi ...
        - c*r(i)/l(i)^2*(gp*gl' + gl*gp') + c*r(i)^2/l(i)^3*(gl*gl') + fl(i)*Hl;
end
Hq = pagemtimes(pagemtimes(T', Hab), T);   % 6-by-6 per vertex
[I, J, V] = local_triplets(Hq, C6);
H = sparse(I, J, V, n, n);
end

function [I, J, V] = local_triplets(Hq, C6)
% Scatter 6-by-6 local Hessians into global (row, column, value) triplets.
nvx = size(C6,1);
Ii = repmat(reshape(C6', 6, 1, nvx), 1, 6, 1);
Jj = repmat(reshape(C6', 1, 6, nvx), 6, 1, 1);
keep = Ii > 0 & Jj > 0;
I = Ii(keep);  J = Jj(keep);  V = Hq(keep);
end
