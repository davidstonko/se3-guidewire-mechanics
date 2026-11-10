function [Pv, Pw] = unpack_nodes3d(prob, x)
%UNPACK_NODES3D  Node positions of vessel and wire from the vector of unknowns.
%   [Pv, Pw] = UNPACK_NODES3D(prob, x) restores the held nodes: the vessel ends
%   (or the whole vessel if prob.rigid) at State 1 and the first wire node at
%   the sheath.
if prob.rigid
    Pv = prob.P;
else
    Pv = [prob.P(1,:); reshape(x(1:prob.nv), 3, [])'; prob.P(end,:)];
end
Pw = [prob.a; reshape(x(prob.nv+1:end), 3, [])'];
end
