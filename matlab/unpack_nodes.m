function [Pv, Pw] = unpack_nodes(prob, x)
%UNPACK_NODES  Node positions of vessel and wire from the vector of unknowns.
%   [Pv, Pw] = UNPACK_NODES(prob, x) restores the held nodes: the vessel ends at
%   their State 1 positions and the first wire node at the sheath.
Pv = [prob.P(1,:); reshape(x(1:prob.nv), 2, [])'; prob.P(end,:)];
Pw = [prob.a;      reshape(x(prob.nv+1:end), 2, [])'];
end
