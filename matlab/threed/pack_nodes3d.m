function x = pack_nodes3d(prob, Pv, Pw)
%PACK_NODES3D  Stack the free vessel and wire nodes into the vector of unknowns.
%   x = PACK_NODES3D(prob, Pv, Pw) with Pv N-by-3 and Pw n_w-by-3 node positions.
%   The held nodes (vessel ends, sheath node of the wire, and the whole vessel
%   when prob.rigid is true) are dropped.
if prob.rigid
    xv = zeros(0,1);
else
    xv = reshape(Pv(2:end-1,:)', [], 1);
end
x = [xv; reshape(Pw(2:end,:)', [], 1)];
end
