function x = pack_nodes(prob, Pv, Pw)
%PACK_NODES  Stack the free vessel and wire nodes into the vector of unknowns.
%   x = PACK_NODES(prob, Pv, Pw) with Pv N-by-2 and Pw n_w-by-2 node positions.
%   The held nodes (vessel ends, sheath node of the wire) are dropped.
x = [reshape(Pv(2:end-1,:)', [], 1); reshape(Pw(2:end,:)', [], 1)];
end
