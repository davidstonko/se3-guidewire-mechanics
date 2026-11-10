function mesh_strict(h)
%MESH_STRICT  Base case at node spacing h with CONTINUATION_SOLVE_STRICT.
%   Prints the continuation log and the summary of POSTPROCESS.
B = 158.4e3*pi*0.889^4/64;  kbar = 0.0158;
G = evar_geometry(h);
fprintf('%s   h = %.2f mm (actual %.4f), %d vessel nodes\n', char(datetime('now')), h, G.hv, size(G.P,1));
[x, out, prob, ~, rigid] = continuation_solve_strict(G, G.S_il + 120, B, kbar);
Sr = postprocess(rigid.prob, rigid.x, rigid.out);
S  = postprocess(prob, x, out);
fprintf('  rigid: wall reaction %.3f N\n', Sr.F_kkt);
fprintf('  deformable: flag %d, E = %.4f, wall reaction %.3f N (line load %.3f N), peak %.4f mm at %.1f, |H| < %.4f N\n', ...
    S.flag, S.E, S.F_kkt, S.F_line, S.peak, S.peak_at_s, S.H_absmax_after_first_edge);
for a = 1:numel(S.arcs)
    fprintf('    s = %6.1f to %6.1f mm   %-6s %.3f N\n', S.arcs(a).s0, S.arcs(a).s1, S.arcs(a).side, S.arcs(a).force);
end
save(fullfile(fileparts(mfilename('fullpath')), 'results', sprintf('mesh_strict_h%s.mat', strrep(sprintf('%g', h), '.', 'p'))), 'x', 'out', 'prob', 'S', 'Sr');
fprintf('done %s\n', char(datetime('now')));
end
