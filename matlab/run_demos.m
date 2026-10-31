%RUN_DEMOS  Computations reported in Section 4 of the paper.
%
%   Planar test vessel (EVAR_GEOMETRY), node spacing about 1.6 mm, 0.035 inch
%   Lunderquist Extra Stiff wire, B = 158.4 GPa * pi d^4/64 (Harrison et al.
%   2011), vessel bending stiffness 200 N mm^2, foundation 0.0158 N/mm^2.
%
%   1. Base case: tip 120 mm into the aorta (Section 4 and Figure 5).
%      The same continuation gives the wire in a nearly rigid vessel. Both
%      contact patterns are printed and the rigid one is saved to rigid.mat.
%   2. Wire length: tip 90, 60, 30, 0 mm and 150, 180 mm into the aorta,
%      each started from the previous equilibrium cut at the new length or
%      continued straight beyond the old tip (Section 4 and Figure 6).
%
%   Results are written to results/*.mat and results/summary.json, and the
%   command window output to results/run_demos.log. AUDIT_SOLUTION then
%   checks each saved equilibrium independently of the solver.

here = fileparts(mfilename('fullpath'));
cd(here);  if ~exist('results', 'dir'), mkdir results; end
if exist('results/run_demos.log', 'file'), delete('results/run_demos.log'); end
diary('results/run_demos.log');
fprintf('%s   MATLAB %s\n', char(datetime('now')), version);
h = 1.6;  kbar = 0.0158;
B = 158.4e3*pi*0.889^4/64;               % N mm^2
G = evar_geometry(h);

fprintf('Derivative check\n');  check_derivatives();

fprintf('Base case, tip 120 mm into the aorta\n');
[x, out, prob, steps, rigid] = continuation_solve(G, G.S_il + 120, B, kbar);
base = postprocess(prob, x, out);
rig  = postprocess(rigid.prob, rigid.x, rigid.out);
save('results/base.mat', '-struct', 'base');
summary.base = rmfield(base, {'Pv','Pw','P','R','U','s_v','s_w','fn','kap','cum_kkt','cum_line','He','N_edge','halfBk2','contact_nodes'});
summary.base.steps = steps;
summary.rigid = struct('F_kkt', rig.F_kkt, 'E', rig.E, 'flag', rig.flag, 'arcs', rig.arcs);
save('results/rigid.mat', '-struct', 'rig');
fprintf('  rigid vessel: wall reaction %.3f N, contacts:\n', rig.F_kkt);
for a = 1:numel(rig.arcs)
    fprintf('    s = %6.1f to %6.1f mm   %-6s %.3f N\n', rig.arcs(a).s0, rig.arcs(a).s1, rig.arcs(a).side, rig.arcs(a).force);
end
fprintf('  deformable vessel: wall reaction %.3f N, contacts:\n', base.F_kkt);
for a = 1:numel(base.arcs)
    fprintf('    s = %6.1f to %6.1f mm   %-6s %.3f N\n', base.arcs(a).s0, base.arcs(a).s1, base.arcs(a).side, base.arcs(a).force);
end

tips = {[90 60 30 0], [150 180]};
summary.length = struct('tip', {}, 'flag', {}, 'E', {}, 'peak', {}, 'F_kkt', {}, 'tip_force', {}, 's_last_contact', {}, 'time', {});
for q = 1:numel(tips)
    xx = x;  Pw_old = base.Pw;
    for tp = tips{q}
        Mn = round((G.S_il + tp)/G.hv) + 1;
        Pw = retail(Pw_old, Mn);
        pr = build_problem(G, (Mn - 1)*G.hv, B, kbar*G.hv, Pw);
        [Pv, ~] = unpack_nodes(prob, xx);
        x0 = pack_nodes(pr, Pv, Pw);
        [xt, oo] = solve_equilibrium(pr, x0, true, 600);
        if oo.flag ~= 1                     % report only solves that met the optimality tolerance
            fprintf('  tip %4d mm   flag %d (no converged solution)   %.1f s\n', tp, oo.flag, oo.time);
            summary.length(end+1) = struct('tip', tp, 'flag', oo.flag, 'E', oo.E, 'peak', NaN, ...
                'F_kkt', NaN, 'tip_force', NaN, 's_last_contact', NaN, 'time', oo.time);
            continue
        end
        xx = xt;
        S = postprocess(pr, xx, oo);
        fprintf('  tip %4d mm   flag %d   E = %.4f   peak %.4f   F = %.3f   tip force %.3f   %.1f s\n', ...
            tp, S.flag, S.E, S.peak, S.F_kkt, S.tip_force, S.time);
        summary.length(end+1) = struct('tip', tp, 'flag', S.flag, 'E', S.E, 'peak', S.peak, ...
            'F_kkt', S.F_kkt, 'tip_force', S.tip_force, 's_last_contact', S.s_last_contact, 'time', S.time);
        save(sprintf('results/tip%d.mat', tp), '-struct', 'S');
        Pw_old = S.Pw;  prob = pr;
    end
end
fid = fopen('results/summary.json', 'w');  fprintf(fid, '%s', jsonencode(summary, 'PrettyPrint', true));  fclose(fid);
fprintf('\nAudit of the saved equilibria\n');
for f = ["base", "tip90", "tip60", "tip150", "tip180"]
    audit_solution(sprintf('results/%s.mat', f));
end
fprintf('done\n');
diary('off');

function Pw = retail(Pw, Mn)
% Cut the wire to Mn nodes, or continue it straight beyond the tip.
M = size(Pw, 1);
if Mn <= M
    Pw = Pw(1:Mn, :);
else
    e = Pw(end,:) - Pw(end-1,:);
    Pw = [Pw; Pw(end,:) + (1:Mn-M)'.*e];
end
end
