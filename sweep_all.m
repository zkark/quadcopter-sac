function out = sweep_all(opt)
% SWEEP_ALL  Model-mismatch sweep for SAC (run_sac.m) and the baselines PID / LQR / LQI.
%
%   out = sweep_all;                                  % default sweep
%   out = sweep_all(struct('delta',[0 0.2]));         % quick test of the run time
%
% For every case (wrong mass / inertia / drag) and every level delta the PLANT
% parameters are nominal*(1+delta) while every controller keeps its NOMINAL design.
% All four controllers see the same plan, the same initial offset and the same
% actuator limits, and are scored with the same metrics (calc_metrics.m).
%
% This is a FUNCTION on purpose: run_sac.m starts with clearvars, which would wipe
% the variables of a script. The function keeps its own workspace, reads the
% pipeline variables from the base workspace and calls run_sac.m there.
% The mismatched parameters reach run_sac.m through setappdata(0,'mm_cfg').
%
% Needs (base workspace): A_num B_num Ck Xr ur idx_lin H_global E_u psi_k (+ tt or tt_lin)
% Needs (path):           run_sac.m (patched), baseline_sim.m, calc_metrics.m, f_fun_gen.m
%
% opt fields (all optional): delta, baseline_opt (struct passed to baseline_sim)
% Result: out.rms / out.fin / out.crash (case x delta x controller), out.Kmax (case x delta)
%         also saved to sweep_results.mat and to the base workspace as sweep_out.

if nargin < 1 || isempty(opt), opt = struct(); end
delta = getopt(opt, 'delta', [0 0.1 0.2 0.4]);
cases = {'mass', 1; 'inertia', 3:5; 'drag', 6:8};      % name, index in [m g Ixx Iyy Izz Ax Ay Az]
ctrls = {'SAC','PID','LQR','LQI'};
p0    = [1.0 9.81 0.01 0.01 0.02 0.25 0.25 0.40];      % nominal, must match run_sac.m

%% --- Checks ---------------------------------------------------------------
need = {'A_num','B_num','Ck','Xr','ur','idx_lin','H_global','E_u','psi_k'};
for i = 1:numel(need)
    if ~hasb(need{i})
        error('sweep_all:missing', ...
            'Variable "%s" not found in the base workspace. Run the pipeline first.', need{i});
    end
end
if isempty(which('f_fun_gen'))
    error('sweep_all:nof', 'f_fun_gen.m not found. Run linearize.m with lz_check = true.');
end
if isempty(which('run_sac'))
    error('sweep_all:norun', 'run_sac.m not found. Rename run.m to run_sac.m and apply the patch.');
end

%% --- Pack the data for the baselines ----------------------------------------
W.A_num = getb('A_num');  W.B_num = getb('B_num');  W.Ck = getb('Ck');
W.Xr = getb('Xr');        W.ur = getb('ur');        W.idx_lin = getb('idx_lin');
if hasb('tt_lin')
    tv = getb('tt_lin');
elseif hasb('tt')
    tv = getb('tt');  tv = tv(W.idx_lin);
else
    error('sweep_all:time', 'Time vector (tt or tt_lin) not found.');
end
W.t_vec = tv(:).';
W.x_p0 = zeros(12,1);
W.x_p0([1 3 5]) = [0.10; -0.10; 0.05];                 % same offset as run_sac.m
W.p_nom = p0;
tau_lim = 1.0;  tau_psi_lim = 0.5;                     % same limits as run_sac.m
W.u_min = [0;                  -tau_lim; -tau_lim; -tau_psi_lim];
W.u_max = [2*p0(1)*p0(2);       tau_lim;  tau_lim;  tau_psi_lim];

C.W = W;
C.opt = getopt(opt, 'baseline_opt', struct());
C.crash_pos = 10;   C.crash_ang = 1.2;                 % crash envelope, same for all controllers
C.opt.crash_pos = C.crash_pos;  C.opt.crash_ang = C.crash_ang;

cleanupObj = onCleanup(@cleanup_mm);                   %#ok<NASGU>  always remove the mismatch hook

%% --- Sweep ----------------------------------------------------------------
nC = size(cases,1);  nD = numel(delta);  nK = numel(ctrls);
rms   = nan(nC, nD, nK);
fin   = nan(nC, nD, nK);
crash = false(nC, nD, nK);
Kmax  = nan(nC, nD);
nom   = [];
tstart = tic;

for ic = 1:nC
    for id = 1:nD
        p = p0;
        g = cases{ic,2};
        p(g) = p0(g)*(1 + delta(id));
        if delta(id) == 0 && ~isempty(nom)
            R = nom;                                   % nominal plant: reuse
        else
            fprintf('[%6.1f s] %-8s delta = %3.0f %% ...\n', toc(tstart), cases{ic,1}, 100*delta(id));
            R = eval_all(p, C);
            C.opt.G = R.G;                             % keep the LQ gains for the next runs
            if delta(id) == 0, nom = R; end
        end
        rms(ic,id,:)   = R.rms;
        fin(ic,id,:)   = R.fin;
        crash(ic,id,:) = R.crash;
        Kmax(ic,id)    = R.Kmax;
    end
end
fprintf('Sweep finished in %.1f s.\n', toc(tstart));

%% --- Tables ---------------------------------------------------------------
for ic = 1:nC
    fprintf('\nWrong %s: rms position error [m]   (CRASH = left the safe envelope)\n', cases{ic,1});
    fprintf('%8s', 'delta');  fprintf('%9s', ctrls{:});  fprintf('%14s\n', 'SAC max|Ki|');
    for id = 1:nD
        fprintf('%7.0f%%', 100*delta(id));
        for ik = 1:nK
            if crash(ic,id,ik), fprintf('%9s', 'CRASH');
            else, fprintf('%9.4f', rms(ic,id,ik)); end
        end
        fprintf('%14.3g\n', Kmax(ic,id));
    end
end

%% --- Degradation curves -----------------------------------------------------
figure('Name', 'Degradation under model mismatch');
for ic = 1:nC
    for row = 1:2
        subplot(2, nC, (row-1)*nC + ic);  hold on;
        for ik = 1:nK
            if row == 1, y = rms(ic,:,ik); else, y = fin(ic,:,ik); end
            lw = 1.3;  if ik == 1, lw = 2.8; end       % SAC drawn thicker
            plot(100*delta, y(:).', '-o', 'LineWidth', lw);
        end
        grid on;  xlabel('mismatch (%)');
        if row == 1
            title(['wrong ' cases{ic,1}]);  ylabel('rms position error (m)');
        else
            ylabel('final position error (m)');
        end
        if ic == 1 && row == 1, legend(ctrls, 'Location', 'northwest'); end
    end
end

%% --- Sanity check that the mismatch really reached the SAC plant -------------
x = rms(1,:,1);  x = x(isfinite(x));
if numel(x) > 1 && max(x) - min(x) < 1e-12
    warning('sweep_all:nochange', ...
        ['The SAC result does not change with the mass mismatch. Patch B (plant parameters from mm_cfg) ' ...
         'is probably missing in run_sac.m, or use_nl is false.']);
end

%% --- Output -----------------------------------------------------------------
out.ctrls = ctrls;  out.cases = cases(:,1);  out.delta = delta;
out.rms = rms;  out.fin = fin;  out.crash = crash;  out.Kmax = Kmax;
save('sweep_results.mat', 'out');
assignin('base', 'sweep_out', out);
end

%% =========================================================================
function R = eval_all(p, C)
% All four controllers on the plant with parameters p.
R.rms = nan(1,4);  R.fin = nan(1,4);  R.crash = false(1,4);  R.Kmax = nan;

% --- SAC: run_sac.m reads the plant parameters from the app data hook ---
setappdata(0, 'mm_cfg', struct('p', p));
evalc('evalin(''base'',''run_sac'');');                 % evalc hides the printouts
res = getappdata(0, 'mm_res');
if isempty(res)
    error('sweep_all:nores', 'run_sac.m did not return a result. Patch D (setappdata mm_res) is missing.');
end
ii   = 1:res.last;
Xr_s = C.W.Xr(:, C.W.idx_lin(ii));
ur_s = C.W.ur(:, C.W.idx_lin(ii));
E    = res.X_act - Xr_s;
ep   = abs(E([1 3 5],:));   ea = abs(E([7 9],:));
R.crash(1) = logical(res.aborted) || any(ep(:) > C.crash_pos) || any(ea(:) > C.crash_ang);
if ~R.crash(1)
    m = calc_metrics(res.X_act, res.U_app, Xr_s, ur_s, false);
    R.rms(1) = m.rms_pos;
    R.fin(1) = m.final_pos;
end
R.Kmax = max(res.K_norm);

% --- Baselines ---
o = C.opt;
names = {'PID','LQR','LQI'};
for k = 1:3
    r = baseline_sim(names{k}, p, C.W, o);
    o.G = r.G;
    R.crash(k+1) = r.diverged;
    if ~r.diverged
        R.rms(k+1) = r.metrics.rms_pos;
        R.fin(k+1) = r.metrics.final_pos;
    end
end
R.G = o.G;
end

%% =========================================================================
function v = getb(name)
v = evalin('base', name);
end

function tf = hasb(name)
tf = logical(evalin('base', ['exist(''' name ''',''var'')']));
end

function v = getopt(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

function cleanup_mm()
if isappdata(0, 'mm_cfg'), rmappdata(0, 'mm_cfg'); end
if isappdata(0, 'mm_res'), rmappdata(0, 'mm_res'); end
end