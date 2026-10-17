%% =========================================================================
% RUN_ALL.M  -  Run the whole quadcopter pipeline in the right order
% =========================================================================
%   1. Nonlinear_Sym          symbolic model, A_sym, B_sym, f_sym   (skipped if cached)
%   2. trajectory_generation  flat-output plan: Xr, ur, tt
%   3. linearize              LTV model along the plan: A_num, B_num, Ck, idx_lin
%   4. new_compensator        E_u, H_global (parallel feedforward) + ASPR checks
%   5. run_sac                one SAC run on the nonlinear plant (+ plots)
%   6. sweep_all              model-mismatch sweep SAC / PID / LQR / LQI (slow, optional)
%
% Usage: open MATLAB in the repo root and run   run_all
%
% This is a SCRIPT on purpose: the pipeline scripts exchange data through the
% base workspace (assignin / exist checks), so they must run in it.
% Options are kept in appdata because Nonlinear_Sym.m and run_sac.m call
% clear / clearvars, which would wipe ordinary variables.
% =========================================================================

%% --- Options ---------------------------------------------------------------
clear; close all; clc;

ra.force_regen = false;              % true: delete the generated *_fun_gen.m files and rebuild
ra.run_single  = true;               % step 5: one SAC run with plots
ra.run_sweep   = true;              % step 6: mismatch sweep (minutes; needs steps 1-4)
ra.sweep_delta = [0 0.1 0.2 0.4];    % mismatch levels for the sweep
ra.save_figs   = true;               % save all open figures as PNG at the end
ra.fig_dir     = fullfile('figures', 'results');

%% --- Setup -------------------------------------------------------------------
ra_root = fileparts(mfilename('fullpath'));
if ~isempty(ra_root), cd(ra_root); end          % generated files are written to the cwd

if ~license('test', 'Symbolic_Toolbox')
    warning('run_all:toolbox', 'Symbolic Math Toolbox not found: step 1 will fail unless the generated files exist.');
end
if ~license('test', 'Control_Toolbox')
    error('run_all:toolbox', 'Control System Toolbox is required (tf, ss, lqr).');
end

% Remove a stale mismatch hook left by an interrupted sweep (would force run_sac into sweep mode)
if isappdata(0, 'mm_cfg'), rmappdata(0, 'mm_cfg'); end
if isappdata(0, 'mm_res'), rmappdata(0, 'mm_res'); end

ra.t_start = now;                    % wall-clock start (tic/toc would be reset by run_sac)
setappdata(0, 'run_all_cfg', ra);

%% --- 1. Symbolic model (skipped when the cached numeric functions are valid) --
ra_files = {'A_fun_gen', 'B_fun_gen', 'f_fun_gen'};
if ra.force_regen
    for ra_i = 1:numel(ra_files)
        ra_f = which([ra_files{ra_i} '.m']);
        if ~isempty(ra_f), delete(ra_f); end
    end
    rehash;
end

ra_cache_ok = true;
ra_src = which('Nonlinear_Sym.m');
if isempty(ra_src)
    error('run_all:nosym', 'Nonlinear_Sym.m not found. Run from the repo root.');
end
ra_tsrc = dir(ra_src);
for ra_i = 1:numel(ra_files)
    ra_f = which([ra_files{ra_i} '.m']);
    if isempty(ra_f)
        ra_cache_ok = false;
    else
        ra_tf = dir(ra_f);
        if ra_tf.datenum < ra_tsrc.datenum, ra_cache_ok = false; end
    end
end

if ra_cache_ok
    disp('[1/5] Generated functions are up to date: skipping Nonlinear_Sym.');
else
    disp('[1/5] Building the symbolic model (Nonlinear_Sym) ...');
    Nonlinear_Sym;                               % NOTE: this script calls clear
end
ra = getappdata(0, 'run_all_cfg');               % restore options (Nonlinear_Sym may have cleared them)

%% --- 2. Trajectory ---------------------------------------------------------------
disp('[2/5] Generating the reference trajectory ...');
trajectory_generation;

%% --- 3. Linearization along the plan ------------------------------------------------
% lz_check = true in linearize.m is required: it also generates f_fun_gen.m (nonlinear plant)
disp('[3/5] Linearizing along the trajectory ...');
linearize;

%% --- 4. Compensator + ASPR checks --------------------------------------------------
disp('[4/5] Designing the parallel feedforward compensator ...');
new_compensator;
if exist('overall_aspr', 'var') && ~overall_aspr
    warning('run_all:aspr', 'ASPR check failed at some frozen samples (see the summary above).');
end

%% --- 5. Single SAC run ------------------------------------------------------------------
if ra.run_single
    disp('[5/5] Running SAC on the nonlinear plant ...');
    run_sac;                                     % NOTE: this script calls clearvars -except <pipeline vars>
    ra = getappdata(0, 'run_all_cfg');
end

%% --- 6. Model-mismatch sweep (optional) ---------------------------------------------
if ra.run_sweep
    disp('[6] Running the model-mismatch sweep (this takes a while) ...');
    ra_out = sweep_all(struct('delta', ra.sweep_delta));
    ra = getappdata(0, 'run_all_cfg');
    assignin('base', 'sweep_out', ra_out);
end

%% --- Save figures ---------------------------------------------------------------------
if ra.save_figs
    if ~exist(ra.fig_dir, 'dir'), mkdir(ra.fig_dir); end
    ra_figs = findall(0, 'type', 'figure');
    for ra_i = 1:numel(ra_figs)
        ra_nm = get(ra_figs(ra_i), 'Name');
        if isempty(ra_nm), ra_nm = sprintf('figure_%d', get(ra_figs(ra_i), 'Number')); end
        ra_nm = regexprep(ra_nm, '[^A-Za-z0-9]+', '_');
        print(ra_figs(ra_i), fullfile(ra.fig_dir, [ra_nm '.png']), '-dpng', '-r150');
    end
    fprintf('Saved %d figure(s) to %s\n', numel(ra_figs), ra.fig_dir);
end

fprintf('\nrun_all finished in %.1f s.\n', (now - ra.t_start)*86400);
rmappdata(0, 'run_all_cfg');