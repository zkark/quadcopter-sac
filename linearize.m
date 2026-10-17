%% === Linearize the quadrotor model along a reference trajectory (efficient) ===
% Run order:  Nonlinear_Sym.m  ->  trajectory_generation.m  ->  linearize.m
%
% Efficiency: the symbolic -> numeric conversion (matlabFunction) is done ONCE and
% saved as A_fun_gen.m, B_fun_gen.m (and f_fun_gen.m for the optional check).
% Later runs reuse these files, for ANY trajectory, N, Tfinal or parameter values,
% as long as the model (Nonlinear_Sym.m) has not changed. They are regenerated
% automatically when Nonlinear_Sym.m is newer than the files, or when forced.
% When the cache is valid, no symbolic variables are needed in the workspace.
%
% Needs:    Xr (12xN), ur (4xN), tt (1xN)           (trajectory_generation.m)
%           A_sym, B_sym, f_sym, X_sym, u_sym        (only when regenerating)
% Produces: A_num (12x12xN), B_num (12x4xN), Ck, Dk, N, Tfinal
%           (+ idx_lin, tt_lin: the trajectory samples that were linearized)
%
% NOTE: drag below MUST equal tg_D in trajectory_generation.m.

%% --- Options ---
lz_force_regen = false;   % true: rebuild the generated function files
lz_check       = true;    % true: verify f(Xr,ur) = d(Xr)/dt (adds one call per sample)
lz_stride      = 1;       % linearize every lz_stride-th sample (1 = all samples)

%% --- Check trajectory data ---
lz_need = {'Xr','ur','tt'};
for lz_i = 1:numel(lz_need)
    if ~exist(lz_need{lz_i},'var')
        error('linearize:missing', ...
            'Variable "%s" not found. Run trajectory_generation.m first.', lz_need{lz_i});
    end
end

%% --- Numeric parameters (must match model + trajectory script) ---
lz_m   = 1.0;
lz_g   = 9.81;
lz_Ixx = 0.01;
lz_Iyy = 0.01;
lz_Izz = 0.02;
lz_Ax  = 0.25;  lz_Ay = 0.25;  lz_Az = 0.40;

% Output / feedthrough matrices (outputs: x, y, z, psi)
Ck = zeros(4,12);  Ck(1,1) = 1;  Ck(2,3) = 1;  Ck(3,5) = 1;  Ck(4,11) = 1;
Dk = zeros(4,4);

%% --- Consistency check with the trajectory script ---
if exist('tg_D','var')
    if norm(diag(tg_D) - [lz_Ax; lz_Ay; lz_Az]) > 1e-12
        warning('linearize:drag', ...
            ['Drag in trajectory_generation.m (tg_D = diag([%g %g %g])) differs from the drag used here.\n' ...
             'Set tg_D = diag([%g %g %g]) in trajectory_generation.m and rerun it.'], ...
            tg_D(1,1), tg_D(2,2), tg_D(3,3), lz_Ax, lz_Ay, lz_Az);
    end
end

%% --- Cached numeric functions: decide whether to (re)generate ---
lz_files = {'A_fun_gen','B_fun_gen'};
if lz_check, lz_files{end+1} = 'f_fun_gen'; end

lz_regen = lz_force_regen;
for lz_i = 1:numel(lz_files)
    if isempty(which([lz_files{lz_i} '.m']))
        lz_regen = true;                       % file missing
    end
end
if ~lz_regen
    lz_src = which('Nonlinear_Sym.m');         % model source newer than cache?
    if ~isempty(lz_src)
        lz_tsrc = dir(lz_src);
        for lz_i = 1:numel(lz_files)
            lz_tf = dir(which([lz_files{lz_i} '.m']));
            if lz_tf.datenum < lz_tsrc.datenum
                lz_regen = true;
            end
        end
    end
end

if lz_regen
    disp('linearize: generating numeric functions from the symbolic model (one-time cost)...');
    lz_symneed = {'A_sym','B_sym','f_sym','X_sym','u_sym'};
    for lz_i = 1:numel(lz_symneed)
        if ~exist(lz_symneed{lz_i},'var')
            error('linearize:missingsym', ...
                'Variable "%s" not found. Run Nonlinear_Sym.m first (needed to generate the functions).', lz_symneed{lz_i});
        end
    end
    % Same names AND same assumptions as Nonlinear_Sym.m
    syms m g real
    syms Ixx Iyy Izz real
    syms Ax Ay Az real
    lz_vars = {X_sym, u_sym, m, g, Ixx, Iyy, Izz, Ax, Ay, Az};

    matlabFunction(A_sym, 'File', 'A_fun_gen', 'Vars', lz_vars);
    matlabFunction(B_sym, 'File', 'B_fun_gen', 'Vars', lz_vars);
    if lz_check
        matlabFunction(f_sym, 'File', 'f_fun_gen', 'Vars', lz_vars);
    end
    rehash;
    for lz_i = 1:numel(lz_files), clear(lz_files{lz_i}); end
else
    disp('linearize: using cached numeric functions (A_fun_gen, B_fun_gen).');
end

A_fun = @A_fun_gen;
B_fun = @B_fun_gen;
if lz_check, f_fun = @f_fun_gen; end

%% --- Sample selection and time grid ---
lz_Nfull = size(Xr,2);
Tfinal   = tt(end);                    % exported: test.m / 1.m read Tfinal and N
lz_dt    = tt(2) - tt(1);              % = Tfinal/(N-1)
idx_lin  = 1:lz_stride:lz_Nfull;
tt_lin   = tt(idx_lin);
N        = numel(idx_lin);             % number of linearization points

%% --- Evaluate A, B (and optionally f) along the trajectory ---
A_num = zeros(12,12,N);
B_num = zeros(12,4,N);
if lz_check, lz_res = zeros(12,N); end

for lz_k = 1:N
    lz_Xk = Xr(:,idx_lin(lz_k));
    lz_uk = ur(:,idx_lin(lz_k));

    A_num(:,:,lz_k) = A_fun(lz_Xk, lz_uk, lz_m, lz_g, lz_Ixx, lz_Iyy, lz_Izz, lz_Ax, lz_Ay, lz_Az);
    B_num(:,:,lz_k) = B_fun(lz_Xk, lz_uk, lz_m, lz_g, lz_Ixx, lz_Iyy, lz_Izz, lz_Ax, lz_Ay, lz_Az);
    if lz_check
        lz_res(:,lz_k) = f_fun(lz_Xk, lz_uk, lz_m, lz_g, lz_Ixx, lz_Iyy, lz_Izz, lz_Ax, lz_Ay, lz_Az);
    end
end

%% --- Optional feasibility check: f(Xr,ur) should equal d(Xr)/dt ---
if lz_check
    lz_Xdot  = gradient(Xr, lz_dt);                    % derivative along time (columns)
    lz_err   = lz_res - lz_Xdot(:,idx_lin);
    lz_inner = (idx_lin >= 4) & (idx_lin <= lz_Nfull-3);   % skip one-sided edge samples
    lz_maxerr = max(abs(lz_err(:,lz_inner)), [], 2);

    fprintf('Max |f(Xr,ur) - dXr/dt| per state (interior samples):\n');
    fprintf('  %s\n', mat2str(lz_maxerr.', 3));
    if max(lz_maxerr([2 4 6])) > 1e-2
        warning('linearize:infeasible', ...
            ['Reference (Xr,ur) does not satisfy the model (check drag / sign conventions).\n' ...
             'The linearization is about an infeasible trajectory.']);
    end
end

%% --- Export to workspace ---
assignin('base','A_num',A_num);
assignin('base','B_num',B_num);
assignin('base','Ck',Ck);
assignin('base','Dk',Dk);
assignin('base','N',N);
assignin('base','Tfinal',Tfinal);
assignin('base','idx_lin',idx_lin);
assignin('base','tt_lin',tt_lin);

fprintf('Linearized trajectory model: A_num (12x12x%d), B_num (12x4x%d)\n', N, N);