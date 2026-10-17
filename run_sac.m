%% =========================================================================
% RUN_SAC.M  -  Closed-loop Simple Adaptive Control (SAC) simulation
% =========================================================================
% Run order: Nonlinear_Sym -> trajectory_generation -> linearize
%            -> new_compensator -> run_sac
%
% No function definitions in this file -> works in any MATLAB version.
%
% PLANT MODES (use_nl):
%   false : the plant is the linearized LTV model, frozen over each segment k:
%               d(dx)/dt = A_k*dx + B_k*du,    dx = X - Xr,  du = u - ur
%   true  : the plant is the REAL nonlinear model f(X,u) (f_fun_gen.m, made by
%           linearize.m with lz_check = true), driven by u = ur + du.
%           The controller only sees the deviation  X - Xr  (rotated outputs),
%           exactly as in the linear case. This is the validation that
%           matters: it shows whether the linearized design really works.
%
% Outputs are yaw-rotated: y = [x_b; y_b; z; psi] = rotC(Ck,psi_ref(k))*dx.
% Plant input du = E_u*u_tilde; the parallel feedforward H is driven by u_tilde.
% Commands are deviations from the plan (yaw frame); r = 0 means "follow Xr".
% Integration: RK4 with sub-steps on the whole loop.
% =========================================================================

clearvars -except A_num B_num Ck Dk N Tfinal H_global H_ss E_u psi_k ...
                  Xr ur idx_lin tt tt_lin;
mm = getappdata(0,'mm_cfg');     % set by sweep_all.m, empty when run by hand
%% --- Workspace checks ----------------------------------------------------
need = {'A_num','B_num','Ck','N','H_global','E_u','psi_k','Xr','ur'};
for i = 1:numel(need)
    if ~exist(need{i},'var')
        error('run_sac:missing', 'Variable "%s" not found. Run the earlier scripts first (new_compensator.m defines E_u, H_global, psi_k).', need{i});
    end
end
if ~exist('idx_lin','var'), idx_lin = 1:N; end
if exist('tt_lin','var')
    t_vec = tt_lin(:).';
elseif exist('tt','var')
    t_vec = tt(idx_lin);
else
    t_vec = linspace(0, Tfinal, N);
end
if ~exist('H_ss','var'), H_ss = ss(H_global); end

%% --- Options -------------------------------------------------------------
use_nl     = true;                   % true: nonlinear plant, false: linearized LTV plant
n_sub      = 20;                     % RK4 sub-steps per sample interval
gain_scale = 10000;                      % global multiplier on Gamma_p, Gamma_i
sigma_leak = 0.01;                    % leakage on K_i
blow_up    = 1e4;                    % abort threshold
match_ic   = true;                    % start the reference model at the plant's initial output
                                     %   (removes the initial error e_y(0) -> no initial kick)
q_max      = Inf;                    % cap on the proportional gain scalar rx'*Gamma_p*rx
                                     %   (the P-term is ey*q, i.e. quadratic in the error)

% Initial deviation from the plan, dx0 = X(0) - Xr(0)  (12 states)
x_p0 = zeros(12,1);
x_p0([1 3 5]) = [0.10; -0.10; 0.05];   % m offsets in x, y, z

% Command = world-frame offset from the plan [dx; dy; dz; dpsi], cosine ramp.
% Zeros -> the loop must follow Xr itself (best for comparing against Xr, ur).
cmd_amp  = [0; 0; 0; 0];             % e.g. [0.3; -0.3; 0.2; 0.15]
t_on     = 1.0;
t_ramp   = 1.0;

% Actuator limits on the physical input u = [T; tau_phi; tau_theta; tau_psi]
% Set these from your real actuators. Do NOT derive them from ur: the planned
% torques are only ~1e-3 N*m, far too small to allow recovery from an error.
tau_lim = 1.0;  tau_psi_lim = 0.5;                 % N*m
u_min = [0;           -tau_lim; -tau_lim; -tau_psi_lim];
u_max = [2*1.0*9.81;   tau_lim;  tau_lim;  tau_psi_lim];   % thrust up to 2*m*g

% Model parameters (must match Nonlinear_Sym / linearize.m)
pm = 1.0;  pg = 9.81;  pIxx = 0.01;  pIyy = 0.01;  pIzz = 0.02;
pAx = 0.25;  pAy = 0.25;  pAz = 0.40;

if use_nl && isempty(which('f_fun_gen'))
    warning('run_sac:nof', 'f_fun_gen.m not found (run linearize.m with lz_check = true). Falling back to the linear plant.');
    use_nl = false;
end
if ~isempty(mm)
    use_nl = true;
    pm = mm.p(1);  pg = mm.p(2);  pIxx = mm.p(3);  pIyy = mm.p(4);  pIzz = mm.p(5);
    pAx = mm.p(6); pAy = mm.p(7); pAz = mm.p(8);
end
%% --- Helpers (anonymous functions) ---------------------------------------
rotC   = @(C,p) [ cos(p)*C(1,:) + sin(p)*C(2,:); ...
                 -sin(p)*C(1,:) + cos(p)*C(2,:); C(3,:); C(4,:) ];
cmdf   = @(d,p) [ cos(p)*d(1) + sin(p)*d(2); -sin(p)*d(1) + cos(p)*d(2); d(3); d(4) ];
ramp   = @(t) 0.5*(1 - cos(pi*min(max((t - t_on)/t_ramp, 0), 1)));
d_world = @(t) cmd_amp*ramp(t);

%% 1. Reference model G_m (yaw-frame deviation outputs) ---------------------
s = tf('s');
wn_xy = 2.0;  zeta_xy = 0.95;
wn_z  = 4.0;  zeta_z  = 0.95;
wn_p  = 4.0;  zeta_p  = 0.95;
% Same second-order models as before, written with states [y; ydot] per
% channel so the model can be started from the plant's initial output.
wn_v = [wn_xy wn_xy wn_z wn_p];
zt_v = [zeta_xy zeta_xy zeta_z zeta_p];
sm = cell(1,4);
for i = 1:4
    sm{i} = ss([0 1; -wn_v(i)^2 -2*zt_v(i)*wn_v(i)], [0; wn_v(i)^2], [1 0], 0);
end
sys_m = blkdiag(sm{1}, sm{2}, sm{3}, sm{4});   % 8 states: [x_b xd_b y_b yd_b z zd psi psid]

Am = sys_m.A;  Bm = sys_m.B;  Cm = sys_m.C;
Ah = H_ss.A;   Bh = H_ss.B;   Ch = H_ss.C;

%% 2. Adaptive gains --------------------------------------------------------
% Regressor rx = [ey (4); x_m (8); r (4)] -> 16 entries.
wp = [ 5*ones(1,4),  0.5*ones(1,8),  0.5*ones(1,4)];
wi = [20*ones(1,4),  1.0*ones(1,8),  1.0*ones(1,4)];
Gamma_p = diag(wp)*10000;
Gamma_i = diag(wi)*100000;

%% 3. Sizes, state packing, storage -----------------------------------------
nx_p = size(A_num,1);   ny = size(Ck,1);   nu = size(B_num,2);
nxm  = size(Am,1);      nxh = size(Ah,1);
nrx  = ny + nxm + ny;

i_xp = 1:nx_p;                       % plant state (absolute X if use_nl, else dx)
i_xm = nx_p + (1:nxm);
i_xh = nx_p + nxm + (1:nxh);
i_Ki = nx_p + nxm + nxh + (1:nu*nrx);

xp_init = x_p0;
if use_nl, xp_init = Xr(:,idx_lin(1)) + x_p0; end
xm0 = zeros(nxm,1);
if match_ic
    c1 = cos(psi_k(1));  s1 = sin(psi_k(1));
    yp0  = rotC(Ck, psi_k(1))*x_p0;                                   % initial deviation output
    ypd0 = [c1*x_p0(2)+s1*x_p0(4); -s1*x_p0(2)+c1*x_p0(4); x_p0(6); x_p0(12)];
    xm0([1 3 5 7]) = yp0;                                            % model outputs
    xm0([2 4 6 8]) = ypd0;                                           % model output rates
end
z = [xp_init; xm0; zeros(nxh,1); zeros(nu*nrx,1)];
nl_flag = double(use_nl);            % 1: subtract Xr to get the deviation

Y_plant = zeros(ny,N);  Y_ref = zeros(ny,N);  Y_h = zeros(ny,N);
E_aug   = zeros(ny,N);  U_dev = zeros(nu,N);  U_tilde = zeros(nu,N);
K_norm  = zeros(1,N);   R_cmd = zeros(ny,N);  Y_world = zeros(ny,N);
Y_des_world = zeros(ny,N);  Y_ref_world = zeros(ny,N);
X_act = zeros(nx_p,N);  U_app = zeros(nu,N);  U_nom = zeros(nu,N);

c_rk = [0 0.5 0.5 1];                % RK4 tableau
w_rk = [1 2 2 1]/6;
Kst  = zeros(numel(z),4);

%% 4. Simulation loop -------------------------------------------------------
if use_nl, pmode = 'NONLINEAR'; else, pmode = 'linearized LTV'; end
fprintf('Simulating SAC closed loop (%s plant): %d samples x %d RK4 sub-steps...\n', pmode, N, n_sub);
tic;
aborted = false;
last = N;
for k = 1:N
    tk   = t_vec(k);
    Ck_k = rotC(Ck, psi_k(k));
    rk   = cmdf(d_world(tk), psi_k(k));
    Xr_k = Xr(:,idx_lin(k));
    ur_k = ur(:,idx_lin(k));

    % ---- log sample k (state before stepping) ----
    xp   = z(i_xp);  xm = z(i_xm);  xh = z(i_xh);
    Ki   = reshape(z(i_Ki), nu, nrx);
    xdev = xp - nl_flag*Xr_k;
    ym = Cm*xm;  yh = Ch*xh;  yp = Ck_k*xdev;
    ey = ym - (yp + yh);
    rx = [ey; xm; rk];
    ut = ey*min(rx.'*Gamma_p*rx, q_max) + Ki*rx;
    uc = min(max(ur_k + E_u*ut, u_min), u_max);

    Y_plant(:,k) = yp;  Y_ref(:,k) = ym;  Y_h(:,k) = yh;
    E_aug(:,k)   = ey;  U_tilde(:,k) = ut;  U_dev(:,k) = uc - ur_k;
    R_cmd(:,k)   = rk;  K_norm(k) = norm(Ki,'fro');
    X_act(:,k)   = Xr_k + xdev;
    U_app(:,k)   = uc;  U_nom(:,k) = ur_k;

    Y_world(:,k)     = Ck*X_act(:,k);
    Y_des_world(:,k) = Ck*Xr_k + d_world(tk);
    cps = cos(psi_k(k));  sps = sin(psi_k(k));
    ymw = ym;  ymw(1:2) = [cps -sps; sps cps]*ym(1:2);     % yaw frame -> world
    Y_ref_world(:,k) = Ck*Xr_k + ymw;

    if ~all(isfinite(z)) || norm(xdev,inf) > blow_up || K_norm(k) > blow_up
        warning('run_sac:diverged', 'Diverged at k=%d (t=%.2f s). Lower gain_scale or raise sigma_leak.', k, tk);
        aborted = true;  last = k;  break;
    end
    if k == N, last = N; break; end

    % ---- advance one sample interval with RK4 ----
    Ak   = A_num(:,:,k);  Bk = B_num(:,:,k);
    Xr_1 = Xr(:,idx_lin(k+1));  ur_1 = ur(:,idx_lin(k+1));
    Tk   = t_vec(k+1) - t_vec(k);
    h    = Tk/n_sub;
    for j = 1:n_sub
        t0 = tk + (j-1)*h;
        for st = 1:4
            ts = t0 + c_rk(st)*h;
            fr = (ts - tk)/Tk;                          % position inside the segment
            Xr_t = (1-fr)*Xr_k + fr*Xr_1;               % nominal state/input (linear interp.)
            ur_t = (1-fr)*ur_k + fr*ur_1;
            if st == 1, zs = z; else, zs = z + c_rk(st)*h*Kst(:,st-1); end
            rs = cmdf(d_world(ts), psi_k(k));

            xp   = zs(i_xp);  xm = zs(i_xm);  xh = zs(i_xh);
            Ki   = reshape(zs(i_Ki), nu, nrx);
            xdev = xp - nl_flag*Xr_t;
            ey = Cm*xm - Ck_k*xdev - Ch*xh;
            rx = [ey; xm; rs];
            ut = ey*min(rx.'*Gamma_p*rx, q_max) + Ki*rx;     % virtual input
            u_app = min(max(ur_t + E_u*ut, u_min), u_max);   % physical input
            dKi = (ey*rx.')*Gamma_i - sigma_leak*Ki;

            if use_nl
                dxp = f_fun_gen(xp, u_app, pm, pg, pIxx, pIyy, pIzz, pAx, pAy, pAz);
                dxp = dxp(:);
            else
                dxp = Ak*xp + Bk*(u_app - ur_t);
            end
            Kst(:,st) = [dxp; Am*xm + Bm*rs; Ah*xh + Bh*ut; dKi(:)];
        end
        z = z + h*(Kst*w_rk.');
    end
end
if aborted, status = 'ABORTED'; else, status = 'completed'; end
fprintf('Simulation %s in %.3f s.\n', status, toc);

%% 5. Diagnostics -----------------------------------------------------------
idx = 1:last;
true_err = Y_ref(:,idx) - Y_plant(:,idx);        % model vs real plant output
X_err    = X_act(:,idx) - Xr(:,idx_lin(idx));    % actual state minus planned state
U_err    = U_app(:,idx) - U_nom(:,idx);
pos_i    = [1 3 5];  ang_i = [7 9 11];
fprintf('Max |ey| (augmented error) per channel : %s\n', mat2str(max(abs(E_aug(:,idx)),[],2).',3));
fprintf('Max |ym - yp| (true error) per channel : %s\n', mat2str(max(abs(true_err),[],2).',3));
fprintf('Position error  X - Xr  [x y z]  max   : %s m\n',   mat2str(max(abs(X_err(pos_i,:)),[],2).',3));
fprintf('Position error  X - Xr  [x y z]  RMS   : %s m\n',   mat2str(sqrt(mean(X_err(pos_i,:).^2,2)).',3));
fprintf('Angle error [phi theta psi]      max   : %s rad\n', mat2str(max(abs(X_err(ang_i,:)),[],2).',3));
fprintf('Input deviation u - ur  [T tph tth tps] max: %s\n', mat2str(max(abs(U_err),[],2).',3));
fprintf('Final ||K_i||_F = %.3g\n', K_norm(last));
if isempty(mm)
%% 6. Plots -----------------------------------------------------------------
t = t_vec(idx);
names = {'x_b','y_b','z','\psi'};
units = {'m','m','m','rad'};

figure('Name','SAC: deviation tracking (yaw frame)');
for i = 1:4
    subplot(2,2,i);
    plot(t, Y_ref(i,idx), 'k--', t, Y_plant(i,idx), 'LineWidth', 1.5); hold on;
    plot(t, R_cmd(i,idx), ':', 'Color', [0.5 0.5 0.5]);
    title([names{i} ' deviation']); xlabel('Time (s)'); ylabel(units{i});
    legend('ref model','plant','command','Location','best'); grid on;
end

figure('Name','SAC: actual state X vs planned Xr');
sidx  = [1 3 5 7 9 11];
slab  = {'x (m)','y (m)','z (m)','\phi (rad)','\theta (rad)','\psi (rad)'};
for i = 1:6
    subplot(3,2,i);
    plot(t, Xr(sidx(i),idx_lin(idx)), 'k--', t, X_act(sidx(i),idx), 'r', 'LineWidth', 1.5);
    title(slab{i}); xlabel('Time (s)'); grid on;
    legend('X_r (plan)','X (actual)','Location','best');
end

figure('Name','SAC: applied input u vs planned ur');
ulab = {'T (N)','\tau_\phi','\tau_\theta','\tau_\psi'};
for i = 1:4
    subplot(2,2,i);
    plot(t, U_nom(i,idx), 'k--', t, U_app(i,idx), 'b', 'LineWidth', 1.5);
    title(ulab{i}); xlabel('Time (s)'); grid on;
    legend('u_r (plan)','u (applied)','Location','best');
end

figure('Name','SAC: tracking errors X - Xr');
subplot(1,2,1); plot(t, X_err(pos_i,:)'); title('Position error X - X_r');
xlabel('Time (s)'); ylabel('m'); legend('x','y','z','Location','best'); grid on;
subplot(1,2,2); plot(t, X_err(ang_i,:)'); title('Angle error X - X_r');
xlabel('Time (s)'); ylabel('rad'); legend('\phi','\theta','\psi','Location','best'); grid on;

figure('Name','SAC: world-frame output vs desired');
wn = {'x','y','z','\psi'};
for i = 1:4
    subplot(2,2,i);
    plot(t, Y_des_world(i,idx), 'k--', t, Y_ref_world(i,idx), 'g:', ...
         t, Y_world(i,idx), 'r', 'LineWidth', 1.5);
    title(wn{i}); xlabel('Time (s)'); ylabel(units{i});
    legend('desired','ref model','plant','Location','best'); grid on;
end

figure('Name','SAC: errors and adaptive gains');
subplot(1,3,1); plot(t, E_aug(:,idx)'); title('e_y = y_m - (y_p + y_h)');
xlabel('Time (s)'); legend(names,'Location','best'); grid on;
subplot(1,3,2); plot(t, true_err'); title('y_m - y_p (true tracking error)');
xlabel('Time (s)'); legend(names,'Location','best'); grid on;
subplot(1,3,3); plot(t, K_norm(idx)); title('||K_i||_F');
xlabel('Time (s)'); grid on;
end
%% --- Export ---------------------------------------------------------------
assignin('base','X_act',X_act);       assignin('base','U_app',U_app);
assignin('base','X_err',X_err);       assignin('base','U_dev',U_dev);
assignin('base','Y_plant',Y_plant);   assignin('base','Y_ref',Y_ref);
assignin('base','Y_world',Y_world);   assignin('base','t_vec',t_vec);
setappdata(0,'mm_res', struct('X_act',X_act(:,1:last), 'U_app',U_app(:,1:last), ...
                              'K_norm',K_norm(1:last), 'aborted',aborted, 'last',last));