function res = baseline_sim(ctrl, p, W, opt)
% BASELINE_SIM  Closed-loop simulation of a baseline controller on the NONLINEAR
% quadrotor (f_fun_gen.m), tracking the same plan (Xr, ur) as the SAC loop.
%
%   res = baseline_sim(ctrl, p, W, opt)
%
%   ctrl : 'PID' | 'LQR' | 'LQI'
%            PID : cascaded loops (position PID -> desired thrust vector ->
%                  desired roll/pitch -> attitude PD -> torques), plan feedforward
%            LQR : time-varying (gain-scheduled) LQR along the trajectory,
%                  one gain per linearization sample, u = ur - K*(X - Xr)
%            LQI : same as LQR plus integral action on the outputs [x y z psi]
%   p    : parameters the PLANT really has  [m g Ixx Iyy Izz Ax Ay Az]
%          (the controllers are always designed with W.p_nom, so p ~= W.p_nom
%           IS the model mismatch)
%   W    : struct with fields
%            A_num B_num Ck Xr ur idx_lin t_vec   (from the earlier scripts)
%            x_p0   initial deviation X(0)-Xr(0), 12x1
%            p_nom  nominal parameters used for the design, 1x8
%            u_min, u_max   actuator limits, 4x1
%   opt  : optional struct (all fields optional)
%            n_sub (20)  RK4 sub-steps per sample interval
%            crash_pos (10 m), crash_ang (1.2 rad)   run counts as crashed beyond these
%            xi_max (1)  clamp on the integrator states (anti-windup)
%            Qx, R, Qi   LQR weights
%            Kp_pos, Kd_pos, Ki_pos (3x1), Kp_att, Kd_att (3x1)   PID gains
%            G           cached LQ gains (returned in res.G, pass back to skip the design)
%
%   res fields: X_act, U_app, t, metrics (see calc_metrics.m), diverged, last, G

if nargin < 4 || isempty(opt), opt = struct(); end
n_sub   = getopt(opt, 'n_sub', 20);
crash_p = getopt(opt, 'crash_pos', 10);
crash_a = getopt(opt, 'crash_ang', 1.2);
xi_max  = getopt(opt, 'xi_max', 1);

N   = size(W.A_num, 3);
idx = W.idx_lin;
tv  = W.t_vec;

%% --- Controller parameters ------------------------------------------------
P.Ck     = W.Ck;
P.p_nom  = W.p_nom;
P.u_min  = W.u_min;
P.u_max  = W.u_max;
% PID gains (second-order loops, zeta ~ 0.95; outer ~2 rad/s, z 4 rad/s, attitude ~15 rad/s)
P.Kp_pos = getopt(opt, 'Kp_pos', [ 4;   4;  16  ]);
P.Kd_pos = getopt(opt, 'Kd_pos', [ 3.8; 3.8; 7.6]);
P.Ki_pos = getopt(opt, 'Ki_pos', [ 2;   2;   8  ]);
P.Kp_att = getopt(opt, 'Kp_att', [225; 225; 36  ]);
P.Kd_att = getopt(opt, 'Kd_att', [ 27;  27; 11  ]);

%% --- LQ gain design (nominal model, frozen at every sample) ---------------
G = getopt(opt, 'G', struct());
if ~strcmp(ctrl, 'PID') && ~isfield(G, ctrl)
    Qx = getopt(opt, 'Qx', diag([100 20 100 20 100 20  5 1  5 1  100 5]));
    Rm = getopt(opt, ['R_' ctrl], getopt(opt, 'R', diag([1 50 50 50])));   % R_LQR / R_LQI override R
    Qi = getopt(opt, 'Qi', diag([20 20 20 20]));
    Kx = zeros(4, 12, N);
    Ki = zeros(4, 4,  N);
    for k = 1:N
        Ak = W.A_num(:,:,k);  Bk = W.B_num(:,:,k);
        if strcmp(ctrl, 'LQI')
            Aa = [Ak, zeros(12,4); W.Ck, zeros(4)];
            Ba = [Bk; zeros(4,4)];
            Ka = lqr(Aa, Ba, blkdiag(Qx, Qi), Rm);
            Kx(:,:,k) = Ka(:, 1:12);
            Ki(:,:,k) = Ka(:, 13:16);
        else
            Kx(:,:,k) = lqr(Ak, Bk, Qx, Rm);
        end
    end
    G.(ctrl).Kx = Kx;
    G.(ctrl).Ki = Ki;
end
P.G = G;

%% --- Simulation -----------------------------------------------------------
% z = [X (12); xi (4)],  xi' = Ck*(X - Xr)  (integral of the output deviation)
z = [W.Xr(:, idx(1)) + W.x_p0; zeros(4,1)];

X_act = zeros(12, N);
U_app = zeros(4,  N);
c_rk  = [0 0.5 0.5 1];
w_rk  = [1 2 2 1]/6;
diverged = false;
last = N;

for k = 1:N
    tk   = tv(k);
    Xr_k = W.Xr(:, idx(k));
    ur_k = W.ur(:, idx(k));

    X_act(:,k) = z(1:12);
    U_app(:,k) = ctrl_u(z(1:12), z(13:16), Xr_k, ur_k, k, ctrl, P);

    dx = z(1:12) - Xr_k;
    if ~all(isfinite(z)) || any(abs(dx([1 3 5])) > crash_p) || any(abs(dx([7 9])) > crash_a)
        diverged = true;  last = k;  break;
    end
    if k == N, break; end

    Xr_1 = W.Xr(:, idx(k+1));
    ur_1 = W.ur(:, idx(k+1));
    Tk   = tv(k+1) - tk;
    h    = Tk/n_sub;
    for j = 1:n_sub
        t0 = (j-1)*h;                                  % time inside the segment
        Ks = zeros(16, 4);
        for st = 1:4
            fr   = (t0 + c_rk(st)*h)/Tk;
            Xr_t = (1-fr)*Xr_k + fr*Xr_1;              % plan, linear interpolation
            ur_t = (1-fr)*ur_k + fr*ur_1;
            if st == 1, zs = z; else, zs = z + c_rk(st)*h*Ks(:,st-1); end
            u  = ctrl_u(zs(1:12), zs(13:16), Xr_t, ur_t, k, ctrl, P);
            dX = f_fun_gen(zs(1:12), u, p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8));
            Ks(:,st) = [dX(:); W.Ck*(zs(1:12) - Xr_t)];
        end
        z = z + h*(Ks*w_rk.');
        z(13:16) = min(max(z(13:16), -xi_max), xi_max);   % anti-windup clamp
    end
end

%% --- Output ---------------------------------------------------------------
ii = 1:last;
res.X_act    = X_act(:, ii);
res.U_app    = U_app(:, ii);
res.t        = tv(ii);
res.diverged = diverged;
res.last     = last;
res.G        = G;
res.metrics  = calc_metrics(res.X_act, res.U_app, W.Xr(:, idx(ii)), W.ur(:, idx(ii)), diverged);
end

%% =========================================================================
function u = ctrl_u(X, xi, Xr_t, ur_t, k, ctrl, P)
% Saturated physical input u = [T; tau_phi; tau_theta; tau_psi]
if strcmp(ctrl, 'PID')
    u = pid_u(X, xi, Xr_t, ur_t, P);
else
    dx = X - Xr_t;
    du = -P.G.(ctrl).Kx(:,:,k)*dx;
    if strcmp(ctrl, 'LQI')
        du = du - P.G.(ctrl).Ki(:,:,k)*xi;
    end
    u = ur_t + du;
end
u = min(max(u, P.u_min), P.u_max);
end

%% =========================================================================
function u = pid_u(X, xi, Xr_t, ur_t, P)
% Cascaded PID. Uses only NOMINAL parameters (no knowledge of the mismatch).
m  = P.p_nom(1);  g = P.p_nom(2);
Iv = P.p_nom(3:5).';                 % inertias (column)
Dv = P.p_nom(6:8).';                 % drag coefficients (column)

% Acceleration implied by the plan (nominal model): m*a = T*b3 - m*g*e3 - D*v
etar = Xr_t([7 9 11]);   vr = Xr_t([2 4 6]);
b3r  = [ cos(etar(3))*sin(etar(2))*cos(etar(1)) + sin(etar(3))*sin(etar(1)); ...
         sin(etar(3))*sin(etar(2))*cos(etar(1)) - cos(etar(3))*sin(etar(1)); ...
         cos(etar(2))*cos(etar(1)) ];
a_ref = (ur_t(1)/m)*b3r - [0; 0; g] - (Dv/m).*vr;

% Position PID (xi = integral of X - Xr, so the error integral is -xi)
pos = X([1 3 5]);   vel = X([2 4 6]);
a_cmd = a_ref + P.Kp_pos.*(Xr_t([1 3 5]) - pos) + P.Kd_pos.*(vr - vel) - P.Ki_pos.*xi(1:3);

% Desired thrust vector and magnitude
Fd  = m*(a_cmd + [0; 0; g]) + Dv.*vel;
T   = norm(Fd);
b3d = Fd/max(T, 1e-6);

% Desired roll/pitch for the model's Euler convention R = Rz(psi)*Ry(theta)*Rx(phi)
psir = Xr_t(11);
phid = asin(max(min(sin(psir)*b3d(1) - cos(psir)*b3d(2), 1), -1));
thed = atan2(cos(psir)*b3d(1) + sin(psir)*b3d(2), b3d(3));

% Attitude PD (+ planned rates and planned torque as feedforward)
eta  = X([7 9 11]);   etad = X([8 10 12]);
acc  = P.Kp_att.*([phid; thed; psir] - eta) + P.Kd_att.*(Xr_t([8 10 12]) - etad);
u    = [T; ur_t(2:4) + Iv.*acc];
end

%% =========================================================================
function v = getopt(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end