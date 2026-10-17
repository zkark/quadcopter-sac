%% Flat-output trajectory generation for quadrotor (consistent with Nonlinear_Sym.m)
% Run AFTER Nonlinear_Sym.m. This script does NOT clear the workspace, and it
% uses unique variable names (tg_ prefix) so it does not overwrite the symbolic
% variables (m, g, I, R, A, B, T, phi, ...) created by Nonlinear_Sym.m.
%
% States: X = [x xdot y ydot z zdot phi phidot theta thetadot psi psidot]'
% Inputs: u = [T tau_phi tau_theta tau_psi]'
% Model convention: z up, grav = [0;0;-g], R = Rz(psi)*Ry(theta)*Rx(phi) (body->inertial),
% states 8/10/12 are EULER RATES, torques from body-rate dynamics incl. gyro term.
%
% Euler angles are computed DIRECTLY from the thrust direction b3 and the yaw psi
% so that R(phi,theta,psi)*e3 = b3 holds exactly for the model's Euler convention.
% Rates and torques come from differentiating the Euler angles on a fine internal
% grid (accurate also at the end points), then sub-sampling to N points.

tg_Tf = 10;   % total time

%% --- Boundary conditions ---
tg_p0   = [0; 0; 0.5];   % start position
tg_pT   = [2; 5; 1.0];   % end position
tg_psi0 = 0;             % yaw start
tg_psiT = pi/2;          % yaw end
% pos/vel/acc/jerk are zero at start and end (7th order polynomial)

%% --- Polynomial coefficients (ascending order c0..c7) ---
% Boundary-condition matrix (same for every axis): derivatives 0..3 at t=0 and t=Tf
tg_M = zeros(8,8);
for tg_d = 0:3
    for tg_i = tg_d:7
        tg_coef = factorial(tg_i)/factorial(tg_i-tg_d);
        tg_M(tg_d+1, tg_i+1) = tg_coef * 0^(tg_i-tg_d);       % at t = 0
        tg_M(tg_d+5, tg_i+1) = tg_coef * tg_Tf^(tg_i-tg_d);   % at t = Tf
    end
end
tg_bc = @(a0,aT) [a0;0;0;0;aT;0;0;0];
tg_cx   = (tg_M \ tg_bc(tg_p0(1), tg_pT(1))).';
tg_cy   = (tg_M \ tg_bc(tg_p0(2), tg_pT(2))).';
tg_cz   = (tg_M \ tg_bc(tg_p0(3), tg_pT(3))).';
tg_cpsi = (tg_M \ tg_bc(tg_psi0,  tg_psiT )).';

% polyval needs descending order
tg_px = fliplr(tg_cx);   tg_py = fliplr(tg_cy);
tg_pz = fliplr(tg_cz);   tg_ps = fliplr(tg_cpsi);

%% --- Parameters (must match the model and linearize.m) ---
tg_m  = 1.0;  tg_g = 9.81;
tg_I  = diag([0.01 0.01 0.02]);     % Ixx Iyy Izz
tg_D  = diag([0.25 0.25 0.40]);     % drag coeffs Ax Ay Az (same as linearize.m)
tg_e3 = [0;0;1];

%% --- Time grids: output grid (N points) and fine internal grid ---
N       = 200;                      % output samples
tg_fine = 10;                       % internal refinement factor
tg_Nf   = (N-1)*tg_fine + 1;
tg_tf   = linspace(0, tg_Tf, tg_Nf);
tg_dt   = tg_tf(2) - tg_tf(1);
tg_sl   = 1:tg_fine:tg_Nf;          % indices of the output samples
tt      = tg_tf(tg_sl);             % 1 x N

%% --- Flat outputs on the fine grid (vectorized) ---
tg_pos = [polyval(tg_px,tg_tf); polyval(tg_py,tg_tf); polyval(tg_pz,tg_tf)];
tg_vel = [polyval(polyder(tg_px),tg_tf); ...
          polyval(polyder(tg_py),tg_tf); ...
          polyval(polyder(tg_pz),tg_tf)];
tg_acc = [polyval(polyder(polyder(tg_px)),tg_tf); ...
          polyval(polyder(polyder(tg_py)),tg_tf); ...
          polyval(polyder(polyder(tg_pz)),tg_tf)];
tg_psi = polyval(tg_ps, tg_tf);                 % 1 x Nf

%% --- Thrust and thrust direction: (T/m)*b3 = a + g*e3 + (D/m)*v ---
tg_b3v = tg_acc + tg_g*(tg_e3*ones(1,tg_Nf)) + (tg_D/tg_m)*tg_vel;
tg_nb  = sqrt(sum(tg_b3v.^2, 1));               % 1 x Nf
tg_Thr = tg_m * tg_nb;
tg_b3  = bsxfun(@rdivide, tg_b3v, tg_nb);       % unit thrust direction

%% --- Euler angles (model convention: R = Rz(psi)*Ry(theta)*Rx(phi)) ---
% b3 = R*e3 = [cpsi*sth*cphi + spsi*sphi; spsi*sth*cphi - cpsi*sphi; cth*cphi]
tg_Phi = asin( sin(tg_psi).*tg_b3(1,:) - cos(tg_psi).*tg_b3(2,:) );
tg_The = atan2( cos(tg_psi).*tg_b3(1,:) + sin(tg_psi).*tg_b3(2,:), tg_b3(3,:) );
tg_eta = [tg_Phi; tg_The; tg_psi];              % 3 x Nf

%% --- Euler rates, body rates, torques ---
tg_etad = zeros(3,tg_Nf);
for tg_i = 1:3
    tg_etad(tg_i,:) = gradient(tg_eta(tg_i,:), tg_dt);
end

% Body rates nu = Weta*eta_dot (same Weta as the model)
tg_nu = [ tg_etad(1,:) - sin(tg_The).*tg_etad(3,:); ...
          cos(tg_Phi).*tg_etad(2,:) + cos(tg_The).*sin(tg_Phi).*tg_etad(3,:); ...
         -sin(tg_Phi).*tg_etad(2,:) + cos(tg_The).*cos(tg_Phi).*tg_etad(3,:) ];

tg_nud = zeros(3,tg_Nf);
for tg_i = 1:3
    tg_nud(tg_i,:) = gradient(tg_nu(tg_i,:), tg_dt);
end

% Torques (model: I*nudot + nu x (I*nu) = tau)
tg_tau = tg_I*tg_nud + cross(tg_nu, tg_I*tg_nu);   % cross works column-wise on 3xN

%% --- Assemble Xr, ur at the N output samples ---
Xr = [ tg_pos(1,tg_sl); tg_vel(1,tg_sl); ...
       tg_pos(2,tg_sl); tg_vel(2,tg_sl); ...
       tg_pos(3,tg_sl); tg_vel(3,tg_sl); ...
       tg_eta(1,tg_sl); tg_etad(1,tg_sl); ...
       tg_eta(2,tg_sl); tg_etad(2,tg_sl); ...
       tg_eta(3,tg_sl); tg_etad(3,tg_sl) ];
ur = [ tg_Thr(tg_sl); tg_tau(:,tg_sl) ];

%% --- Export ---
assignin('base','Xr',Xr);
assignin('base','ur',ur);
assignin('base','tt',tt);

disp('Generated trajectory: Xr (12xN), ur (4xN), tt (1xN) in workspace');