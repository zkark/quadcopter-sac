%% Quadcopter symbolic model and linearization (gyro term neglected)
% State:  X = [x  xdot  y  ydot  z  zdot  phi  phidot  theta  thetadot  psi  psidot]'
% Input:  u = [T  tau_phi  tau_theta  tau_psi]'

clear; 
clc;

%% ---------- Symbols ----------
% States
syms x xdot y ydot z zdot real
syms phi phidot theta thetadot psi psidot real

% Inputs
syms T tau_phi tau_theta tau_psi real

% Parameters
syms m g real                 % mass, gravity
syms Ixx Iyy Izz real         % inertias (principal)
syms Ax Ay Az real            % linear drag coeffs (optional, set 0)

%% ---------- Convenience ----------
Sphi = sin(phi);   Cphi = cos(phi);
Sthe = sin(theta); Cthe = cos(theta);
Spsi = sin(psi);   Cpsi = cos(psi);

% Rotation matrix (body->inertial)
R = [ Cpsi*Cthe,  Cpsi*Sthe*Sphi - Spsi*Cphi,  Cpsi*Sthe*Cphi + Spsi*Sphi;
      Spsi*Cthe,  Spsi*Sthe*Sphi + Cpsi*Cphi,  Spsi*Sthe*Cphi - Cpsi*Sphi;
      -Sthe,      Cthe*Sphi,                    Cthe*Cphi                 ];

% Mapping Euler rates -> body rates
Weta = [ 1,        0,         -Sthe;
         0,      Cphi,   Cthe*Sphi;
         0,     -Sphi,   Cthe*Cphi ];

Weta_inv = [ 1, Sphi*Sthe/Cthe,  Cphi*Sthe/Cthe;
             0,        Cphi,            -Sphi;
             0,  Sphi/Cthe,       Cphi/Cthe ];

% State and input vectors
X  = [x; xdot; y; ydot; z; zdot; phi; phidot; theta; thetadot; psi; psidot];
u  = [T; tau_phi; tau_theta; tau_psi];

%% ---------- Translational dynamics ----------
grav = [0; 0; -g];
e3 = [0;0;1];
a_thrust = (T/m) * R * e3;

% Linear drag in inertial frame
D = diag([Ax Ay Az]);
v_inertial = [xdot; ydot; zdot];
a_drag = - (D/m) * v_inertial;

acc_lin = grav + a_thrust + a_drag;   % [xddot; yddot; zddot]
xddot = acc_lin(1);
yddot = acc_lin(2);
zddot = acc_lin(3);

%% ---------- Rotational dynamics (gyro term neglected) ----------
nu = Weta * [phidot; thetadot; psidot];   % body rates [p;q;r]
p = nu(1); q = nu(2); r = nu(3);

I = diag([Ixx Iyy Izz]);

% Body-rate dynamics: I*nu_dot + nu x (I*nu) = tauB
skew = @(v)[   0   -v(3)  v(2);
              v(3)   0   -v(1);
             -v(2)  v(1)   0  ];

tauB = [tau_phi; tau_theta; tau_psi];
nudot = I \ ( tauB - skew(nu)*(I*nu) );   % no gyro term

% Map to Euler accelerations
qAng   = [phi; theta; psi];
qAng_d = [phidot; thetadot; psidot];

dWeta_inv_dt = sym(zeros(3));
for k = 1:3
    dWeta_inv_dqk = jacobian(Weta_inv(:), qAng(k));
    dWeta_inv_dqk = reshape(dWeta_inv_dqk, 3, 3);
    dWeta_inv_dt  = dWeta_inv_dt + dWeta_inv_dqk * qAng_d(k);
end

etaddot  = dWeta_inv_dt * nu + Weta_inv * nudot;
phiddot    = simplify(etaddot(1));
thetaddot  = simplify(etaddot(2));
psiddot    = simplify(etaddot(3));

%% ---------- Nonlinear state derivative f(X,u) ----------
f = [ xdot;
      xddot;
      ydot;
      yddot;
      zdot;
      zddot;
      phidot;
      phiddot;
      thetadot;
      thetaddot;
      psidot;
      psiddot ];

f = simplify(f);

%% ---------- Symbolic Jacobians ----------
A = jacobian(f, X);
B = jacobian(f, u);

A = simplify(A);
B = simplify(B);

%% ---------- Push to workspace ----------
assignin('base','f_sym',f);
assignin('base','A_sym',A);
assignin('base','B_sym',B);
assignin('base','X_sym',X);
assignin('base','u_sym',u);

disp('Symbolic model (no gyro term) built: f_sym, A_sym, B_sym, X_sym, u_sym');
