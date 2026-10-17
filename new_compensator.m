%% =========================================================================
% NEW_COMPENSATOR.M  (aligned channels + yaw-rotated outputs + analytic H)
% =========================================================================
% Global Parallel Feedforward Compensator H(s) for the quadrotor.
% Run order: Nonlinear_Sym -> trajectory_generation -> linearize -> new_compensator
%
% Frozen-time plant at sample k (outputs y = [x_b; y_b; z; psi]):
%   x_b,y_b = yaw-frame position  [x_b; y_b] = [c s; -s c]*[x; y],  psi = psi_ref(k)
%   inputs  u = E_u * u'  ->  u' = [tau_theta; -tau_phi; T; tau_psi]
%   relative degrees of Gp: x_b,y_b = 4 ; z,psi = 2
% Ga = Gp + H, H diagonal with relative degree 1  ->  C*B of Ga = C_H*B_H.
%
% ASPR criterion checked at every frozen sample:
%   (1) Ga minimum phase (all transmission zeros in the open LHP)
%   (2) C*B of Ga symmetric positive definite
% These are checks of FROZEN LTI plants only. They do not prove stability of the
% time-varying / nonlinear closed loop.
% =========================================================================

%% --- Options -------------------------------------------------------------
q_xy      = 6.0;     % zeros of Ga for x,y channels are all placed at -q_xy
d_z       = 0.50;    sigma_z   = 2.0;     % H_z   = d_z  /(s+sigma_z)
d_psi     = 0.50;    sigma_psi = 2.0;     % H_psi = d_psi/(s+sigma_psi)
kappa     = 4*q_xy;  % closed-loop stability diagnostic: K = kappa*inv(C*B)
do_spr    = true;    % optional grid check of SPR of the closed loop (slower)
kappa_spr = 100;     % gain used for the SPR grid check
w_freq    = logspace(-2, 3, 200);

%% --- Workspace checks ----------------------------------------------------
if ~exist('A_num','var') || ~exist('B_num','var') || ~exist('Ck','var')
    error('Missing workspace variables! Run linearize.m first.');
end
if ~exist('Xr','var')
    error('Xr not found (needed for the reference yaw). Run trajectory_generation.m and linearize.m first.');
end
if ~exist('Dk','var'), Dk = zeros(size(Ck,1), size(B_num,2)); end
N = size(A_num,3);
if exist('idx_lin','var'), psi_k = Xr(11, idx_lin); else, psi_k = Xr(11, 1:N); end
if numel(psi_k) ~= N
    error('Number of linearization points (%d) does not match the yaw samples (%d).', N, numel(psi_k));
end

s = tf('s');

% Output rotation into the yaw frame (rows of Ck: x, y, z, psi)
rotC = @(C,p) [ cos(p)*C(1,:) + sin(p)*C(2,:); ...
               -sin(p)*C(1,:) + cos(p)*C(2,:); ...
                C(3,:); C(4,:) ];

%% -------------------------------------------------------------------------
% 1. INPUT ALIGNMENT E_u from the high-frequency gain matrix
% -------------------------------------------------------------------------
% Row i uses the first non-zero Markov parameter C*A^(r-1)*B:
%   x_b,y_b: r = 4 -> C*A^3*B ;  z,psi: r = 2 -> C*A*B
A0 = A_num(:,:,1);  B0 = B_num(:,:,1);  C0 = rotC(Ck, psi_k(1));
K_hf = [ C0(1,:)*A0^3*B0; C0(2,:)*A0^3*B0; C0(3,:)*A0*B0; C0(4,:)*A0*B0 ];

[kmax, col_idx] = max(abs(K_hf), [], 2);
if numel(unique(col_idx)) < 4 || any(kmax < 1e-9)
    disp('K_hf ='); disp(K_hf);
    error('Could not build a signed permutation E_u from K_hf (a channel is not driven by a unique input).');
end
E_u = zeros(4,4);
for i = 1:4
    E_u(col_idx(i), i) = sign(K_hf(i, col_idx(i)));
end
Kal = K_hf * E_u;                                  % aligned HF gains (diagonal expected)
offdiag = Kal - diag(diag(Kal));
if max(abs(offdiag(:))) > 0.1*min(abs(diag(Kal)))
    warning('new_compensator:coupling', 'Aligned high-frequency gain matrix is not diagonally dominant.');
end

fprintf('\nAligned HF gains at k=1  [x_b y_b z psi] = %s\n', mat2str(diag(Kal).', 4));
disp('E_u (u = E_u*u'', u''=[tau_theta; -tau_phi; T; tau_psi]) ='); disp(E_u);

%% -------------------------------------------------------------------------
% 2. ANALYTIC DESIGN OF H FOR THE x,y CHANNELS
% -------------------------------------------------------------------------
% Plant (aligned):  Gp ~ a / ( s^3 (s + c) ).  Choose H = n(s)/d(s) such that
%   a*d(s) + s^3 (s+c) n(s) = (s + q)^6      (all zeros of Ga at -q)
a_xy = mean([Kal(1,1), Kal(2,2)]);
c_xy = -mean([A0(2,2), A0(4,4)]);                  % translational drag / m
if a_xy <= 0 || q_xy <= c_xy
    error('Need a_xy > 0 and q_xy > c_xy (a_xy = %.4g, c_xy = %.4g, q_xy = %.4g).', a_xy, c_xy, q_xy);
end

n1 = 6*q_xy - c_xy;
n0 = 15*q_xy^2 - c_xy*n1;
d3 = (20*q_xy^3 - c_xy*n0)/a_xy;
d2 = 15*q_xy^4/a_xy;
d1 =  6*q_xy^5/a_xy;
d0 =    q_xy^6/a_xy;
num_xy = [1 n1 n0];
den_xy = [d3 d2 d1 d0];
if d3 <= 0 || any(real(roots(den_xy)) >= 0)
    error('Designed H_xy is not stable / has non-positive leading denominator coefficient.');
end

h_x   = tf(num_xy, den_xy);
h_y   = tf(num_xy, den_xy);
h_z   = tf(d_z,   [1 sigma_z]);
h_psi = tf(d_psi, [1 sigma_psi]);
h_z   = tf([d_z   0], conv([1 sigma_z],  [1 sigma_z]));    % was tf(d_z,[1 sigma_z])
h_psi = tf([d_psi 0], conv([1 sigma_psi],[1 sigma_psi]));  % was tf(d_psi,[1 sigma_psi])
H_global = [ h_x,   0,   0,     0;
               0, h_y,   0,     0;
               0,   0, h_z,     0;
               0,   0,   0, h_psi ];
H_ss = blkdiag(ss(h_x), ss(h_y), ss(h_z), ss(h_psi));   % 8 states

fprintf('\n=======================================================\n');
fprintf('  GLOBAL PARALLEL FEEDFORWARD COMPENSATOR H(s) CREATED\n');
fprintf('=======================================================\n');
fprintf('Design values: a_xy = %.4g, c_xy = %.4g, q_xy = %.4g\n', a_xy, c_xy, q_xy);
disp('H(s) Transfer Functions (diagonal):');
disp(h_x); disp(h_z); disp(h_psi);

%% -------------------------------------------------------------------------
% 3. VERIFICATION ACROSS ALL TRAJECTORY SAMPLES (frozen-time)
% -------------------------------------------------------------------------
maxRe_zero = -inf(N,1);     % largest real part of the zeros of Ga
minEig_CB  =  inf(N,1);     % min eigenvalue of sym(C*B) of Ga
asym_CB    =  zeros(N,1);   % |CB - CB'|
maxRe_cl   = -inf(N,1);     % closed loop with K = kappa*inv(CB)
minHerm    =  inf(N,1);     % SPR grid check of closed loop
gain_xy    =  zeros(N,2);   % aligned HF gains over k (to see how much they vary)

fprintf('\nVerifying ASPR conditions across %d trajectory points...\n', N);
tic;
for k = 1:N
    Ak = A_num(:,:,k);
    Bk = B_num(:,:,k) * E_u;
    Ck_k = rotC(Ck, psi_k(k));

    Gp_k = ss(Ak, Bk, Ck_k, Dk);
    Ga_k = Gp_k + H_ss;
    [Aa, Ba, Ca, ~] = ssdata(Ga_k);
    n = size(Aa,1);
    CBa = Ca*Ba;

    % (2) CB symmetric positive definite
    asym_CB(k)   = norm(CBa - CBa', 1);
    minEig_CB(k) = min(eig((CBa + CBa')/2));

    % (1) Transmission zeros (CB invertible): zeros = eig of V'*P*A*V, V = null(C)
    if rcond(CBa) > 1e-12
        P = eye(n) - Ba*(CBa\Ca);
        V = null(Ca);
        maxRe_zero(k) = max(real(eig(V'*P*Aa*V)));

        % Diagnostic: static output feedback u = -K*y, K = kappa*inv(CB)
        Kg  = kappa*inv(CBa);
        Acl = Aa - Ba*Kg*Ca;
        maxRe_cl(k) = max(real(eig(Acl)));

        if do_spr
            Gcl = ss(Aa - Ba*(kappa_spr*inv(CBa))*Ca, Ba, Ca, zeros(4));
            Gf  = freqresp(Gcl, w_freq);
            mh  = inf;
            for iw = 1:numel(w_freq)
                Gv = Gf(:,:,iw);
                mh = min(mh, min(real(eig((Gv + Gv')/2))));
            end
            minHerm(k) = mh;
        end
    end

    % aligned HF gains for reporting
    Kk = [Ck_k(1,:)*Ak^3*Bk; Ck_k(2,:)*Ak^3*Bk];
    gain_xy(k,:) = [Kk(1,1), Kk(2,2)];
end
elapsed_time = toc;

%% -------------------------------------------------------------------------
% 4. DIAGNOSTICS & SUMMARY
% -------------------------------------------------------------------------
all_min_phase = all(maxRe_zero < 0);
all_CB_ok     = all(minEig_CB > 0) && all(asym_CB < 1e-8*max(1,max(abs(minEig_CB))));
overall_aspr  = all_min_phase && all_CB_ok;

fprintf('Verification completed in %.3f seconds.\n\n', elapsed_time);
fprintf('--- ASPR Diagnostic Summary (frozen-time, all %d samples) ---\n', N);
fprintf('Aligned x,y HF gain range over k : [%.1f, %.1f]  (design a_xy = %.1f)\n', ...
        min(gain_xy(:)), max(gain_xy(:)), a_xy);
if all_min_phase, str1 = 'PASSED'; else, str1 = 'FAILED'; end
fprintf('Minimum phase (zeros of Ga in LHP): %s (worst zero real part = %.4f)\n', str1, max(maxRe_zero));
if all_CB_ok, str2 = 'PASSED'; else, str2 = 'FAILED'; end
fprintf('C*B symmetric positive definite   : %s (min eig = %.4g, max asymmetry = %.2e)\n', ...
        str2, min(minEig_CB), max(asym_CB));
fprintf('Static feedback K=%.1f*inv(CB)    : worst closed-loop real part = %.4f\n', kappa, max(maxRe_cl));
if do_spr
    fprintf('Closed-loop SPR grid check (kappa=%g, w<=%g): min Hermitian eig = %.3e (grid check, not a proof)\n', ...
            kappa_spr, w_freq(end), min(minHerm));
end

if overall_aspr
    fprintf('\n=======================================================\n');
    fprintf('  ASPR criterion (min-phase + CB>0) holds at all %d frozen samples.\n', N);
    fprintf('  Note: this does not prove stability of the time-varying/nonlinear loop.\n');
    fprintf('=======================================================\n\n');
else
    fprintf('\n=======================================================\n');
    fprintf('  WARNING: ASPR criterion fails at some samples.\n');
    fprintf('=======================================================\n\n');
end

assignin('base', 'E_u', E_u);
assignin('base', 'H_global', H_global);
assignin('base', 'H_ss', H_ss);
assignin('base', 'h_x', h_x);
assignin('base', 'h_y', h_y);
assignin('base', 'h_z', h_z);
assignin('base', 'h_psi', h_psi);
assignin('base', 'psi_k', psi_k);
assignin('base', 'maxRe_zero', maxRe_zero);