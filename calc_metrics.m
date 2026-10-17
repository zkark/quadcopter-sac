function m = calc_metrics(X_act, U_app, Xr_s, ur_s, diverged)
% CALC_METRICS  Tracking metrics used to compare SAC, PID, LQR, LQI.
%   X_act, Xr_s : 12 x K   actual and planned state at the same samples
%   U_app, ur_s : 4  x K   applied and planned input at the same samples
%   diverged    : true if the run was aborted (crash / blow-up)
%
% The same function is used for every controller so the numbers are comparable.

E   = X_act - Xr_s;
pos = [1 3 5];
ang = [7 9 11];
K   = size(E,2);
tail = max(1, round(0.8*K)):K;                 % last 20 % of the run

epos = sqrt(sum(E(pos,:).^2, 1));              % 3-D position error norm per sample
eang = sqrt(sum(E(ang,:).^2, 1));              % angle error norm per sample

m.rms_pos   = sqrt(mean(epos.^2));             % m
m.max_pos   = max(epos);                       % m
m.final_pos = mean(epos(tail));                % m, steady-state error (shows bias)
m.rms_ang   = sqrt(mean(eang.^2));             % rad
m.max_ang   = max(eang);                       % rad
m.rms_dT    = sqrt(mean((U_app(1,:) - ur_s(1,:)).^2));              % N
m.rms_dtau  = sqrt(mean(sum((U_app(2:4,:) - ur_s(2:4,:)).^2, 1)));  % N*m
m.diverged  = logical(diverged);
end
