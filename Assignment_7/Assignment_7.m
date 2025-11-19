% Author: Hashaam Zafar
% Date: 13/11/2025

clear all; clc;
SN = 10078020;
A = 11;
B = 10;
C = 10;
D = 17;
E = 18;
F = 10;
G = 12;
H = 10;
s = tf('s');

%% Calculate pole/zero locations and gains
s1 = A * -1;        % -11
s2 = B * -2;        % -20
w2 = C * 2;         % 20
s3 = D * -3;        % -51
s4 = E * -5;        % -90
s5 = F * -12;       % -120
Kp = A^-1 * 3;      % 3/11
Ks = B^-1 * 5;      % 0.5 (not explicitly used)

%% Parameters
CF    = A * 10;           % Control Frequency in Hz
DC    = 30;               % Duty Cycle in %
tau_f = (B + C) * 1e-3;   % Filter time constant in seconds

%% Q1: Open-Loop TFs

Hs = zpk([], s5, 1);      % 1/(s+120)
Hs = Hs / dcgain(Hs);

Gp = zpk(s4, [s3, s1, s2 - 1j*w2, s2 + 1j*w2], 1);
Gp = Gp / dcgain(Gp);

Q1.G = Kp * Gp;           % m/V
Q1.H = Hs;                % m/m

%% Q2: Delay Overhead

Ts    = 1 / CF;
Q2.N  = DC/100 + 0.5;
Q2.wp = 1 / (Q2.N * Ts);
Q2.Dh = zpk([], -Q2.wp, Q2.wp);

%% Q3: FIR Filter

k = 0;
coeffs_raw = [];
w0 = 1.0;
while true
    w_k = exp(-k * Ts / tau_f);
    coeffs_raw = [coeffs_raw, w_k];
    k = k + 1;
    w_next = exp(-k * Ts / tau_f);
    if w_next < 0.02 * w0
        break;
    end
end
coeffs  = coeffs_raw / sum(coeffs_raw);
Q3.NC   = length(coeffs);
Q3.W    = coeffs;
Q3.Nf   = (tau_f / Ts) * exp(-Ts / tau_f);
Q3.N    = Q2.N + Q3.Nf;
Q3.wp   = 1 / (Q3.N * Ts);

%% Q4: P-Control

Q4.Dh = tf(Q3.wp, [1, Q3.wp]);  % Dh(s) with FIR delay

OS_target = 2;   % %
OSmin = 1;
OSmax = 3;

K_low  = 0.1;
K_high = 10.0;

evalOS = @(K) stepinfo(feedback(K * Q1.G, Q1.H * Q4.Dh)).Overshoot;

OS_low  = evalOS(K_low);
OS_high = evalOS(K_high);
iter_expand = 0;
while ~((OS_low  <= OS_target && OS_high >= OS_target) || ...
        (OS_low  >= OS_target && OS_high <= OS_target)) && ...
        iter_expand < 20
    K_low  = K_low  / 2;
    K_high = K_high * 2;
    OS_low  = evalOS(K_low);
    OS_high = evalOS(K_high);
    iter_expand = iter_expand + 1;
end

for it = 1:40
    K_mid = 0.5 * (K_low + K_high);
    OS_mid = evalOS(K_mid);
    if OS_mid > OS_target
        K_high = K_mid;
        OS_high = OS_mid;
    else
        K_low = K_mid;
        OS_low = OS_mid;
    end
end

Q4.K   = K_mid;
CL_P   = feedback(Q4.K * Q1.G, Q1.H * Q4.Dh);
info_P = stepinfo(CL_P);
Q4.Ts  = info_P.SettlingTime;
yfinal = dcgain(CL_P);
Q4.Ess = (1 - yfinal) * 100;

%% Q5: PD Dynamics

sys_openloop  = Q1.G * (Q1.H * Q4.Dh);
poles_sys     = pole(sys_openloop);
[~, idx]      = max(real(poles_sys));
dominant_pole = poles_sys(idx);
wz            = -dominant_pole;
Q5.D          = tf([1/wz, 1], 1);

%% Q6: PD Control

K_range = logspace(-2, 4, 2000);
K_master_PD = NaN;
OSu_PD      = NaN;
Ts_PD       = NaN;
for K_test = K_range
    CL_PD   = feedback(K_test * Q1.G * Q5.D, Q1.H * Q4.Dh);
    info_PD = stepinfo(CL_PD);
    OSu     = info_PD.Overshoot;
    if OSu > OSmin && OSu < OSmax
        K_master_PD = K_test;
        OSu_PD      = OSu;
        Ts_PD       = info_PD.SettlingTime;
        break;
    end
end
Q6.Kp = 1;
Q6.Kd = 1 / wz;

%% Q7: PD Metrics

CL_PD   = feedback(K_master_PD * Q1.G * Q5.D, Q1.H * Q4.Dh);
info_PD = stepinfo(CL_PD);
Q7.Ts   = info_PD.SettlingTime;
yfinal_PD = dcgain(CL_PD);
Q7.Ess  = (1 - yfinal_PD) * 100;

%% Q8: PI Dynamics

wi   = wz;
Q8.D = tf([1/wi, 1], [1, 0]);

%% Q9: PI Control

K_range = logspace(-2, 4, 2000);
K_master_PI = NaN;
OSu_PI      = NaN;
Ts_PI       = NaN;
for K_test = K_range
    CL_PI   = feedback(K_test * Q1.G * Q8.D, Q1.H * Q4.Dh);
    info_PI = stepinfo(CL_PI);
    OSu     = info_PI.Overshoot;
    if OSu > OSmin && OSu < OSmax
        K_master_PI = K_test;
        OSu_PI      = OSu;
        Ts_PI       = info_PI.SettlingTime;
        break;
    end
end
Q9.Kp = 1;
Q9.Ki = wi;

%% Q10: PI Metrics

CL_PI   = feedback(K_master_PI * Q1.G * Q8.D, Q1.H * Q4.Dh);
info_PI = stepinfo(CL_PI);
Q10.Ts  = info_PI.SettlingTime;
yfinal_PI = dcgain(CL_PI);
Q10.Ess = (1 - yfinal_PI) * 100;

%% Summary
fprintf('\n========== SUMMARY ==========\n');
fprintf('P-Control:  Ts=%.3fs, Ess=%.2f%%\n', Q4.Ts,  Q4.Ess);
fprintf('PD-Control: Ts=%.3fs, Ess=%.2f%%\n', Q7.Ts,  Q7.Ess);
fprintf('PI-Control: Ts=%.3fs, Ess=%.2f%%\n', Q10.Ts, Q10.Ess);

%% Submit
a7Submit