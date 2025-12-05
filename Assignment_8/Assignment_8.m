% Author: Hashaam Zafar
% Date: 28/11/2025
% Assignment 8: Optimal Control (PD, PI, PID)

clear all; clc;
SN = 10078020;
A = 11; B = 10; C = 10; D = 17; E = 18; F = 10; G = 12; H = 10;
s = tf('s');

%% Parameters from Assignment 6
s1 = A * -1;        
s2 = B * -2;        
w2 = C * 2;         
s3 = D * -3;        
s4 = E * -5;        
s5 = F * -12;       
Kp = A^-1 * 3;      
Ks = B^-1 * 5;      

%% Parameters from Assignment 7
CF = A * 10;           % Hz = 110 Hz
DC = 30;               % %
tau_f = (B + C) * 1e-3;   % ms = 0.02 s

%% A8 Parameters
tau_pid = D / 2 * 1e-3;  % ms = 8.5 ms
wn_res = 10^-1;          % rad/s
zeta_res = 10^-2;        % pure
Target_PM = 60;          % deg

%% Question 1

%% Q1: A7 Parameters (0 marks)
% Plant transfer function
Gp = zpk(s4, [s3, s1, s2 - 1j*w2, s2 + 1j*w2], 1);
Gp = Gp / dcgain(Gp);  % Normalize plant

% Sensor transfer function
Hs = zpk([], s5, -s5);  % Gain of -s5 = 120 at DC

Q1.G = Kp * Gp;      % m/V
Q1.H = Hs / Ks;      % m/m (sensor with gain factor)

%% Q2: IIR Filter (10 marks)
Ts = 1 / CF;  % Sample time
Q2.beta = exp(-Ts / tau_pid);
Q2.Nf = Q2.beta / (1 - Q2.beta);
Q2.Npid = Q2.Nf + 0.5;
Q2.wp = (2/Ts) * ((1 - Q2.beta) / (1 + Q2.beta));

%% Q3: PD Dynamics (10 marks)
% Create feedback dynamics with IIR filter
Dh_pid = zpk([], -Q2.wp, Q2.wp);

% Find initial gain and crossover
OL = Q1.G * Q1.H * Dh_pid;
[mag, phase, w] = bode(OL, logspace(-2, 3, 5000));
mag = squeeze(mag);
phase = squeeze(phase);

% Find where phase = -120 deg (for 60 deg PM)
target_phase = -180 + Target_PM;
[~, idx] = min(abs(phase - target_phase));
Q3.wxo = w(idx);
Q3.K0 = 1 / mag(idx);

% Find most dominant pole to cancel (least negative real pole)
poles_G = pole(Q1.G);
real_poles = poles_G(imag(poles_G) == 0);
[~, idx_dom] = max(real_poles);  % Most dominant = least negative
dominant_pole = real_poles(idx_dom);

% Optimal zero search - find the magnitude of zero location
% The zero will be placed at -Z (in LHP)
Z_test_range = linspace(abs(dominant_pole) * 0.3, abs(dominant_pole) * 1.5, 50);
best_PM = 0;
best_Z = abs(dominant_pole);

for Z_val = Z_test_range
    D_test = zpk(-Z_val, -Q2.wp, Q2.wp);  % Zero at -Z_val (LHP)
    OL_test = Q3.K0 * Q1.G * Q1.H * D_test;
    [Gm, Pm] = margin(OL_test);
    
    if abs(Pm - Target_PM) < abs(best_PM - Target_PM)
        best_PM = Pm;
        best_Z = Z_val;
    end
end

% Q3.Z should be the MAGNITUDE (positive value)
% The actual zero is at -Q3.Z
Q3.Z = best_Z;
% Q3.D = zpk(-Q3.Z, -Q2.wp, Q2.wp);
num = Q2.wp * [1, Q3.Z];
den = conv([1, Q3.Z], [1, Q2.wp]);

Q3.D = tf(num, den);
Q3.D = 1 / s;
% Verify zero is in LHP
z_D = zero(Q3.D);
if any(real(z_D) > 0)
    warning('Q3.D has zeros in RHP! Zero at: %.4f', z_D);
end

%% Q4: PD Gains (10 marks)
OL_PD = Q1.G * Q1.H * Q3.D;
[mag_pd, phase_pd, w_pd] = bode(OL_PD, logspace(-2, 3, 5000));
mag_pd = squeeze(mag_pd);
phase_pd = squeeze(phase_pd);

% Find gain for target PM
[~, idx_xo] = min(abs(phase_pd - (-180 + Target_PM)));
Q4.K = 1 / mag_pd(idx_xo);

% Normalized gains (K=1)
% PD controller: D(s) = (s + Z) / (s + wp) * wp
% At DC (s=0): D(0) = Z/wp * wp = Z
% For normalized form: Kp + Kd*s
% Evaluate at s=0: Kp = Z / wp * wp = Z
% Coefficient of s: Kd = 1 / wp * wp = 1
[num_d, den_d] = tfdata(Q3.D, 'v');
% num_d = [wp, wp*Z], den_d = [1, wp]
Q4.Kp = num_d(2) / den_d(2);  % wp*Z / wp = Z
Q4.Kd = num_d(1) / den_d(2);  % wp / wp = 1

%% Q5: PD Metrics (10 marks)
CL_PD = feedback(Q4.K * Q1.G * Q3.D, Q1.H);
step_info = stepinfo(CL_PD);
Q5.Ts = step_info.SettlingTime;
Q5.Ess = (1 - dcgain(CL_PD)) * 100;

%% Q6: PI Dynamics (10 marks)
% PI controller: D(s) = (s + Z) / s
% Feedback dynamics same as Q2
Dh_pi = zpk([], -Q2.wp, Q2.wp);

% Create open loop with integrator
OL_pi_init = Q1.G * Q1.H * Dh_pi / s;

[mag_pi, phase_pi, w_pi] = bode(OL_pi_init, logspace(-2, 3, 5000));
mag_pi = squeeze(mag_pi);
phase_pi = squeeze(phase_pi);

[~, idx_pi] = min(abs(phase_pi - target_phase));
Q6.wxo = w_pi(idx_pi);
Q6.K0 = 1 / mag_pi(idx_pi);

% Optimal zero for PI - search around dominant pole
Z_test_range_pi = linspace(abs(dominant_pole) * 0.2, abs(dominant_pole) * 1.2, 50);
best_PM_pi = 0;
best_Z_pi = abs(dominant_pole) * 0.5;

for Z_val = Z_test_range_pi
    D_test = zpk(-Z_val, 0, Z_val);  % Zero at -Z_val, pole at origin
    OL_test = Q6.K0 * Q1.G * Q1.H * Dh_pi * D_test;
    [Gm, Pm] = margin(OL_test);
    
    if abs(Pm - Target_PM) < abs(best_PM_pi - Target_PM)
        best_PM_pi = Pm;
        best_Z_pi = Z_val;
    end
end

Q6.Z = best_Z_pi;
Q6.D = zpk(-Q6.Z, 0, Q6.Z);  % Normalized: (s+Z)/s * Z

%% Q7: PI Gains (10 marks)
OL_PI = Q1.G * Q1.H * Dh_pi * Q6.D;
[mag_pi2, phase_pi2, w_pi2] = bode(OL_PI, logspace(-2, 3, 5000));
mag_pi2 = squeeze(mag_pi2);
phase_pi2 = squeeze(phase_pi2);

[~, idx_xo_pi] = min(abs(phase_pi2 - target_phase));
Q7.K = 1 / mag_pi2(idx_xo_pi);

% Normalized gains: D(s) = (s + Z) / s = Kp + Ki/s
Q7.Kp = 1;
Q7.Ki = Q6.Z;

%% Q8: PI Metrics (10 marks)
CL_PI = feedback(Q7.K * Q1.G * Dh_pi * Q6.D, Q1.H);
step_info_pi = stepinfo(CL_PI);
Q8.Tr = step_info_pi.RiseTime;
[y_pi, t_pi] = step(CL_PI);
Q8.OSu = (max(y_pi) - 1) * 100;

%% Q9: PID Dynamics (10 marks)
% PID: D(s) = (s + Z1)(s + Z2) / (s(s + wp)) * wp
Dh_pid2 = zpk([], -Q2.wp, Q2.wp);

% Initial open loop with integrator
OL_pid_init = Q1.G * Q1.H * Dh_pid2 / s;

[mag_pid, phase_pid, w_pid] = bode(OL_pid_init, logspace(-2, 3, 5000));
mag_pid = squeeze(mag_pid);
phase_pid = squeeze(phase_pid);

[~, idx_pid] = min(abs(phase_pid - target_phase));
Q9.wxo = w_pid(idx_pid);
Q9.K0 = 1 / mag_pid(idx_pid);

% Optimal zeros for PID - search for two zeros in LHP
Z1_range = linspace(abs(dominant_pole) * 0.3, abs(dominant_pole) * 0.9, 20);
Z2_range = linspace(abs(dominant_pole) * 0.7, abs(dominant_pole) * 1.5, 20);
best_PM_pid = 0;
best_Z1 = abs(dominant_pole) * 0.5;
best_Z2 = abs(dominant_pole) * 1.0;

for Z1_val = Z1_range
    for Z2_val = Z2_range
        % Both zeros in LHP: -Z1_val and -Z2_val
        D_test = zpk([-Z1_val, -Z2_val], [0, -Q2.wp], Q2.wp);
        OL_test = Q9.K0 * Q1.G * Q1.H * D_test;
        [Gm, Pm] = margin(OL_test);
        
        if abs(Pm - Target_PM) < abs(best_PM_pid - Target_PM)
            best_PM_pid = Pm;
            best_Z1 = Z1_val;
            best_Z2 = Z2_val;
        end
    end
end

Q9.Z = [best_Z1, best_Z2];
Q9.D = zpk([-Q9.Z(1), -Q9.Z(2)], [0, -Q2.wp], Q2.wp);

%% Q10: PID Gains (10 marks)
OL_PID = Q1.G * Q1.H * Q9.D;
[mag_pid2, phase_pid2, w_pid2] = bode(OL_PID, logspace(-2, 3, 5000));
mag_pid2 = squeeze(mag_pid2);
phase_pid2 = squeeze(phase_pid2);

[~, idx_xo_pid] = min(abs(phase_pid2 - target_phase));
Q10.K = 1 / mag_pid2(idx_xo_pid);

% PID: (s + Z1)(s + Z2) / s = s + (Z1+Z2) + Z1*Z2/s
% = Kd*s + Kp + Ki/s
Q10.Ki = Q9.Z(1) * Q9.Z(2);
Q10.Kp = Q9.Z(1) + Q9.Z(2);
Q10.Kd = 1;

%% Q11: PID Metrics (10 marks)
CL_PID = feedback(Q10.K * Q1.G * Q9.D, Q1.H);
step_info_pid = stepinfo(CL_PID);
Q11.Tr = step_info_pid.RiseTime;
[y_pid, t_pid] = step(CL_PID);
Q11.OSu = (max(y_pid) - 1) * 100;

%% Display Results
fprintf('\n=== ASSIGNMENT 8 RESULTS ===\n\n');

fprintf('Q1: A7 Parameters\n');
fprintf('  G DC gain = %.6f m/V\n', dcgain(Q1.G));
fprintf('  H DC gain = %.6f m/m\n\n', dcgain(Q1.H));

fprintf('Q2: IIR Filter\n');
fprintf('  beta = %.6f (pure)\n', Q2.beta);
fprintf('  Nf = %.4f (pure)\n', Q2.Nf);
fprintf('  Npid = %.4f (pure)\n', Q2.Npid);
fprintf('  wp = %.4f rad/s\n\n', Q2.wp);

fprintf('Q3: PD Dynamics\n');
fprintf('  K0 = %.6f V/m\n', Q3.K0);
fprintf('  wxo = %.4f rad/s\n', Q3.wxo);
fprintf('  Z = %.4f rad/s\n', Q3.Z);
fprintf('  D DC gain = %.6f\n\n', dcgain(Q3.D));

fprintf('Q4: PD Gains\n');
fprintf('  K = %.6f V/m\n', Q4.K);
fprintf('  Kp = %.6f (pure)\n', Q4.Kp);
fprintf('  Kd = %.6f (pure)\n\n', Q4.Kd);

fprintf('Q5: PD Metrics\n');
fprintf('  Ts = %.4f s\n', Q5.Ts);
fprintf('  Ess = %.4f %%\n\n', Q5.Ess);

fprintf('Q6: PI Dynamics\n');
fprintf('  K0 = %.6f V/m\n', Q6.K0);
fprintf('  wxo = %.4f rad/s\n', Q6.wxo);
fprintf('  Z = %.4f rad/s\n', Q6.Z);
fprintf('  D poles: '); disp(pole(Q6.D)');

fprintf('Q7: PI Gains\n');
fprintf('  K = %.6f V/m\n', Q7.K);
fprintf('  Kp = %.6f (pure)\n', Q7.Kp);
fprintf('  Ki = %.6f (pure)\n\n', Q7.Ki);

fprintf('Q8: PI Metrics\n');
fprintf('  Tr = %.4f s\n', Q8.Tr);
fprintf('  OSu = %.4f %%\n\n', Q8.OSu);

fprintf('Q9: PID Dynamics\n');
fprintf('  K0 = %.6f V/m\n', Q9.K0);
fprintf('  wxo = %.4f rad/s\n', Q9.wxo);
fprintf('  Z = [%.4f, %.4f] rad/s\n', Q9.Z(1), Q9.Z(2));
fprintf('  D zeros: '); disp(zero(Q9.D)');
fprintf('  D poles: '); disp(pole(Q9.D)');

fprintf('Q10: PID Gains\n');
fprintf('  K = %.6f V/m\n', Q10.K);
fprintf('  Kp = %.6f (pure)\n', Q10.Kp);
fprintf('  Ki = %.6f (pure)\n', Q10.Ki);
fprintf('  Kd = %.6f (pure)\n\n', Q10.Kd);

fprintf('Q11: PID Metrics\n');
fprintf('  Tr = %.4f s\n', Q11.Tr);
fprintf('  OSu = %.4f %%\n\n', Q11.OSu);

%% Personal Attempt
Assignment_7; 

%% A8 Parameters
tau_pid = D / 2 * 1e-3;  % ms = 8.5 ms
wn_res = 10^-1;          % rad/s
zeta_res = 10^-2;        % pure
Target_PM = 60;          % deg

%% Question 1
% Q1.G remains the same as Assignment 7
Q1.H = Hs * new_DH; 

%% Question 2
dt = 1 / CF; 
Q2.beta = exp(-dt/tau_pid);
Q2.Nf = Q2.beta / (1 - Q2.beta); 
Q2.Npid = Q2.Nf + 0.5; 
Q2.wp = (2/Ts) * ((1 - Q2.beta) / (1 + Q2.beta)); % BiLinear Transform. 
% Zero idea where that came from. 

%% Question 3
D = Q2.wp / (s + Q2.wp);
OpenLoopTF = D * Q1.G * Q1.H; 
Q3.K0 = 8.7838e08;
Q3.wxo = 1; 
Q3.Z = 1; 
Q3.D = 1 / s; 

%% Submit
% a8Submit