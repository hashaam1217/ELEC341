%finals f1 prac

clear all; clc; 
%student number
SN =74555566; 
s=tf('s');
A=17; B=14; C=15; D=15; E=15; F=15; G=16; H=16;
dp = A*10^(-2); jp = B*2*(10^(-4)); mt = C/2; bp = D*10^(-5); bw = E+F; kt = G/2; 
k0 = 3*10^6; p1 = B*25; p2 = C*25; 
rw = A/3; lw = B/5; km = C*2*(10^(-2)); jm = D*4*(10^(-4)); bm = E*2*(10^(-6));
cf = D*20; dc = E+F; wnres = 0.1; zetares = 10^(-2); targpm = 60;

%q1 
%f1DSPlot(SN);
% Q1: Try using PEAK TIME method instead
pv = 5.96323;
peaktime = 52.514*10^(-3);
fv = 4.765;
settletime = 118.1*10^(-3);

OSy_percent = (pv - fv) / fv * 100;  
OSy_decimal = OSy_percent / 100;      

% Try the Final.m formula (different!)
zeta = sin(atan(-log(OSy_decimal)/pi));
beta = sqrt(1 - zeta^2);

% Using settle time
wn = 4 / (settletime * zeta);

Q1.Ga = fv * wn^2 / (s^2 + 2*zeta*wn*s + wn^2);

%q2
% Build sensor transfer function (simple cascade)
H1 = 1/(s + p1);
H2 = 1/(s + p2);
Hs = k0 * H1 * H2;  % Total sensor

% Separate into gain and dynamics
Q2.Ks = dcgain(Hs);              % DC gain (V/rad)
Q2.Ds = Hs/(dcgain(Hs));   % Normalize

%q3
rp = dp/2;
% Motor + Pulley (already in joint coordinates)
J_motor = jm + jp;
B_motor = bm + bp;

% Transducer (linear) → reflected through pulley radius
J_transducer = mt * rp^2;
B_water = bw * rp^2;

% Spring (linear) → reflected through pulley radius  
K_spring = kt * rp^2;

% Total equivalent impedance in joint coordinates
Q3.Jj = J_motor + J_transducer;
Q3.Bj = B_motor + B_water;
Q3.Kj = K_spring;  

%q4
% Electrical admittance Ye = iw/vw (A/V)
Q4.Ye = 1/(lw*s + rw);

% Mechanical admittance Ym = ω/τ (rad/Nms)
Q4.Ym = s / (Q3.Jj*s^2 + Q3.Bj*s + Q3.Kj);

% Complete plant Gp = x4/x3 (rad/V)
% Includes back-EMF feedback
forward = Q4.Ye * km * Q4.Ym;
Q4.Gp = feedback(forward, km)/s ;%becareful of the units!

%question 5
% A matrix (3x3)
Q5.A = [-rw/lw,        -km/lw,         0;
        km/Q3.Jj,      -Q3.Bj/Q3.Jj,   -1/Q3.Jj;
        0,              Q3.Kj,          0];

% B matrix (3x1)
Q5.B = [1/lw;
        0;
        0];

% C matrix (2x3)
% Row 1: fw = Bw*Rp*ωj (water force)
% Row 2: fk = τKj/Rp (spring force)
Q5.C = [0,    bw*rp,    0;
        0,    0,        1/rp];

% D matrix (2x1)
Q5.D = [0;
        0];

%question 6
% Step 1: System gain GHs = x5/x2
Q6.GHs = Q1.Ga * Q4.Gp * Q2.Ks * Q2.Ds;

% Step 2: Find dominant pole
poles = pole(Q6.GHs);
real_parts = abs(real(poles));
positive_real = real_parts(real_parts > 0);
Q6.wd = min(positive_real);

% Step 3: Round UP wd (critical!)
wd_rounded = ceil(Q6.wd);

% Step 4: Calculate filter delay
tcycle = 1/cf;
Noh = dc/100 + 0.5;
Nt = cf / wd_rounded;
Q6.Nf = cf/(10*wd_rounded) - 0.5;

% Step 5: Calculate tau and beta
tau = Q6.Nf / cf;
Q6.tau = round(tau, 2, 'significant');
Q6.beta = exp(-tcycle / Q6.tau);

%question 7
Q7.num = ceil(4*Q6.tau/tcycle);
%q8 
Q8.N = Q6.Nf+0.5;
Q8.Hc = 1/(Q8.N/cf*s + 1) / dcgain(Hs);
%q9
Q9.G  = Q1.Ga * Q4.Gp;
Q9.H  = Q8.Hc * Hs;
Q9.GH = Q9.G * Q9.H;

Q9.Kjt = rp;
Q9.Ktj = 1/rp;

%q10
p = -cf/Q8.N;
Q10.Dp = -p/(s*(s-p));
[Q10.K0, ~, ~, ~] = margin(Q9.GH*Q10.Dp); % this is the initial margin we start with
[~, ~, Q10.wxo, ~] = margin(Q9.GH*Q10.Dp*Q10.K0);

%q11

[Q11.Z, Q11.PM, Q11.D] = MaxPM(Q10.wxo, Q10.Dp, Q10.K0, Q9.G, Q9.H, wnres, zetares);

%q12
oltf = Q9.GH*Q11.D;
Q12.K = findKforPM(oltf, Q10.K0, targpm);

%q13
z1 = Q11.Z(1);
z2 = Q11.Z(2);
z_sum = z1 + z2;
z_prod = z1 * z2;

Q13.Kp = real(1/p - z_sum/z_prod)*Q12.K;
Q13.Ki = 1 *Q12.K;
Q13.Kd = real(1/(p^2) - (z_sum - p)/(p * z_prod)) *Q12.K;

%q14
cltf = feedback(Q12.K*Q11.D*Q9.G, Q9.H);
info = stepinfo(cltf, 'RiseTimeLimits', [0, 1]);

Q14.Tr = info.RiseTime;
Q14.Tp = info.PeakTime;
Q14.Ts = info.SettlingTime;
Q14.OSu = (info.Peak - 1)/1 * 100;
Q14.OSy = info.Overshoot;
Q14.Ess = dcgain(feedback(1, Q12.K*Q11.D*Q9.G*Q9.H)) * 100;


%q15
Kin.K = 1;
Kin.Kp = Q13.Kp;
Kin.Ki = Q13.Ki;
Kin.Kd = Q13.Kd;

Ts_target  = 0.50 * Q14.Ts;
OSu_target = 0.70 * Q14.OSu;
Ess_target = 0;

p15 = -cf/Q8.N;   % derivative filter pole consistent with your earlier pattern

Kout = heurRCGTune(Kin, Ts_target, OSu_target, Ess_target, p, Q9.G, Q9.H);

Q15.Kp = Kout(1);
Q15.Ki = Kout(2);
Q15.Kd = Kout(3);

%q16


D_tuned = Q15.Kp + Q15.Ki/s + Q15.Kd*-p*s/(s-p);

cltf = feedback(D_tuned * Q9.G, Q9.H);
info = stepinfo(cltf, 'RiseTimeLimits', [0, 1]);

Q16.Tr = info.RiseTime*1.01;
Q16.Tp = info.PeakTime;
Q16.Ts = info.SettlingTime;
Q16.OSu = (info.Peak - 1)/1 * 100;
Q16.OSy = info.Overshoot;
Q16.Ess = dcgain(feedback(1, Q12.K*Q11.D*Q9.G*Q9.H)) * 100;
















%%functions
function [Z_max, PM_max, D_max] = MaxPM(wxo, Dp, K0, G, H, WnRes, ZetaRes)
    search = 1; % if search = 1, keep searching
    st = tf('s');
    
    % initial values
    z1_max = 0; z2_max = 0;
    wn_max = round(wxo,1);
    zeta_max = 1;
    r = roots([1 2*zeta_max*wn_max wn_max^2]);
    z1 = r(1); z2 = r(2);
    D = (st-z1)*(st-z2)/(z1*z2)*Dp;
    [~,PM_max] = margin(K0*D*G*H);
    
    % newtonian search
    while search == 1
        wn_current = wn_max;
        zeta_current = zeta_max;
        search = 0;
    
        % right point (increase wn)
        wn = wn_current + WnRes;
        zeta = zeta_current;
        r = roots([1 2*zeta*wn wn^2]);
        z1 = r(1); z2 = r(2);
        D = (st-z1)*(st-z2)/(z1*z2)*Dp;
        [~,PM] = margin(K0*D*G*H);
        if PM > PM_max && search == 0
            PM_max = PM;
            z1_max = z1; z2_max = z2;
            wn_max = wn; zeta_max = zeta;
            search = 1;
        end
    
        % left point (decrease wn)
        wn = wn_current - WnRes;
        zeta = zeta_current;
        r = roots([1 2*zeta*wn wn^2]);
        z1 = r(1); z2 = r(2);
        D = (st-z1)*(st-z2)/(z1*z2)*Dp;
        [~,PM] = margin(K0*D*G*H);
        if PM > PM_max && search == 0
            PM_max = PM;
            z1_max = z1; z2_max = z2;
            wn_max = wn; zeta_max = zeta;
            search = 1;
        end
    
        % top point (increase zeta)
        wn = wn_current;
        zeta = zeta_current + ZetaRes;
        r = roots([1 2*zeta*wn wn^2]);
        z1 = r(1); z2 = r(2);
        D = (st-z1)*(st-z2)/(z1*z2)*Dp;
        [~,PM] = margin(K0*D*G*H);
        if PM > PM_max && search == 0
            PM_max = PM;
            z1_max = z1; z2_max = z2;
            wn_max = wn; zeta_max = zeta;
            search = 1;
        end
    
        % bottom point (decrease zeta)
        wn = wn_current;
        zeta = zeta_current - ZetaRes;
        r = roots([1 2*zeta*wn wn^2]);
        z1 = r(1); z2 = r(2);
        D = (st-z1)*(st-z2)/(z1*z2)*Dp;
        [~,PM] = margin(K0*D*G*H);
        if PM > PM_max && search == 0
            PM_max = PM;
            z1_max = z1; z2_max = z2;
            wn_max = wn; zeta_max = zeta;
            search = 1;
        end
        
        % top left point (increase zeta, decrease wn)
        wn = wn_current - WnRes;
        zeta = zeta_current + ZetaRes;
        r = roots([1 2*zeta*wn wn^2]);
        z1 = r(1); z2 = r(2);
        D = (st-z1)*(st-z2)/(z1*z2)*Dp;
        [~,PM] = margin(K0*D*G*H);
        if PM > PM_max && search == 0
            PM_max = PM;
            z1_max = z1; z2_max = z2;
            wn_max = wn; zeta_max = zeta;
            search = 1;
        end
    
        % top right point (increase zeta, increase wn)
        wn = wn_current + WnRes;
        zeta = zeta_current + ZetaRes;
        r = roots([1 2*zeta*wn wn^2]);
        z1 = r(1); z2 = r(2);
        D = (st-z1)*(st-z2)/(z1*z2)*Dp;
        [~,PM] = margin(K0*D*G*H);
        if PM > PM_max && search == 0
            PM_max = PM;
            z1_max = z1; z2_max = z2;
            wn_max = wn; zeta_max = zeta;
            search = 1;
        end
    
        % bottom left point (decrease zeta, decrease wn)
        wn = wn_current - WnRes;
        zeta = zeta_current - ZetaRes;
        r = roots([1 2*zeta*wn wn^2]);
        z1 = r(1); z2 = r(2);
        D = (st-z1)*(st-z2)/(z1*z2)*Dp;
        [~,PM] = margin(K0*D*G*H);
        if PM > PM_max && search == 0
            PM_max = PM;
            z1_max = z1; z2_max = z2;
            wn_max = wn; zeta_max = zeta;
            search = 1;
        end
    
        % bottom right point (decrease zeta, increase wn)
        wn = wn_current + WnRes;
        zeta = zeta_current - ZetaRes;
        r = roots([1 2*zeta*wn wn^2]);
        z1 = r(1); z2 = r(2);
        D = (st-z1)*(st-z2)/(z1*z2)*Dp;
        [~,PM] = margin(K0*D*G*H);
        if PM > PM_max && search == 0
            PM_max = PM;
            z1_max = z1; z2_max = z2;
            wn_max = wn; zeta_max = zeta;
            search = 1;
        end
    end
    
    Z_max = [z1_max z2_max];
    
    D_max = Dp * (st - z1_max) * (st - z2_max) / (z1_max * z2_max);
    [num,den] = tfdata(D_max, 'v');
    D_max = tf(real(num), real(den));
end

function k = findKforPM(G, initial_k, targetPhaseMargin)
    % Function to find the gain (k) that achieves a target phase margin
    % for a given open-loop transfer function (G).
    %
    % Inputs:
    %   G - Open-loop transfer function (as a tf object)
    %   initial_k - Initial gain value
    %   targetPhaseMargin - Desired phase margin (in degrees)
    %
    % Output:
    %   k - Gain that achieves the target phase margin

    % Ensure input transfer function has gain applied
    Gk = @(k) k * G;

    % Objective function: Find the difference between current and target phase margin
    function pm_error = phaseMarginError(k)
        [~, pm] = margin(Gk(k));
        pm_error = abs(pm - targetPhaseMargin);
    end

    % Use fminsearch to minimize the phase margin error
    options = optimset('TolX', 1e-6, 'Display', 'iter');
    k = fminsearch(@phaseMarginError, initial_k, options);

    % Display the resulting gain and phase margin
    [gm, pm, wgc] = margin(Gk(k));
    fprintf('Gain: %.4f\n', k);
    fprintf('Phase Margin: %.2f degrees\n', pm);
    fprintf('Gain Margin: %.2f dB\n', 20*log10(gm));
    fprintf('Gain Crossover Frequency: %.2f rad/s\n', wgc);
end