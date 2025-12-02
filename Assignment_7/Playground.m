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
CF = A * 10;             % Hz
DC = 30;                 % %
tau_f = (B + C) * 1e-3;  % ms

%% Q1: Question from Assignment 6 (Zero Marks) 

Hs = zpk([], s5, 1);          % 1/(s+120)
Hs = Hs / dcgain(Hs);         % normalize dc gain to 1

Gp = zpk(s4, [s3, s1, s2 - 1j*w2, s2 + 1j*w2], 1);
Gp = Gp / dcgain(Gp);         % normalize dc gain to 1
Kh = 1 / Ks;

Q1.G = Kp * Gp;               % m/V
Q1.H = Hs;                    % m/m

%% Question 2

Ts = 1 / CF;                  % Sample period
Q2.N  = DC / 100 + 0.5;
Q2.wp = CF / Q2.N;
Q2.Dh = zpk([], -Q2.wp, CF / Q2.N);   % CF / (Q2.N * s + CF)

%% Question 3

Q3.NC = 9; 
coefficients = zeros(1, Q3.NC);

for c = 1:Q3.NC
   coefficients(c) = exp(-Ts*(c - 1) / tau_f);
end

sum_c = 0;
for c = 1:Q3.NC
    sum_c = sum_c + coefficients(c);
end
Q3.W = coefficients / sum_c;

Bcoef   = coefficients(2);
Q3.Nf   = 0.92 * Bcoef / (1 - Bcoef) - 0.2;  % 0.92 B / (1 - B) - 0.2
Q3.N    = Q2.N + Q3.Nf; 
Q3.wp   = CF / Q3.N; 

%% Question 4  -- P-Control design to meet 1% < OS_r < 3%

% Re-calculated delay dynamics with the FIR delay included
new_DH = zpk([], -Q3.wp, CF / Q3.N);

% Search for a K such that 1% < Overshoot < 3%
K_candidates = linspace(0.5, 2.0, 301);     % search range for K
K_star = NaN;
info_star = [];

for K_try = K_candidates
    TF_try = feedback(K_try * Q1.G, Q1.H * new_DH);
    info   = stepinfo(TF_try);
    OS     = info.Overshoot;   % percent overshoot

    if (OS > 1) && (OS < 3)
        K_star   = K_try;
        info_star = info;
        break;                 % first K that satisfies the RCG
    end
end

if isnan(K_star)
    error('No K found in search range that satisfies 1%% < OS_r < 3%%.');
end

TF = feedback(K_star * Q1.G, Q1.H * new_DH);

% Fill in Q4
Q4.K   = K_star;                      % proportional gain that meets RCG
Q4.Ts  = info_star.SettlingTime;      % settling time from stepinfo
Q4.Ess = (1 - dcgain(TF)) * 100;      % steady-state error in percent

%% Submit
a7Submit
