% Author: Hashaam Zafar
% Date: 13/11/2025

clear all; clc;
SN = 10078020;
A = 11;
B = 10;
C = 10;n
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
CF = A * 10;           % Hz
DC = 30;               % %
tau_f = (B + C) * 1e-3;   % ms

%% Q1: Question from Assignment 6 (Zero Marks) 

Hs = zpk([], s5, 1);      % 1/(s+120)
Hs = Hs / dcgain(Hs);

Gp = zpk(s4, [s3, s1, s2 - 1j*w2, s2 + 1j*w2], 1);
Gp = Gp / dcgain(Gp);
Kh = 1 / Ks;

Q1.G = Kp * Gp;           % m/V
Q1.H = Hs;                % m/m

%% Question 2

Ts = 1 / CF; % Finding time
Q2.N  = DC / 100 + 0.5;
Q2.wp = CF / Q2.N; 
Q2.Dh = zpk([], - Q2.wp, CF / Q2.N); % CF / (Q2.N * s + CF);

%% Question 3
Q3.NC = 9; 
coefficients = zeros(1, Q3.NC);

for c = 1:Q3.NC
   coefficients(c) = exp(-Ts*(c - 1) / tau_f);
end
sum = 0;
for c = 1:Q3.NC
    sum = sum + coefficients(c);
end
Q3.W = coefficients / sum;


Q3.Nf = 0.92 * coefficients(2) / (1 - coefficients(2)) - 0.2; % 0.92 B / (1 - B) - 0.2
Q3.N = Q2.N + Q3.Nf; 
Q3.wp = CF / Q3.N; 

%% Question 4
new_DH = zpk([], - Q3.wp, CF / Q3.N); % Setting new Dynamics
K = 0.01; 
TF = feedback(K * Q1.G, Q1.H * new_DH);
Q4.K = 0.8;
Q4.Ts = 1; 
Q4.Ess = 1; 


%% Submit

a7Submit


