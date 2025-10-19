% Date: 18/10/2025
% Name: Hashaam Zafar
% Assignment: Project Part 1

clear all; clc; 
%student number
SN =10078020; 
s=tf('s');

A = 11;
B = 10; 
C = 10; 
D = 17; 
E = 18;
F = 10; 
G = 12; 
H = 10; 

% Table 1
Js = (A / 5) * 1e-7; % gcm^2
Ms = B / 4 * 1e-3; % g
Bs = C / 3; % Ns/m
ns = D / (2*E) * 1e2; % turn/cm
Jf = F / 3 * 1e-7; % gcm^2
Bf = G * 1e-3; % mNms
nf = H * 1e2; % deg/cm
Bt = A *1e-3; % mNms
Kt = B*C* 1e-3; % mNm
% nt? 
Bl = D / 5; % Ns/m
Kl = E*F; % N/m
L6 = 4*(G+H)*1e-3; % mm

% Table 2
P1 = A*7; 
P2 = B*800; 
P3 = C*8;
P4 = D*700; 
P5 = E*600; 
P6 = F*50;
P7 = G*500; 
P8 = H*5; 

% Table 3
Rw = A / 2; % Ohms
Lw = B*30 * 1e-6; % uH
Km = C * 1e-3; % mNm/A
Jr = D/15 * 1e-7; % gcm^2
Br = E / 30 * 1e-6; % uNms
Mm = F+G; % Kg

% Table 4
CF = A*B*5; % Hz
DC = C+D+E+F; % %

% Table 5
WnRes = 0.1;
ZetaRes = 10^-2; 
TargPM = 40; % deg
OSu = 60; % Caution: these are funky
Ts = 75; 
% Tr less than final Ts
Ess = 0; 

% Question 1
% Figure 2
G1 = P1 / (s + P2);
G2 = P3 / (s + P4);
G3 = 10^5 / (s + P5);
G4 = P6 / (s + P7);
H1 = 4 / (s + P8); 

Q1.Tr = 0.75 *1e-3;
Q1.Tp = 0.00165;
Q1.Ts = 4.65 * 1e-3;
Q1.OSy = 25; 

p1Submit