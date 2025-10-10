% Author: Hashaam Zafar
% Date: 03/10/2025

% Setup
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

% Figure 1 Parameters
M0 = A / 5; 
M1 = B / 10; 
M2 = C / 10; 
M3 = D / 5; 

B20 = E / 2; 
B21 = F / 3; 
B31 = G / 4; 

K0 = A; 
K1 = B; 
K20 = C; 
K32 = D / 3; 

% Figure 2 Parameters
Rw = A / 3; 
Lw = B / 1000; 

Jr = C / 10 * 1 / 1000; 
Br = (D + E) * 1 / 1000; 
Ke = G * 50 / 1000; 
Km = G * 50 / 1000; 


syms s;   

% Figure 1 Impedances
% C = M 
Zc0 = (s * M0)^-1; 
Zc1 = (s * M1)^-1; 
Zc2 = (s * M2)^-1; 
Zc3 = (s * M3)^-1; 

% R = 1 / B
Zr20 = 1 / B20;
Zr21 = 1 / B21; 
Zr31 = 1 / B31; 

% L = 1 / K
% Z_L = s*L = s / K
Zl0 = s / K0; 
Zl1 = s / K1;
Zl20 = s / K20; 
Zl32 = s / K32; 