% Author: Hashaam Zafar
% Date: 12/09/2025

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

s = tf('s');                % Laplace operator

% Setup For Q1
R1 = A * 10;
R2 = D * 10;
L1 = B * 1e-3;
C1 = C * 1e-6; 
L2 = E * 1e-3; 
C2 = F * 1e-6; 

ZR1 = R1;
ZR2 = R2;
ZL1 = s * L1;
ZL2 = s * L2;
ZC1 = 1 / (s * C1);
ZC2 = 1 / (s * C2);
