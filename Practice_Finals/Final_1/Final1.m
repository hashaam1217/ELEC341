% Practice Final 1
% Author: Hashaam Zafar
% Date: 19/12/2025

%% Setup
SN = 10078020; 
A = 11; 
B = 10;
C = 10; 
D = 17; 
E = 18; 
F = 10; 
G = 12; 
H = 10; 

tf = 's';

% Table 1
Dp = A * 1e-2; % cm
Jp = B * 2e-4; % Nms^2/rad
Mt = C / 2; % Kg
Bp = D * 1e-5;


H1 = 1 / (s + P1);
H2 = 1 / (s + P2);
