%%%%%%%%%%%%%%%%%%%%
% Nodal Analysis (Matrix Method) Tutorial
%
% This is a stand-alone script for solving the RLC circuit problem
% that is used to introduce the Matrix Method for solving 
% nodal analysis problems in ELEC 341.
%
% Calling Syntax:
% nodal
%
% Note: All variables cleared when this is run
%
% Author: Leo Stocco
%%%%%%%%%%%%%%%%%%%%

clear all; clc;

% Parameters
R1 = 200;       % (Ohm)
R2 = 500;       % (Ohm)
R3 = 300;       % (Ohm)
L1 = 1;         % (H)
L2 = 2;         % (H)
C1 = 1e-3;      % (F)
C2 = 2e-3;      % (F)

I  = 5;         % (A)

%%%%%%%%%%%%%%%%
% Solve
%%%%%%%%%%%%%%%%
s = tf('s');                % Laplace operator

% Impedances
zr1   = R1;
zr2   = R2;
zr3   = R3;

zl1   = s*L1;
zl2   = s*L2;

zc1   = 1/(s*C1);
zc2   = 1/(s*C2);

% Admittances
yr1   = 1/zr1;
yr2   = 1/zr2;
yr3   = 1/zr3;

yl1   = 1/zl1;
yl2   = 1/zl2;

yc1   = 1/zc1;
yc2   = 1/zc2;

% Y matrix
% Avoid using the transpose operator in any matrix containing the
% Laplace operator 's' or you may get an unexpected negative sign
Y = [ yr1+yr3+yl1+yc1   -yr3            0
     -yr3                yr2+yr3+yc2   -yr2
      0                 -yr2            yr2+yl2];

% I vector
% Semi-colons allows you to fit a 3x1 vector onto 1 line
% This re-defines I from a scalar to a vector
% Avoid doing this because you may need the scalar later
% If the scalar was named Iin instead of I, problem solved
% I = [0; Iin; 0];
I = [0; I; 0];

% Solve for the voltage vector using matrix multiplication
V = inv(Y) * I;

% Display the 3 transfer functions that result
% Could have just omitted the ";" from the previous command
V

% Compute voltages in the time domain
% v and t are both matrices (2 vectors combined)
[v t] = step(V);

%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%
% Plot all voltages together
figure(1); clf; hold on; grid on;
plot(t, v(:,1), 'k-',  'Linewidth', 3);
plot(t, v(:,2), 'b:',  'Linewidth', 3);
plot(t, v(:,3), 'r--', 'Linewidth', 3);
title('All Node Voltages');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Node 1', 'Node 2', 'Node 3', 'Location', 'East');
set(gca, 'FontSize', 14);

% The y-axis is automatically scaled to include all plots
% Because the voltage at Node 2 is so high, you can't
% see the detail of Node voltages 1&3 so plot them
% again without Node 2 to improve resolution for Nodes 1&3
figure(2); clf; hold on; grid on;
plot(t, v(:,1), 'k-',  'Linewidth', 3);
plot(t, v(:,3), 'r--', 'Linewidth', 3);
title('Node Voltages 1&3 Only');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Node 1', 'Node 3', 'Location', 'NorthEast');
set(gca, 'FontSize', 14);
