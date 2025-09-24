%%%%%%%%%%%%%%%%%%%%
% Op-Amp Circuit Tutorial
%
% This is a stand-alone script for solving the Op-Amp RLC circuit problem
% that is used to introduce Matlab in ELEC 341.
%
% Calling Syntax:
% opamp
%
% Note: All variables cleared when this is run
%
% Author: Leo Stocco
%%%%%%%%%%%%%%%%%%%%

clear all; clc;

% Parameters
R1 = 200;
R2 = 500;
L  = 1;
C  = 1e-3;

% Solve 3 different ways
xfer = useTF(R1, R2, L, C);   % store the xfer fun
useZPK(R1, R2, L, C);
useAlg(R1, R2, L, C);

% Bode Plot
figure(4)
bode(xfer);
grid on;

% Bode Plot with extra info
figure(5)
margin(xfer);
grid on;

%%%%%%%%%%%%%%%%%%
% Step Response
%%%%%%%%%%%%%%%%%%
[ys ts] = step(xfer);

figure(6);
clf;
plot(ts, ys, 'b-', 'LineWidth', 3);
grid on;
title('Step Response');
xlabel('Time (sec)');
ylabel('OP Voltage (V)');

%%%%%%%%%%%%%%%%%%
% Create PNG files
%%%%%%%%%%%%%%%%%%
if 0                % saves you the trouble of commenting out many lines
% Save figures to PNG format
print -dpng -f1 fig1
print -dpng -f2 fig2
print -dpng -f3 fig3
print -dpng -f4 fig4
print -dpng -f5 fig5
end

%%%%%%%%%%%%%%%%%%
% Internal Funcs
%%%%%%%%%%%%%%%%%%
% Compute numerator & denominator
function [num den] = numDen(R1, R2, L, C)
  num = [L*C*R2 L+C*R1*R2 R1+R2];
  den = [L*C*R2 L+C*R1*R2 R1];
end % function

% Solve using tf()
function xfer = useTF(R1, R2, L, C)
  % Compute xfer func
  [num den] = numDen(R1, R2, L, C);
  xfer = tf(num, den);

  % Cancel common numerator & denominator terms
  % This function may reduce accuracy because
  % It DOES NOT always preserve DC gain
  % Use with caution
  xfer = minreal(xfer);

  % Plot impulse response
  figure(1); clf;
  impulse(xfer);
  grid on;
end % function

% Solve using zpk()
function xfer = useZPK(R1, R2, L, C)
  % Compute xfer func
  [num den] = numDen(R1, R2, L, C);
  z = roots(num)
  p = roots(den)
  xfer = zpk(z, p, 1);

  [y t] = impulse(xfer);

  % Plot impulse response
  figure(2); clf;
  h = plot(t, y, 'r--');
  set(h, 'LineWidth', 3);
  grid on;
  title('Dotted Red Impulse Response');
  xlabel('Time (sec)');
  ylabel('OP Voltage (V)');
  set(gca, 'FontSize', 14);
end % function

% Solve algebraically
function xfer = useAlg(R1, R2, L, C)
  % Compute xfer func
  s = tf('s');
  Zr1 = R1;
  Zr2 = R2;
  Zl  = s*L;
  Zc  = 1/(s*C);

  Z1 = Zr1+Zl;
  Z2 = 1/(1/Zr2 + 1/Zc);
  xfer = minreal(1 + Z2/Z1);

  % Plot impulse response
  figure(3); clf;
  impulse(xfer);
  grid on;
end % function
