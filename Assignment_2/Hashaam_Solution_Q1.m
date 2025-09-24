% Author: Hashaam Zafar
% Date: 24/09/2025

t = 0:0.001:020;  
Z = 0.3; 
B = (1 - Z^2)^(0.5);
wn = 42;
Kdc = 1330; 
y = tf(wn, [1, 2 * wn * Z, wn^2]) * Kdc;
step(y, t)

Q1.FV = Kdc; 
Q1.OSy =;
Q1.zeta = Z; 
Q1.Tr = 1; 
Q1.Tp = 79 / 1000;
Q1.wn = wn; 
Q1.Ts = 0;