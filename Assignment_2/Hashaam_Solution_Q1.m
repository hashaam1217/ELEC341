% Author: Hashaam Zafar
% Date: 24/09/2025

t = 0:0.001:0.20;  
Z = 0.3; 
wn = 42;
Kdc = 1330; 
y = tf(wn, [1, 2 * wn * Z, wn^2]) * Kdc;
% step(y, t);
Q1.FV = 31.7; 
Q1.OSy = 43.3/ 31.7 * 100 - 100;

syms Zeta;
eq1 = (log(Q1.OSy / 100) == - pi *  Zeta / sqrt(1 - Zeta^2));
Zeta = solve(eq1, Zeta);
Zeta = double(Zeta);
Beta = sqrt(1 - Zeta^2);



Q1.zeta = Zeta; 
Q1.Tr = 1 / (Beta * wn) * (pi - atan(Beta / Zeta)); 
Q1.Tp = 79 / 1000;
Q1.wn = wn; 
Q1.Ts = 1 / (Zeta * wn) * log(50 / Beta);


