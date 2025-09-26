% Author: Hashaam Zafar
% Date: 24/09/2025



Q2.OSy = 31.1 / 22.35 * 100 - 100;
Q2.Tr = 130e-3; 

syms Zeta;
eq1 = (log(Q2.OSy / 100) == - pi *  Zeta / sqrt(1 - Zeta^2));
Zeta = solve(eq1, Zeta);
Zeta = double(Zeta);


Q2.zeta = Zeta; 

syms freqr;
eq1 = (Q2.Tr == 1 / (freqr * sqrt(1 - Zeta^2)) * (pi - atan(sqrt(1 - Zeta^2) / Zeta)));
freqr = solve(eq1, freqr);
freqr = double(freqr);

Q2.wnr = freqr;  
Q2.G = tf(Q2.wnr^2, [1, 2 * Q2.zeta * Q2.wnr, Q2.wnr^2]); 
% a2Submit

% t = 0:0.001:1;  
%Z = Zeta; 
%wn = wnr;
Kdc = 1;
y = tf(Q2.wnr^2, [1, 2 * Q2.wnr * Q2.zeta, wn^2]) * Kdc;
step(Q2.G * Kdc, t);
a2DSPlot(10078020, 2)
