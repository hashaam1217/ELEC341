% Author: Hashaam Zafar
% Date: 24/09/2025

Q4.Ts = 0.71;
Q4.wns = 4 / (Q2.zeta * Q4.Ts); 
Q4.G = tf(Kdc * Q4.wns^2, [1, 2 * Q4.wns * Q2.zeta, Q4.wns^2]);
