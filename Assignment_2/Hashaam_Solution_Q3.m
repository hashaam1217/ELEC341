% Author: Hashaam Zafar
% Date: 24/09/2025



Q3.Tp = 0.367;
Q3.wnp = pi / (Q3.Tp * sqrt(1 - Q2.zeta^2));
Q3.G = tf(Kdc * Q3.wnp^2, [1, 2 * Q2.zeta * Q3.wnp , Q3.wnp^2]); 
