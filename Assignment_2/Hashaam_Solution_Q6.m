% Author: Hashaam Zafar
% Date: 24/09/2025

Q6.wn1 = (4.44 * Q5.zeta - 1.15) / Q5.Tr1 - 5.5;
Q6.G = tf(225 * Q6.wn1^2, [1, 2 * Q5.zeta * Q6.wn1, Q6.wn1^2]);
step(Q6.G);
a2Submit