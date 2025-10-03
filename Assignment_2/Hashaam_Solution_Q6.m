% Author: Hashaam Zafar
% Date: 24/09/2025

Q6.wn1 = (4.44 * Q5.zeta - 1.15) / Q5.Tr1 ;
%TempZeta = (Q5.Tr1*(Q6.wn1 + 5.5) + 1.15) / 4.4 + 0.4 ;
TempZeta = Q5.zeta; 
Q6.G = tf(200 * Q6.wn1^2, [1, 2 * TempZeta * Q6.wn1, Q6.wn1^2]);
Q6.wn1 = (4.44 * 1.7 - 1.15) / Q5.Tr1 -5.5;
step(Q6.G);
a2Submit