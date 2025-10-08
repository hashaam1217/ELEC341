% Author: Hashaam Zafar
% Date: 03/10/2025

% Creating 4 equations for each capacitor
% voltage. 
ZT1 = (1 / Zr20 + 1 / Zl20)^-1;

syms V0 V1 V2 V3; 
eq1 = V0 / Zl0 + V0 / Zc0 - 1 + (V0 - V2) / ZT1 == 0;
eq2 = V1 / Zl1 + V1 / Zc1 + (V1 - V3) / Zr31 + (V1 - V2) / Zr21 == 0; 
eq3 = (V2 - V3) / Zl32 + V2 / Zc2 + (V2 - V0) / ZT1 + (V2 - V1) / Zr21 == 0;
eq4 = V3 / Zc3 + (V3 - V1) / Zr31 + (V3 - V2) / Zl32 == 0; 

% sol = solve({eq1, eq2, eq3, eq4}, {V0, V1, V2, V3});
[V0, V1, V2, V3] = solve(eq1, eq2, eq3, eq4, V0, V1, V2, V3);
D1 = V1 / s; 
D3 = V3 / s; 
[num1, den1] = numden(D1/1);
[num3, den3] = numden(D3/1);

Q1.Gd1 = minreal(tf(sym2poly(num1), sym2poly(den1)));
Q1.Gd3 = minreal(tf(sym2poly(num3), sym2poly(den3)));

current1 = -V1 / Zl1; 
current32 = (V2 - V3) / Zl32; 

[num1, den1] = numden(current1/1);
[num3, den3] = numden(current32/1);

Q2.Gf1 = minreal(tf(sym2poly(num1), sym2poly(den1)));
Q2.Gf32 = minreal(tf(sym2poly(num3), sym2poly(den3)));

% Q3
ZBR = 1 / Br; 
ZJr = 1 / (s * Jr);
ZT2 = (1 / ZBR + 1 / ZJr)^-1;
syms iw w; 
Vs = ZT2 * Km *iw; 
eq1 = Vs - w == (Vs - w);

% Y11 = 1 / Zl0 + 1 / Zc0 + 1 / ZT1;
% Y12 = 0;
% Y13 = -1 / ZT1; 
% Y14 = 0; 
% Y21 = 0; 
% Y22 = 1 / Zl1 + 1 / Zc2 + 1 / Zr31 + 1 / Zr21;
% Y23 = - 1 / Zr21; 
% Y24 = - 1 / Zr31; 
% Y31 = 0; 
% Y32 = - 1 / ZT1; 
% Y33 = 1 / Zl32 + 1 / Zc2 + 1 / ZT1;
% Y34 = - 1 / Zl32; 
% Y41 = 0;
% Y42 = - 1 / Zr31; 
% Y43 = - 1 / Zl32; 
% Y44 = 1 / Zc3 + 1 / Zr31 + 1 / Zl32; 

% Y = [Y11 Y12 Y13 Y14;
%      Y21 Y22 Y23 Y24;
%      Y31 Y32 Y33 Y34;
%      Y41 Y42 Y43 Y44];
% I = [1;0;0;0];
% V = Y \ I


