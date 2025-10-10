% Date: 08/10/2025
% Name: Hashaam Zafar
% Student Number: 10078020
ZBr = 1 / Br; 
ZJr = 1 / (s * Jr); 
ZT1 = (1 / ZBr + 1 / ZJr)^-1;

A = [1 Ke/(s*Lw+Rw);
    -Km Jr*s+Br];
B = [1/(s*Lw+Rw); 0];
Temp_Ans = A\B;
%Q3.Gi = Temp_Ans(1);
%Q4.Gw = Temp_Ans(2);


Q3.Ye = Temp_Ans(1) / (Temp_Ans(2) * Ke);
Q3.Ye = 1 / (s * Lw +Rw);
Q3.Ym = 1 / (s*Jr + Br + n^2*(s*Jg +s*Jf+Bg+Bf));
