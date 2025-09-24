% Combining the impedances for the top portion of the circuit (in series).
% Additionally combining the impedances in parallel at the bottom of the
% circuit. 

Z5 = ZR2 + ZL2 + ZC2; 
Z6 = (1 / ZL1 + 1 / (ZR1 + ZC1))^(-1);

MatrixA = [Z6 + ZC2 + ZL2, -ZL2; 
          -ZL2, Z5 + ZR2 + ZL2];
MatrixB = [-Z6; -ZR2];

TempSol = inv(MatrixA) * MatrixB
Q2.G = TempSol(1) + 1