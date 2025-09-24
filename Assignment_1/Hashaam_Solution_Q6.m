TotalImpedance = (1 / ZL1 + 1 / (ZR1 + ZC1));
VoltageAtJunction = TotalImpedance * Q2.G; 


% IOutRight = Q2.G * (ZL1) / ((1 / (ZR1 + ZC1) + 1 / ZL1)^-1);
CapacitorVoltage = VoltageAtJunction * ZC1/(ZR1 + ZC1); 

VOut = CapacitorVoltage * Q1.G % Gives us VOut

% Since I_in = 1, VOut = VOut / I_in
minrealoutput = minreal(VOut, 0.1)
p = pole(minrealoutput)
Q6.mdp = 109; 
Q6.dp = [109, 109, 455, 455, 588];
Q6.ndp = [31094];
a1Submit