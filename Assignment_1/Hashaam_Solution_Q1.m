% Isolating left and right parts of the question. I will find the transfer
% functions independently. 

% For the first part, we will treat the impedances left of the opamp as Z1
% and the impedances above the opamp as Z2. 

% Calculating for the left OpAmp. 
Z1 = (1 / ZC1 + 1 / ZR1)^(-1);
Z2 = (1 / ZC1 + 1 / ZR1 + 1 / ZL1)^(-1); 
TFLEFT = - Z2 / Z1;

% Calculating for the right OpAmp. Assuming Z3 and Z4 replacing Z1 and Z2's
% positions. 
Z3 = (1 / ZR2 + 1 / ZL2)^(-1); 
Z4 = (1 / ZR2 + 1 / ZC2)^(-1);

TFRIGHT = - Z4 / Z3;

TFFINAL = TFLEFT * TFRIGHT;
Q1.G = TFFINAL;


