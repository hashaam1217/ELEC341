% A step function is 1 / s
% A pulse would be a step - delayed step.
Q5.tau = round(-1 / Q4.p * 1000) / 1000;
SquarePulseTime = Q5.tau; 
% SquarePulse = (exp(-s*SquarePulseTime)) / s;

% SquarePulse = -Q4.p / s - (-Q4.p * exp(s * Q4.p))/s;

% square_pulse = (-1/Q5.tau)*(1/s)+(1/Q5.tau)*(exp(s/(Q4.p))/s);
t = [0:1:100];
t = t / 1000;

% impulse(DataSheetCurve)
% y = DataSheetCurve * SquarePulse;
% impulse(y, t)

% pulse = double(t >= 0 & t < SquarePulseTime) * (1 / SquarePulseTime);
pulse = double(t >= 0 & t < Q5.tau) * (1 / Q5.tau);  % 1x101 vector, area = 1

% Step 4: simulate response to that pulse
Q5.y = lsim(DataSheetCurveOriginal, pulse, t);
Q5.y = Q5.y(:)'; 
% Q5.y = impulse(SquarePulse * DataSheetCurve, t)
% Q5.y = Q5.y';
plot(t, y); 
grid on; 
xlabel("Time (ms)")
ylabel("Amplitude (V)")