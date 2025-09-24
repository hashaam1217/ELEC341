OmegaN = 93;
DataSheetCurve = tf(OmegaN^2, [1, -2 * sigma, OmegaN^2]) * Kdc;
% step(DataSheetCurve)
grid on;
complexroots = roots([1 (-2 * sigma) OmegaN^2]);
Q4.p = 10 * real(complexroots(1));
ApproximateCurve = 90 / (s + 90); % Explore later
DataSheetCurveOriginal = DataSheetCurve; 
DataSheetCurve = DataSheetCurve * ApproximateCurve ;
% a1DSPlot(10078020)
% Second part of the question 
t = [0:1:100];
t = t / 1000;
% y = impulse(ApproximateCurve, t)
% step(ApproximateCurve, t)
y = impulse(DataSheetCurve, t);
% step(DataSheetCurve, t)
% plot(t*1000, y');
grid on;
title('DataSheetCurveApprox');
xlabel('Time (ms)');
ylabel('Voltage(V)');

Q4.y = y';

