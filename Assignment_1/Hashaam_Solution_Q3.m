Kdc = 50;
sigma = -9;

%t = 0:0.01:0.2;  
%y = Kdc * (1 + exp(sigma*t));
sys = tf(Kdc,[1 -sigma]) + Kdc / s;
% test = Kdc / (s - sigma); 
impulse(sys)
% step(sys)

Q3.Kdc = Kdc;
Q3.sigma = sigma;

%figure(1); clf; hold on; grid on;
%plot(y, t, 'k-',  'Linewidth', 3);
%title('Step Response (custom)');
%xlabel('Time (s)');
%ylabel('Voltage (V)');

