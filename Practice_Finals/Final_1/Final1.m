% Practice Final 1
% Author: Hashaam Zafar
% Date: 19/12/2025

%% Setup
SN = 10078020; 
A = 11; 
B = 10;
C = 10; 
D = 17; 
E = 18; 
F = 10; 
G = 12; 
H = 10; 

s = tf('s');

% Table 1
Dp = A * 1e-2; % (stored in m)  <-- 11 cm = 0.11 m
Jp = B * 2e-4; % Nms^2/rad
Mt = C / 2; % Kg
Bp = D * 1e-5; % Nms/rad
Bw = E + F;  % Ns/m
Kt = G / 2; % N/m

% Table 2
K0 = 3e6; % V/rad
P1 = B * 25; % rad/s
P2 = C * 25; % rad/s

% Table 3
Rw = A / 3; % Ohm
Lw = B / 5; % H
Km = C * 2e-2; % Nm/A
Jm = D * 4e-4; % Nms^2/rad
Bm = E * 2e-6; % Nms/rad

% Table 4  (IMPORTANT: CF is in Hz; DC is in %)
CF = D * 20;          % Hz
DCpct = E + F;        % %
DC = DCpct / 100;     % fraction

% Table 5
WnRes = 0.1; % rad/s
ZetaRes = 10^-2; % pure
TargPM = 60 * pi / 180; % deg
OSu = 70; 
Ts = 50;
Tr = Ts; 
Ess = 0; 

H1 = 1 / (s + P1);
H2 = 1 / (s + P2);

%% Question 1
Peak = 5.53658;
Steady_State = 4.62361;
Overshoot = (Peak - Steady_State) / Steady_State;

syms Zeta;
equation1 = Overshoot == exp(-Zeta * pi / (1 - Zeta^2)^0.5);
Zeta = solve(equation1, Zeta, 'Real', true);
Zeta = double(Zeta);

Ts = 0.102;
Wn = 4 / (Zeta*Ts);

Q1.Ga = 4.623 * Wn^2 / (s^2 + 2 * Zeta * Wn * s + Wn^2);

%% Question 2
Hs = K0*H1*H2/(1+H1+H2);
Q2.Ks = dcgain(Hs);
Q2.Ds = Hs / dcgain(Hs);

%% Question 3
r = Dp/2;
Q3.Jj = Jm + Jp + Mt*r^2;
Q3.Bj = Bm + Bp + Bw*r^2;
Q3.Kj = Kt*r^2;

%% ---------- Q4: Plant Model ----------
Q4.Ye = 1 / (Lw*s + Rw);
Q4.Ym = s / (Q3.Jj*s^2 + Q3.Bj*s + Q3.Kj);

n = 1;
nKm = n*Km;

omega_over_x3 = minreal( (Q4.Ym*nKm*Q4.Ye) / (1 + Q4.Ym*(nKm^2)*Q4.Ye) );
Q4.Gp = minreal( omega_over_x3 * (1/s) );

%% ---------- Q5: State Space ----------
Jt = Q3.Jj;
Bt = Q3.Bj;
Kj = Q3.Kj;

Q5.A = [ -(Rw/Lw),    -(Km/Lw),     0;
          (Km/Jt),    -(Bt/Jt),   -(1/Jt);
            0,          Kj,         0 ];

Q5.B = [ 1/Lw; 0; 0 ];

Q5.C = [ 0,      (Bw*r),   0;
         0,        0,    (1/r) ];

Q5.D = [0; 0];

%% ---------- Q6: IIR Filter ----------
Q6.GHs = minreal(Q1.Ga * Q4.Gp * Hs);

p = pole(Q6.GHs);
p = p(real(p) < 0);
[~, idx] = max(real(p));
pdom = p(idx);

wd = abs(imag(pdom));
if wd < 1e-12
    wd = abs(real(pdom));
end
Q6.wd = wd;

% ---- FIXED Q6 (uses CF in Hz and DC in %) ----
Tsamp = 1/CF;

wd_up = ceil(Q6.wd);                        % ROUND-UP wd
Q6.Nf = floor((DC * CF + 0.5) / wd_up);           % max delay multiple (integer)
Q6.tau = round_sig(Q6.Nf / CF, 2);          % tau in seconds, round 2 sig digs
Q6.beta = exp(-Tsamp / Q6.tau);             % weighting factor
% ---- END FIX ----

%% ---------- Q7: FIR Filter ----------
Q7.num = 2*Q6.Nf + 1;

%% ---------- Q8: Feedback Path ----------
Q8.N = Q6.Nf;

Fiir = tf(1, [Q6.tau 1], 'InputDelay', Q8.N * Tsamp);
Q8.Hc = minreal( (Fiir / Q2.Ks) );

%% ---------- Q9: System Model (Gc = 1) ----------
Kjt = r;
Ktj = 1/r;
Q9.Kjt = Kjt;
Q9.Ktj = Ktj;

Q9.G  = minreal(Q1.Ga * Q4.Gp);
Q9.H  = minreal(Q8.Hc * Hs);
Q9.GH = minreal(Q9.G * Q9.H);

%% ---------- Q10: Partial Dynamics (FAST FRF VERSION) ----------
Q10.Dp = minreal( Q9.GH / (s*(Q6.tau*s + 1)) );

wPM = logspace(-3, 6, 6000);
[Q10.K0, Q10.wxo] = gain_margin_from_frf(Q10.Dp, wPM);

%% ----------------- Helpers -----------------
function y = round_sig(x, nSig)
    if x == 0
        y = 0;
        return;
    end
    p = floor(log10(abs(x)));
    y = round(x, nSig - 1 - p);
end

function [K0, wxo] = gain_margin_from_frf(Lsys, w)
    L = squeeze(freqresp(Lsys, w));
    mag = abs(L);
    ph  = unwrap(angle(L)) * 180/pi;

    idx = find((ph(1:end-1) > -180) & (ph(2:end) <= -180), 1, 'first');
    if isempty(idx)
        K0 = NaN; wxo = NaN; return;
    end

    p1 = ph(idx); p2 = ph(idx+1);
    w1 = w(idx);  w2 = w(idx+1);
    t = (-180 - p1) / (p2 - p1);
    wxo = 10^(log10(w1) + t*(log10(w2) - log10(w1)));

    m1 = mag(idx); m2 = mag(idx+1);
    mag_at = 10^(log10(m1) + t*(log10(m2) - log10(m1)));

    K0 = 1 / mag_at;
end
