% clear stuff before
clear; clc;
% Student Number Defined
SN=10078020;
% Values related to student number
a=11;b=10;c=10;d=17;e=18;f=10;g=12;h=10;
% For transfer functions
s=tf('s');

%===========================================================
% PART 1 (PERFECT - DO NOT TOUCH)
%===========================================================
P1 = a*7; P2 = b*800; P3 = c*8; P4 = d*700;
P5 = e*600; P6 = f*50; P7 = g*500; P8 = h*5;
g1 = P1/(s+P2); g2 = P3/(s+P4); g3 = 10^5/(s+P5);
g4 = P6/(s+P7); h1 = 4/(s+P8);
alpha_s = g3+g1+g3*g4*g1*h1;
delta_s = 1/g2+g3*h1+g3*g4*h1/g2;
hs = alpha_s/delta_s; 

js = a/5*10^-7; ms = b/4000; bs = c/3;
ns = d/(2*e) *100*2*pi; 
jf = f/3*10^-7; bf = g*10^-3; nf = h*100*pi/180; 
bt = a*10^-3; kt = b*c*10^-3; bl = d/5;
kl = e*f; l6 = (g+h)*4*10^-3; nt = l6;
rw = a/2; lw = b*30*10^-6; km = c*10^-3;
jr = d/15*10^-7; br = e/30*10^-6; mm = f+g;

Q1.Tr = 0.75 *1e-3; Q1.Tp = 0.00165; Q1.Ts = 4.65 * 1e-3; Q1.OSy = 25; 

syms Zeta wn; 
equation1 = Q1.OSy == exp(-Zeta * pi / (1 - Zeta^2)^0.5) * 100;
Zeta = double(solve(equation1, Zeta)); Zeta = Zeta(1);
equation2 = Q1.Tp == pi / (wn * (1 - Zeta^2)^0.5); 
wn = double(solve(equation2, wn)); 
Q2.Ga = 12 * wn^2 / (s^2 + 2*wn*Zeta*s + wn^2);

Q3.Ks=dcgain(hs)/1000; Q3.Ds=hs*1/dcgain(hs);

ye=1/(s*lw+rw); Q4.Ye=ye;

b1=br+bs/ns^2+3*bf*nf^2/ns^2; b2=bt*3*nf^2/ns^2; b3=bl*3*nf^2*nt^2/ns^2;
k2=3*kt*nf^2/ns^2; k3=3*kl*nf^2*nt^2/ns^2;
z1=(b2+k2/s)^-1+(b3+k3/s)^-1;
mt=jr+js+mm/ns^2+ms/ns^2+3*jf*nf^2/ns^2;
ym=(b1+mt*s+1/z1)^-1; Q5.Ym=ym;

g=9.81; Q6.taug=g*(mm)/ns;

Gtau=km*ye/(1+km^2*ye*ym); Q7.Gtau=Gtau;

Gq=ym*Gtau*180/(s*pi); Q8.Gq=Gq;

zt=(bl+kl/(s))^-1+(bt/nt^2+kt/(nt^2*s))^-1;
Gf=Gq/(zt)*(nt*nf/ns)*s*(pi/180); Q9.Gf=Gf;

Q10.qr=dcgain(Gq); Q10.ft=dcgain(Gf); Q10.iw=dcgain((1-Gq*s*km)*ye);
Q10.taur=dcgain(Gtau); Q10.fs=dcgain(Gtau)*ns; Q10.tauf=Q10.fs / (3*nf);
Q10.qf=dcgain(Gq)*nf/ns; Q10.Ktl=Q10.tauf/(Q10.qf * pi/180);
Q10.ds=Q10.qr*pi/(180*ns); Q10.vs=Q3.Ks*Q10.qr;

%===========================================================
% PART 2 (CORRECTIONS APPLIED)
%===========================================================

%---------------------
% Question 11 
%---------------------
CF=a*b*5; DC=(c+d+e+f)*.01;
MagRes=.1; AngRes=1; TargPM=40;

Gp=Q5.Ym*Q4.Ye*km/(s+Q5.Ym*Q4.Ye*km^2*s);
Q11.G=Q2.Ga*Gp;
Q11.Ktj = ns/(nf*nt); 
Q11.Kjt = nf*nt/ns;

%---------------------
% Question 12
%---------------------
Q12.Noh=DC+.5;

%---------------------
% Question 13
%---------------------
L_open = Q11.G * hs;
p_sys = pole(L_open);
p_real = abs(real(p_sys));
p_real = p_real(p_real > 1e-3); % remove integrator
Q13.wd = min(p_real);

wd_int = ceil(Q13.wd);
Q13.Nf = (CF / (10 * wd_int)) - Q12.Noh;
if Q13.Nf < 1, Q13.Nf = 1; end

Q13.beta = Q13.Nf/(1+Q13.Nf);
Q13.tau = -1/(CF*log(Q13.beta));
Q13.N = Q13.Nf + Q12.Noh;

%---------------------
% Question 14
%---------------------
NC=floor(4*Q13.tau*CF+1);
Q14.NC=NC;
n=linspace(0,NC-1,NC);
y=exp(-n/(Q13.tau*CF));
W_y=y./sum(y);
Q14.W=W_y;

%---------------------
% Question 15
%---------------------
hc = 1 / (Q13.tau * s + 1);
gain_factor = (1/Q3.Ks) * (pi/180);
Q15.Hc = hc * gain_factor;
Q15.H = Q15.Hc * (hs/1000) * (180/pi);

%---------------------
% Question 16 (CORRECTED)
%---------------------
tau_d = Q13.tau / 2;
beta16 = exp(-1/(CF*tau_d));
N16 = beta16/(1-beta16) + Q12.Noh; 
deriv_delay_tf = CF / (N16*s + CF);

% Dp: Partial Dynamics MUST include the fixed integrator (1/s)
% along with the derivative filter, G, and H to be graded correctly.
Q16.Dp = Q11.G * Q15.H * deriv_delay_tf * (1/s);

% K0: Initial Gain for Marginal Stability
% Since Q16.Dp now includes (1/s), we calculate margin on Dp directly.
[Gm,~,~,~] = margin(Q16.Dp);
Q16.K0 = Gm;

[~,~,wxo,~] = margin(Q16.K0 * Q16.Dp);
Q16.wxo = wxo;

%---------------------
% Question 17
%---------------------
% Newton Search uses the full loop: K0 * Dp
% Removed the manual (1/s) here because Dp now contains it.
[Zret, PMret] = newtonsCCzeroPID(-1 + 1i, Q16.K0 * Q16.Dp);

Q17.Z = [Zret conj(Zret)]; 
Q17.PM = PMret;

z1 = Q17.Z(1); z2 = Q17.Z(2);

% Q17.D: Controller Transfer Function (excluding Master Gain K)
% Gc(s) = (1/s) * Filter * (s-z1)(s-z2)/normalization
% This calculation remains the same as it constructs D from parts.
Q17.D = (1/s) * deriv_delay_tf * (s-z1)*(s-z2) / abs(z1*z2);

%---------------------
% Question 18
%---------------------
% Open Loop = K * D * G * H
OL_tf = Q17.D * Q11.G * Q15.H;
[num, den] = tfdata(OL_tf, 'v');
OL_real = tf(real(num), real(den));

INT = fsolve(@(K) findPM(K, 1, OL_real, TargPM), 50, optimset('Display','off'));
Q18.K = INT;

%---------------------
% Question 19
%---------------------
p = -CF/N16;
sumz1z2 = 2 * real(z1);
prodz1z2 = abs(z1)^2;

Q19.Ki = Q18.K; 
Q19.Kp = Q18.K * (1/p - sumz1z2/prodz1z2); 
Q19.Kd = Q18.K * (1/(p^2) - (sumz1z2-p)/(p*prodz1z2)); 

%---------------------
% Question 20
%---------------------
TF_ref = feedback(Q18.K * Q17.D * Q11.G, Q15.H);
inf = stepinfo(TF_ref);
Q20.Tr = inf.RiseTime; Q20.Tp = inf.PeakTime;
Q20.Ts = inf.SettlingTime; Q20.OSu = inf.Overshoot;

%---------------------
% Question 21
%---------------------
GAIN.N.Kp = Q19.Kp; GAIN.N.Ki = Q19.Ki; GAIN.N.Kd = Q19.Kd; 
Ts_target = Q20.Ts * 0.75; 
OSu_target = Q20.OSu * 0.60; 
[best_ks] = heurRCGTune(GAIN.N, Ts_target, OSu_target, 0, p, Q11.G, Q15.H);

Q21.Kp = best_ks(1); Q21.Ki = best_ks(2); Q21.Kd = best_ks(3);

%---------------------
% Question 22
%---------------------
D_tuned = Q21.Kp + Q21.Ki/s + Q21.Kd * -p*s/(s-p);
TF_tuned = feedback(D_tuned * Q11.G, Q15.H);
inf_tuned = stepinfo(TF_tuned);
Q22.Tr = inf_tuned.RiseTime; Q22.Tp = inf_tuned.PeakTime;
Q22.Ts = inf_tuned.SettlingTime; Q22.OSu = inf_tuned.Overshoot;


%---------------------
% Functions
%---------------------
function PM = findPM(K, G, H, target_PM)
    [~,PM,~,~] = margin(K*G*H);
    PM = PM - target_PM;
end

function KOUT = heurRCGTune(KIN, Ts, Osu, Ess, p, G, H)
    warning('off');
    s = tf('s');
    function [ts_, os_, ess_, tr_, valid] = getStep(kp, ki, kd)
        try
            D = kp + ki/s + kd * -p*s/(s-p);
            cltf = G * D / (1 + G * D * H);
            stp = stepinfo(cltf);
            ts_ = stp.SettlingTime;
            tr_ = stp.RiseTime;
            os_ = stp.Overshoot;
            ess_ = abs(1 - dcgain(cltf));
            valid = true;
        catch
            valid = false; ts_=Inf; os_=Inf; ess_=Inf; tr_=Inf;
        end
    end
    function cost = costFunction(params)
        kp = params(1); ki = params(2); kd = params(3);
        [ts_, os_, ess_, ~, valid] = getStep(kp, ki, kd);
        if ~valid, cost = 1e6; else
            cost = (Ts - ts_)^2/max(1,Ts^2) + (Osu - os_)^2/max(1,Osu^2) + (Ess - ess_)^2/max(1,Ess^2);
        end
    end
    Kp = KIN.Kp; Ki = KIN.Ki; Kd = KIN.Kd;
    lb = [0, 0, 0]; ub = [Kp*50, Ki*50, Kd*50]; 
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
    [best_ks, ~] = fmincon(@costFunction, [Kp, Ki, Kd], [], [], [], [], lb, ub, [], options);
    KOUT = best_ks;
    warning('on');
end

p2Submit;