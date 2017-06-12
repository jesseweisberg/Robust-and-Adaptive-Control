% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Purpose:
%   This file sets up pitch matrices and 
%   and designs a LQGLTR state feedback controller 
%   with an accel cmd (essentially LQR with loop-transfer recovery)
%
% Author: Jesse Weisberg
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% RSLQR Design 
Za_V = -1.3046;
Ma = 47.711;
Zd_V = -.2142;
Md = -104.83;
V = 886.78; % (fps)
Za = V*Za_V;
Zd = V*Zd_V;
w_act = 2.*pi*11; % actuator natural frequency rps
z_act = 0.707; % actuator damping

%Define plant matrices without actuator
Ap1 = [Za_V, 1; Ma, 0];
Bp1 = [Zd_V; Md];
Cp1 = [Za, 0; 0, 1];
Dp1 = [Zd; 0];

%Define plant matrices with actuator
Ap2 = [Za_V, 1, Zd_V, 0; Ma, 0, Md, 0; 0, 0, 0, 1; 0, 0, -w_act*w_act, -2*z_act*w_act]
Bp2 = [0; 0; 0; w_act*w_act]
Cp2 = [Za, 0, Zd, 0; eye(4,4)]
Dp2 = 0.*Cp2*Bp2

%Define wiggle matrices
o1 = [0; 0];
Aw = [0, Cp1(1,:); o1, Ap1];
Bw = [Dp1(1); Bp1];

%taken from HW3 RSLQR with abs...rdmin < .0005
Kc = [0.0093   -2.5874   -0.2360];

%% Optimal Observer Design
Qo = diag([ .001 .0014 .005]);
Ro = 1000*diag([.025^2 .001^2]);
Co = .00001;

p1 = 10^10; p2 = 10^4; p3 = 10^2;
v1 = 100; v2 = .3162; v3 = .001;

Cw = [1 0 0; 0 0 1];
Dw = [0; 0];
sys_1 = ss(Aw, Bw, Cw, Dw);
zz = [-10];
sys_2 = square_system_inputs_pp(sys_1,zz);
disp(['Bbar MAtrix z = ' num2str(zz) ])
sys_2.b
tzero(sys_2)

%Qf = Qo + (1/p1)*Bw*(Bw');
Qv = Qo + Co*((v1+1)/v1)*(sys_2.b*sys_2.b');
Rv = (v1/(v1+1))*Ro;

[Kf,~,~] = lqr(Aw', Cw', Qv, Rv);
Kf = Kf';


%% Controller
Ac = [ zeros(1,4); Kf(:,1) Aw-Bw*Kc-Kf*Cw]
Bc1 = [1 0 0 0 0; zeros(3,1) zeros(3,1) Kf(:,2) zeros(3,1) zeros(3,1)]
Bc2 = [-1 -1 0 0]'
Cc = [0 -Kc]
Dc1 = zeros(1,5)
Dc2 = zeros(1,4)

%% Forming Closed Loop System (RSLQR and Estimator)
Z = inv(eye(size(Dc1*Dp2))-Dc1*Dp2);
Acl = [ (Ap2+Bp2*Z*Dc1*Cp2) (Bp2*Z*Cc);
    (Bc1*(Cp2+Dp2*Z*Dc1*Cp2)) (Ac+Bc1*Dp2*Z*Cc)];
Bcl = [ Bp2*Z*Dc2;
    (Bc2+Bc1*Dp2*Z*Dc2)];
Bcl = Bcl(:,1);
Ccl = [(Cp2+Dp2*Z*Dc1*Cp2) (Dp2*Z*Cc)];
Dcl =(Dp2*Z*Dc2);
Dcl = Dcl(:,1);
sys_cl = ss(Acl,Bcl,Ccl,Dcl);

Bv = [ Bp2*Z*Dc1;
    (Bc1+Bc1*Dp2*Z*Dc1)];
Cv = [Z*Dc1*Cp2 Z*Cc];
Dv = [Z*Dc1];
sys_clv = ss(Acl,Bv,Cv,Dv);
sys_clvr = ss(Acl, Bv, Cv*Acl,Cv*Bv);

%% Closed Loop System of just RSLQR
% Define Controller Matrices
Acr = 0.;
Bc1r = [1. 0. 0. 0. 0.];
Bc2r = -1;
Ccr = -Kc(1);
Dc1r = [0. -Kc(2:3) 0. 0.];
Dc2r = 0.;

% Form the closed loop system of RSLQR designed controller with full plant (w/ actuator model)  
Zr = inv(eye(size(Dc1r*Dp2))-Dc1r*Dp2);
Aclr = [ (Ap2+Bp2*Zr*Dc1r*Cp2) (Bp2*Zr*Ccr);
    (Bc1r*(Cp2+Dp2*Zr*Dc1r*Cp2)) (Acr+Bc1r*Dp2*Zr*Ccr)];
Bclr = [ Bp2*Zr*Dc2r;
    (Bc2r+Bc1r*Dp2*Zr*Dc2r)];
Cclr = [(Cp2+Dp2*Zr*Dc1r*Cp2) (Dp2*Zr*Ccr)];
Dclr = (Dp2*Zr*Dc2r);
sys_clr = ss(Aclr,Bclr,Cclr,Dclr);

%% State-Space model of loop gain at the plant input Lu
qq = logspace(-3,2, 500); %these are the varying values of q11, evenly spaced on log scale
w = logspace(-3, 3, 500); %varying frequency values
t = linspace(0, 2, 1000); %time scale  (t = 0:.002:2;)
rtd = 180/pi; %radian to degree scale
% for plotting unit circle in nyquist later
dd = 0.:.001:2*pi;
xx1 = cos(dd)-1; 
yy1 = sin(dd);

A_Lur = [ Ap2 0.*Bp2*Ccr; Bc1r*Cp2 Acr];
B_Lur = [ Bp2; Bc1r*Dp2];
C_Lur = -[ Dc1r*Cp2 Ccr]; %change sign for loop gain
D_Lur = -[ Dc1r*Dp2];
sys_Lur = ss(A_Lur,B_Lur,C_Lur,D_Lur);
magdbr = 20*log10(abs(squeeze(freqresp(sys_Lur,w))));
wcr = crosst(magdbr,w); % LGCF, assumes Lu is a scalar

A_Lu = [ Ap2 0.*Bp2*Cc; Bc1*Cp2 Ac];
B_Lu = [ Bp2; Bc1*Dp2];
C_Lu = -[ Dc1*Cp2 Cc];%change sign for loop gain
D_Lu = -[ Dc1*Dp2];
sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);
magdb = 20*log10(abs(squeeze(freqresp(sys_Lu,w))));
wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar

%singular values LQG/LTR
sr = sigma(sys_Lu,w,3); % singular values of stability robustness
srmin = min(abs(sr));   % minimium singular value of stability robustness
rd = sigma(sys_Lu,w,2); % singular values of return difference
rdmin = min(abs(rd));   % minimum singular value of return difference

%singular values RSLQR
srr = sigma(sys_Lur,w,3); % singular values of stability robustness
srminr = min(abs(srr));   % minimium singular value of stability robustness
rdr = sigma(sys_Lur,w,2); % singular values of return difference
rdminr = min(abs(rdr));   % minimum singular value of return difference

%% Frequency Domain Analysis at the Plant Input
Lu = freqresp(sys_Lu,w);
T  = freqresp(sys_cl,w);
S = 1 - T;

Lur = freqresp(sys_Lur,w);
Tr  = freqresp(sys_clr,w);
Sr = 1 - Tr;

% i=1;
% SNRr = zeros(numel(w),5);
% while i<numel(w)
% SNRr(i,:) = (Cv*inv(j*i*eye(size(Acl))-Acl)*Bv+Dv);
% SNRr(i,1) = 1/SNRr(i,1);
% i=i+1;
% end
Luz = freqresp(sys_clv,w);
figure('Name', 'Accel Noise to Control'),
semilogx(w,20*log10(abs(squeeze(Luz(1,1,:)))), 'k');

Luzr = freqresp(sys_clvr,w);
figure('Name', 'Accel Noise to Control Rate'),
semilogx(w,20*log10(abs(squeeze(Luzr))), 'k');


figure('Name','Nyquist Plot at Plant Input'),
plot(xx1,yy1,'r',real(squeeze(Lu)),imag(squeeze(Lu)),'k', real(squeeze(Lur)),imag(squeeze(Lur)), 'b', 'LineWidth',2);grid
axis([-2 2 -2 2]);
legend('Unit Circle','LQG/LTR','RSLQR','Location','Best');
xlabel('Re(L)')
ylabel('Im(L)')
title('Nyquist Plot at Plant Input')

figure('Name','Bode Magnitude at Plant Input'),
semilogx(w,20*log10(abs(squeeze(Lu))),'k', w, 20*log10(abs(squeeze(Lur))), 'b', 'LineWidth',2);grid
legend('LQG/LTR','RSLQR','Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag')
title('Bode at Input')

figure('Name', 'Bode from Matlab'),
margin(sys_Lu)
hold on 
margin(sys_Lur)
hold off

figure('Name','Return Difference at Plant Input'),
semilogx(w,20*log10(abs(rd)),'k', w,20*log10(abs(rdr)),'b', 'LineWidth',2);grid
%legend([' min(I+Lu) = ' num2str(rdmin), '],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input')

figure('Name','Stability Robustness at Plant Input'),
semilogx(w,20*log10(abs(sr)),'k', w,20*log10(abs(srr)),'b','LineWidth',2);grid
%legend([' min(I+invLu) = ' num2str(srmin)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')


% Loop Gain/Phase Crossover Frequency  
disp('Classical Margins')
allmargin(sys_Lu) 

%Singular Value Gain/Phase Margins
disp('\\\ LQG/LTR Singular Value Margins \\\')
RDu_nGM = 1/(1+rdmin);
RDu_pGM = 1/(1-rdmin);
RDu_Pha = 2*asin(rdmin/2);
RDu_nGM_dB = 20*log10(RDu_nGM);
RDu_pGM_dB = 20*log10(RDu_pGM);
RDu_Pha_deg = 180*RDu_Pha/pi ;
disp('Return Difference at Input')
disp([ '  Negative Gain Margin: ' num2str(RDu_nGM_dB) ' dB,  Positive Gain Margin: ' num2str(RDu_pGM_dB) ' dB'])
disp([ '  Phase Margin: +/-' num2str(RDu_Pha_deg)  ' deg'])
disp(' ')
SRu_nGM = 1-srmin;
SRu_pGM = 1+srmin;
SRu_Pha = 2*asin(srmin/2);
SRu_nGM_dB = 20*log10(SRu_nGM);
SRu_pGM_dB = 20*log10(SRu_pGM);
SRu_Pha_deg = 180*SRu_Pha/pi ;
disp('Stability Robustness at Input')
disp([ '  Negative Gain Margin: ' num2str(SRu_nGM_dB) ' dB,  Positive Gain Margin: ' num2str(SRu_pGM_dB) ' dB'])
disp([ '  Phase Margin: +/-' num2str(SRu_Pha_deg)  ' deg'])
disp(' ')
    
%% Frequency Domain Analysis at the Plant Output (Sensitivity, Complementary Sensitivity)

T  = freqresp(sys_cl,w); % Complementary Sensitivity
S = 1 - T; % Sensitivity
T_st = 20*log10(abs(squeeze(T(1,1,:))));
S_st = 20*log10(abs(squeeze(S(1,1,:))));
Tmax = max(T_st); % Inf Norm of T in dB
Smax = max(S_st); % Inf Norm of S in dB
    
figure('Name','Comp Sens T');
semilogx(w,T_st,'b');grid
legend(['Tmax = ' num2str(Tmax) ' (dB)'],'Location','Best');
title('Comp Sens T');
xlabel('Freq (rps)');ylabel('Mag (dB)');  

figure('Name','Sens S');
semilogx(w,S_st,'b');grid
legend(['Smax = ' num2str(Smax) ' (dB)'],'Location','Best');
title('Sens S');
xlabel('Freq (rps)');ylabel('Mag (dB)');

%% RSLQR and LQG/LTR Time Domain Plots for Az, a, q, del, deldot
LQGLTRstep = step(sys_cl, t);
RSLQRstep = step(sys_clr,t); %step response of tuned RSLQR
Azr = RSLQRstep(:,1); Az = LQGLTRstep(:,1);
ar = RSLQRstep(:,2); a = LQGLTRstep(:,2);
qr = RSLQRstep(:,3); q = LQGLTRstep(:,3);
delr = RSLQRstep(:,4); del = LQGLTRstep(:,4);
deldotr = RSLQRstep(:,5); deldot = LQGLTRstep(:,5);

figure('Name','Az (Step Response)')
plot(t, Azr, 'b', t, Az, 'r');
legend('RSLQR', 'LQG/LTR');
grid on;
xlabel('time (s)');
ylabel('Az (f/s^2)');
title('Az (Step Response)'); 
hold on

figure('Name','AOA (Step Response)')
plot(t, ar, 'b', t, a, 'r');
legend('RSLQR', 'LQG/LTR');
grid on;
xlabel('time');
ylabel('AOA (radians)');
title('AOA (Step Response)');
hold on


figure('Name','Pitch (Step Response)')
plot(t, qr, 'b', t, q, 'r');
legend('RSLQR', 'LQG/LTR');
grid on;
xlabel('time');
ylabel('Pitch (radians)');
title('Pitch (Step Response)');
hold on

figure('Name','Elevon Deflection (Step Response)')
plot(t, delr, 'b', t, del, 'r');
legend('RSLQR', 'LQG/LTR');
grid on;
xlabel('time (s)');
ylabel('Elevon Deflection (radians)');
title('Elevon Deflection (Step Response)');
hold on

figure('Name','Elevon Deflection Rate (Step Response)')
plot(t, deldotr, 'b', t, deldot, 'r');
legend('RSLQR', 'LQG/LTR');
grid on;
xlabel('time (s)');
ylabel('Elevon Deflection Rate (radians/s)');
title('Elevon Deflection Rate (Step Response)');
hold on












