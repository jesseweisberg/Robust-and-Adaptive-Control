% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Purpose:
%   This file sets up the pitch matrices 
%   and designs a Hinf SF Accel Cmd Tracking autopilot
%
% Created : 2/8/2017, Kevin A Wise
%
% Modified:
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clc
clear all
close all
format short e

disp('****************************** Program Start ****************');
plot_file_name = 'Homework_4_HINF_2017.ppt';
save_plots = 0; % Flag to bypass saving plots
w = logspace(-3,4,1000);
t = linspace(0,2.5,1000);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);
rtd = 180/pi;

g    = 32.174;

m2ft = 3.2808;    % meters to feet conversion
ft2m = 1/3.2808;  % feet to meters conversion

% Public release airframe parameters
%*************************************************************************
% Plant Model
%*************************************************************************
% State Names 
%  ----------- 
%  AZ fps2         
%  q rps           
%  Dele rad        
%  Dele dot rps    
 
%  Input Names
%  -----------
%  Dele cmd deg    
% Plant.Ap = [ -0.576007     -3255.07            4.88557       9.25796;
%           -0.0410072       -0.488642       -2.03681       0      ;
%                  0                0               0             1      ;
%                  0                0           -8882.64       -133.266  ];
% Plant.Bp = [      0  ;  0   ;     0   ;  8882.64];
% Plant.Cp = [1 0 0 0];
% Plant.Dp = 0.*Plant.Cp*Plant.Bp;
% pstatnam = {'Az fps2' 'q rps' 'dele rad' 'deledot rps'};
% poutnam  = {'Az fps2'};
% pinnam   = {'delecmd rad'};
% plant=ss(Plant.Ap,Plant.Bp,Plant.Cp,Plant.Dp,'statename',pstatnam,'inputname',pinnam,'outputname',poutnam);

%**************************************************************************
% Aircraft Model Data
%**************************************************************************
% Model (taken from example 5.2, the aircraft pitch axis plant data)
% Aero constants used in the aircraft model:
Za_V = -1.3046;
Ma   = 47.711; % Positve Ma = unstable aircraft
Zd_V = -.2142;
Md   = -104.83;
V    = 886.78; % (fps)
Za   = V*Za_V;
Zd   = V*Zd_V;
w_act = 2.*pi*11; % actuator natural frequency rps
z_act = 0.707;    % actuator damping

grav = 32.174; % (fps2)

%**************************************************************************
% Plant model for analysis with actuator 
%**************************************************************************
% States are AOA (rad),  pitch rate (rps), dele (rad), deledot (rps)
disp('Plant Model')
disp('xpdot = Ap*xp + Bp*u')
disp('    y = Cp*xp + Dp*u')
disp('   xp = AOA (rad), pitch rate q (rps), dele (rad), deledot (rps)')
disp('    u =  delec (rad)')
disp('    y =  Az (fps2), AOA (rad), pitch rate q (rps), dele (rad), deledot (rps)')

disp(' ')

Plant.Ap = [Za_V  1.      Zd_V          0.; 
      Ma    0.       Md           0.;
      0.    0.       0.           1.;
      0.    0. -w_act*w_act -2*z_act*w_act];

% Input is dele (rad)
Plant.Bp = [0.; 0.; 0.; w_act*w_act ];

% Outputs are Az (fps2), AOA (rad), pitch rate (rps), dele (rad), deledot (rps)
Plant.Cp = [  Za   0.  Zd  0. ];
Plant.Dp = 0.*Plant.Cp*Plant.Bp;




% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~  Hinf Optimal Control  ~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Target Loop Gain Cross-Over Frequency
% This parameter set the target bandwidth

wlgcf = 1.2; %Hz
disp(' ')
disp(['Target lgcf = ' num2str(wlgcf)])
disp(' ')

%% Sensitivity Ws

Ws.Tau = 1/(2*pi*wlgcf);
Ws.K = 0.5/Ws.Tau;
%Ws.K = 0.8/Ws.Tau;
Ws.num = Ws.K * [Ws.Tau, 1];
Ws.den = [1, 0];
[Ws.A,Ws.B,Ws.C,Ws.D] = tf2ss(Ws.num,Ws.den);

Ws.sys = ss(Ws.A,Ws.B,Ws.C,Ws.D)
% Ws.A
% Ws.B
% Ws.C
% Ws.D
figure(1)
bodemag(Ws.sys)
grid on
hold on

%% Complementary Sensitivity Wt (what are tauN and TauD?)

Wt.TauN = 1/(3.5*2*pi*wlgcf);
Wt.TauD = 0.005;
Wt.K = 0.707;
Wt.num = Wt.K * [Wt.TauN, 1];
Wt.den = [Wt.TauD, 1];
[Wt.A,Wt.B,Wt.C,Wt.D] = tf2ss(Wt.num,Wt.den);

Wt.sys = ss(Wt.A,Wt.B,Wt.C,Wt.D);
% Wt.A
% Wt.B
% Wt.C
% Wt.D
bodemag(Wt.sys)
grid on
hold on

%% Control Activity (what is this)

Wc = 0.1;
Cc = [0, 0, 0, 1];
Wc_sys = ss(0,0,0,Wc);

bodemag(Wc_sys)
xlim([0.1,100])
if(save_plots == 1) saveppt2(plot_file_name); end   
hold off

%% Hinf State-Feedback Matrix Setup (with E,21 and D1,31 fixes)

HinfSF.A = [     Plant.Ap, zeros(4,1), zeros(4,1);...
            Ws.B*Plant.Cp,       Ws.A,          0;...
            Wt.B*Plant.Cp,          0,       Wt.A];

HinfSF.B = [     Plant.Bp;...
            Ws.B*Plant.Dp;...
            Wt.B*Plant.Dp];

HinfSF.E = [zeros(4,1);...
                 -Ws.B;...
                     0];

HinfSF.C = [ Ws.D*Plant.Cp, Ws.C,    0;...
             Wt.D*Plant.Cp,    0, Wt.C;...
            Wc*Cc*Plant.Ap,    0,    0];
            
HinfSF.D1 = [ Ws.D*Plant.Dp;...
              Wt.D*Plant.Dp;...
             Wc*Cc*Plant.Bp];

HinfSF.D2 = [-Ws.D;...
                 0;...
                 0];

HinfSF.B_hat = [HinfSF.B, HinfSF.E];

HinfSF.D_hat = [HinfSF.D1, HinfSF.D2];

HinfSF.Q = HinfSF.C'*HinfSF.C;
HinfSF.S = [HinfSF.C'*HinfSF.D1 HinfSF.C'*HinfSF.D2];

HinfSF.A 
HinfSF.B 
HinfSF.E 
HinfSF.C 
HinfSF.D1 
HinfSF.D2 
HinfSF.B_hat 
HinfSF.D_hat
HinfSF.Q
HinfSF.S


%% Gamma Search Algorithm

gamma_max = 20;
gamma_min = 1;
gamma = gamma_max;
to_test = 1;

while(abs(gamma_max - gamma_min)>1e-10)
    if to_test == 1
        gamma_max = gamma;
    else
        gamma_min = gamma;
    end
    gamma = (gamma_max + gamma_min)/2;   
    disp(['gamma_max = ',num2str(gamma_max,  '%12.10g'),'  gamma = ',num2str(gamma,  '%12.10g'),...
          '  gamma_min = ',num2str(gamma_min,  '%12.10g')])
    
    HinfSF.R = [HinfSF.D1'*HinfSF.D1 HinfSF.D1'*HinfSF.D2;...
                HinfSF.D2'*HinfSF.D1 HinfSF.D2'*HinfSF.D2 - gamma^2];
    
    [P,CL_eig,K] = care(HinfSF.A,HinfSF.B_hat,HinfSF.Q,HinfSF.R,HinfSF.S);
    
    pos_test = 1;
    ev = eig(P);
    for i = 1:6
        if ev(i) < 0
            pos_test = 0;
        end
    end

    HinfSF.Kxopt = -[1 0]*HinfSF.R^-1*(HinfSF.B_hat'*P + HinfSF.S');
    HinfSF.Aref = HinfSF.A + HinfSF.B*HinfSF.Kxopt;
    
    eig_test = 1;
    ev_cl = eig(HinfSF.Aref);
    for i = 1:6
        if ev_cl(i) >= 0
            eig_test = 0;
        end
    end
    to_test = pos_test*eig_test;    
end

disp(' ')
disp('1 ********************')

HinfSF.Kxref = HinfSF.Kxopt;
disp(['Gamma optimal = ',num2str(gamma,  '%12.10g')])
disp(['K optimal = ',num2str(HinfSF.Kxopt)])
Niter = 0;

EE=HinfSF.A'*P+P*HinfSF.A - (P*HinfSF.B_hat+HinfSF.S)*inv(HinfSF.R)*(HinfSF.B_hat'*P+HinfSF.S')+HinfSF.Q;
disp(['Norm EE opt = ',num2str(norm(EE))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set gamm to a value slighly larger than the optimal to get smaller gains
gamma = 1.05 %% gamma chosen such that time response of Hinf matches RSLQR in accel tracking
    HinfSF.R = [HinfSF.D1'*HinfSF.D1 HinfSF.D1'*HinfSF.D2;...
                HinfSF.D2'*HinfSF.D1 HinfSF.D2'*HinfSF.D2 - gamma^2];
    [P,CL_eig,K] = care(HinfSF.A,HinfSF.B_hat,HinfSF.Q,HinfSF.R,HinfSF.S);
    HinfSF.Kxref = -[1 0]*HinfSF.R^-1*(HinfSF.B_hat'*P + HinfSF.S');
    HinfSF.Aref = HinfSF.A + HinfSF.B*HinfSF.Kxref;
disp(['Gamma = ',num2str(gamma, '%12.10g')])
disp(['K ref = ',num2str(HinfSF.Kxref)])
EE=HinfSF.A'*P+P*HinfSF.A - (P*HinfSF.B_hat+HinfSF.S)*inv(HinfSF.R)*(HinfSF.B_hat'*P+HinfSF.S')+HinfSF.Q;
disp(['Norm EE opt = ',num2str(norm(EE))]);
disp(' ')
disp('2 ********************')


%% Closed-Loop System
%*******************************************************
% General Form 
%*******************************************************  
%Close the loop to test the model
% Plant form  xdot = Apx + Bpu; 
%                      y = Cpx +Dpu

  Ap = Plant.Ap
  Bp = Plant.Bp
  Cp = [  Za   0.  Zd  0. ;
          eye(4)]
  Dp = 0*Cp*Bp

% plant size
[nCp,nAp] = size(Cp);
[~,nBp]   = size(Bp);

% Controller xcdot = Acxc + Bc1y + Bc2r
%                      u = Ccxc + Dc1y + Dc2r

  Ac = [ Ws.A       0.;
               0.    Wt.A]
  Bc1 = [ Ws.B*[1 0 0 0 0]
          Wt.B*[1 0 0 0 0]]
  Bc2 = [-Ws.B
                 0.]
  Cc   = [ HinfSF.Kxref(:,5:6) ]
  Dc1 = [ 0. HinfSF.Kxref(:,1:4)]
  Dc2 = [ 0. ]
  
    
%   w=logspace(-1,4,500);
  rtd=180/pi;

  
% Form the closed loop system
     Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
     Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc);
         (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
     Bcl = [       Bp*Z*Dc2;
         (Bc2+Bc1*Dp*Z*Dc2)];
     Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
     Dcl =(Dp*Z*Dc2);
     sys_cl = ss(Acl,Bcl,Ccl,Dcl);


disp('Closed Loop Eigenvalues')
damp(Acl)
[V,D] = eig(Acl);
disp('Closed Loop Eigenvectors')
V
disp('Closed Loop Eigenvalues')
D


% SS model of loop gain at the plant input Lu
     A_Lu = [ Ap 0.*Bp*Cc;  Bc1*Cp Ac];
     B_Lu = [ Bp; Bc1*Dp];
     C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain
     D_Lu = -[ Dc1*Dp];
     sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);
     magdb = 20*log10(abs(squeeze(freqresp(sys_Lu,w))));
     wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
     sr = sigma(sys_Lu,w,3);
     sru_min = min(abs(sr));
     rd = sigma(sys_Lu,w,2);
     rdu_min = min(abs(rd));
     Lu = freqresp(sys_Lu,w);
     
% Analysis at plant output
     T  = freqresp(sys_cl,w); % Complementary Sensitivity
     S = 1 - T; % Sensitivity
     T_Az = 20*log10(abs(squeeze(T(1,1,:))));
     S_Az = 20*log10(abs(squeeze(S(1,1,:))));
     Tmax = max(T_Az); % Inf Norm of T in dB
     Smax = max(S_Az); % Inf Norm of S in dB
     Ws_magdb = 20*log10(abs(squeeze(1/freqresp(Ws.sys,w)))) ;
     Wt_magdb = 20*log10(abs(squeeze(1/freqresp(Wt.sys,w)))) ;



     %Compute singluar value margins
     neg_gm =  min([ (1/(1+rdu_min)) (1-sru_min)]); % in dB
     pos_gm =  max([ (1/(1-rdu_min)) (1+sru_min)]); % in dB
     neg_gmdB = 20*log10( neg_gm ); % in dB
     pos_gmdB = 20*log10( pos_gm ); % in dB
     pm = 180*(max([2*asin(rdu_min/2) 2*asin(sru_min/2)]))/pi;% in deg
     disp('Singular value margins')
     disp(['Min Singular value I+Lu =    ' num2str(rdu_min)])
     disp(['Min Singular value I+invLu = ' num2str(sru_min)])
     disp(['Singular value gain margins = [' ...
         num2str(neg_gmdB) ' dB,' num2str(pos_gmdB) ' dB ]' ])
     disp(['Singular value phase margins = [ +/-' ...
         num2str(pm)  ' deg ]' ])

  
%########################################################################
%LQR Addition
%Define plant matrices without actuator
Ap1 = [Za_V, 1; Ma, 0];
Bp1 = [Zd_V; Md];
Cp1 = [Za, 0; 0, 1];
Dp1 = [Zd; 0]; 

%Define wiggle matrices
o1 = [0; 0];
Aw = [0, Cp1(1,:); o1, Ap1];
Bw = [Dp1(1); Bp1];


%BEGIN LQR
% Setup range of penalties for the LQR
Q=0.*Aw; %sets up a matrix Q that is the size of Aw and is all 0s
R=1;
xeig=[];
%qq is a vector which each element scales the LQR Q matrix
qq=logspace(-6,-1,50); %these are the varying values of q11

xopenloop=eig(Aw);
% Preallocate matrices for speed
az_st = 0.*ones(numel(qq),numel(t));
q_st = 0.*ones(numel(qq),numel(t));
del_st = 0.*ones(numel(qq),numel(t));
deldot_st = 0.*ones(numel(qq),numel(t));


% Loop LQR design increasing the LQR penealty (qq)
% and compute the time domain and frequency domain metrics as the penalty
% is varied
npts = numel(qq);
ii = 1;
%while ii < npts,
% The only penalty is the 1,1 element on the error state
% Desing the gains
Q(1,1)=qq(17);
[Kx_lqr,~,~]=lqr(Aw,Bw,Q,R);

% populate the controller matrices
Ac_act = 0.;
Bc1_act = [1. 0. 0. 0. 0.];
Bc2_act = -1;
Cc_act = -Kx_lqr(1);
Dc1_act = [0. -Kx_lqr(2:3) 0. 0.];
Dc2_act = 0.;

% Form the closed loop system
Z = inv(eye(size(Dc1_act*Dp))-Dc1_act*Dp);
Acl_act = [ (Ap+Bp*Z*Dc1_act*Cp) (Bp*Z*Cc_act);
    (Bc1_act*(Cp+Dp*Z*Dc1_act*Cp)) (Ac_act+Bc1_act*Dp*Z*Cc_act)];
Bcl_act = [ Bp*Z*Dc2_act;
    (Bc2_act+Bc1_act*Dp*Z*Dc2_act)];
Ccl_act = [(Cp+Dp*Z*Dc1_act*Cp) (Dp*Z*Cc_act)];
Dcl_act =(Dp*Z*Dc2_act);
sys_cl_act = ss(Acl_act,Bcl_act,Ccl_act,Dcl_act);
% Compute closed loop eigenvalues for root locus plot
xx=eig(Acl_act);
xeig=[xeig;xx];

% SS model of loop gain at the plant input Lu
A_Lu_act = [ Ap 0.*Bp*Cc_act;  Bc1_act*Cp Ac_act];
B_Lu_act = [ Bp; Bc1_act*Dp];
C_Lu_act = -[ Dc1_act*Cp Cc_act];%change sign for loop gain
D_Lu_act = -[ Dc1_act*Dp];
sys_Lu_act = ss(A_Lu_act,B_Lu_act,C_Lu_act,D_Lu_act);
magdb_act = 20*log10(abs(squeeze(freqresp(sys_Lu_act,w))));
wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
sr_act = sigma(sys_Lu_act,w,3);
sru_min_act = min(abs(sr_act));
rd_act = sigma(sys_Lu_act,w,2);
rdu_min_act = min(abs(rd_act));
Lu_act = freqresp(sys_Lu_act,w);
     
% Analysis at plant output
T_act  = freqresp(sys_cl_act,w); % Complementary Sensitivity
S_act = 1 - T_act; % Sensitivity
T_Az_act = 20*log10(abs(squeeze(T_act(1,1,:))));
S_Az_act = 20*log10(abs(squeeze(S_act(1,1,:))));
Tmax_act = max(T_Az_act); % Inf Norm of T in dB
Smax_act = max(S_Az_act); % Inf Norm of S in dB


Ain = A_Lu_act;
Bin = B_Lu_act;
Cin = C_Lu_act;
Din = D_Lu_act;


% Plant Input Freqeuncy Domain Analysis
for i=1:numel(w),
    s = sqrt(-1)*w(i);
    GG = Cp*inv(s*eye(size(Ap))-Ap)*Bp+Dp;
    KK = Cc_act*inv(s*eye(size(Ac_act))-Ac_act)*Bc1_act+Dc1_act;
    Lu_lqr(i)  = -KK*GG;
    RDu_lqr(i)  = 1. + Lu_lqr(i);
    SRu_lqr(i) = 1. + 1./Lu_lqr(i);
    Lin = Cin*inv(s*eye(size(Ain))-Ain)*Bin+Din;
    % This loops computes SISO freq response at plant input one loop open,
    % rest closed
    for jj = 1:nBp,
        Fyjj = eye(size(Lin));
        Fyjj(jj,jj) = 0.;
        Tujj(jj,i) = inv(eye(size(Lin)) + Lin*Fyjj)*Lin;
    end
end
    
%#######################################################################     


y = step(sys_cl,t);
az = y(:,1); %  acceleration (fps2)
aze = abs(ones(size(az))-az);  % error for az
taur = 0.; taus= 0.; % rise time and settling time
fv = aze(numel(aze)); % final value of the error
e_n = aze - fv*ones(size(aze)) - 0.36*ones(size(aze));
e_n1 = abs(e_n) + e_n;
taur = crosst(e_n1,t); % rise time 
e_n = aze - fv*ones(size(aze)) - 0.05*ones(size(aze));
e_n1 = abs(e_n) + e_n;
taus = crosst(e_n1,t); % settling time

y_act = step(sys_cl_act,t);
az_act = y_act(:,1); %  acceleration (fps2)
aze_act = abs(ones(size(az_act))-az_act);  % error for az
taur_act = 0.; taus_act= 0.; % rise time and settling time
fv_act = aze_act(numel(aze_act)); % final value of the error
e_n_act = aze_act - fv_act*ones(size(aze_act)) - 0.36*ones(size(aze_act));
e_n1_act = abs(e_n_act) + e_n_act;
taur_act = crosst(e_n1,t); % rise time
e_n_act = aze_act - fv_act*ones(size(aze_act)) - 0.05*ones(size(aze_act));
e_n1_act = abs(e_n_act) + e_n_act;
taus_act = crosst(e_n1,t); % settling time

% Scale the linear response by 32.174 for "one gee"
% Plot the time histories
figure('Name','Acceleration Time History')
plot(t,az,'b','LineWidth',2);grid
hold on
plot(t,az_act,'r','LineWidth',2);
hold off
legend('RSLQR response', 'Hinf response');
%legend(['63% Tr = ' num2str(taur) ' 95% Ts = ' num2str(taus)]);
xlabel('Time (sec)');
ylabel('Az (fps2)');
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Pitch Rate Time History')
plot(t,y(:,3)*rtd,'b','LineWidth',2);grid
xlabel('Time (sec)');
ylabel('Pitch Rate (dps)');
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Elevon Time History')
plot(t,y(:,4)*rtd,'b','LineWidth',2);grid
xlabel('Time (sec)');
ylabel('Elevon (deg)');
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Elevon Rate Time History')
plot(t,y(:,5)*rtd,'b','LineWidth',2);grid
xlabel('Time (sec)');
ylabel('Elevon Rate(dps)');
if(save_plots == 1) saveppt2(plot_file_name); end


ngm = -1/neg_gm;
pgm = -1/pos_gm;
figure('Name','Nyquist Plot'),
plot(xx1,yy1,'k:',real(squeeze(Lu)),imag(squeeze(Lu)),'k',...
    [ngm -1.],[0. 0.],'r',[-1 pgm],[0. 0.],'c','LineWidth',2);grid
axis([-3 3 -3 3]);
xlabel('Re(Lu)')
ylabel('Im(Lu)')
legend('Unit Circle at -1,j0','Lu','Neg SV Margin','Pos SV Margin');
if(save_plots == 1) saveppt2(plot_file_name); end

%
figure('Name','Nyquist Plot at Plant Input'),
plot(xx1,yy1,'r',real(squeeze(Lu)),imag(squeeze(Lu)),'k', real(squeeze(Lu_act)),...
imag(squeeze(Lu_act)), 'b', 'LineWidth',2);grid
axis([-2 2 -2 2]);
legend('Unit Circle','Hinf-SF','RSLQR-SF','Location','Best');
xlabel('Re(L)')
ylabel('Im(L)')
title('Nyquist Plot at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Bode Magnitude at Plant Input'),
semilogx(w,20*log10(abs(squeeze(Lu))),'k','LineWidth',2);grid
hold on
semilogx(w,20*log10(abs(squeeze(Lu_act))),'b','LineWidth',2)
hold off
legend('Hinf-SF','RSLQR-SF','Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag')
title('Bode at Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure
margin(sys_Lu);grid
if(save_plots == 1) saveppt2(plot_file_name); end

figure
margin(sys_Lu_act);grid
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Return Difference at Plant Input'),
semilogx(w,20*log10(abs(rd)),'k','LineWidth',2);grid
hold on
semilogx(w,20*log10(abs(rd_act)),'b','LineWidth',2);
hold off
legend(['Hinf min(I+Lu) = ' num2str(rdu_min)], ['RSLQR min(I+Lu) = ' num2str(rdu_min_act)], ...
    'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Stability Robustness at Plant Input'),
semilogx(w,20*log10(abs(sr)),'k', 'LineWidth',2);grid
hold on
semilogx(w,20*log10(abs(sr_act)),'b', 'LineWidth',2);
hold off
legend(['Hinf min(I+invLu) = ' num2str(sru_min)],['RSLQR min(I+invLu) = ' num2str(sru_min_act)], ...
    'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')
if(save_plots == 1) saveppt2(plot_file_name); end

figure('Name','Comp Sens T');
semilogx(w,T_Az,'b',w,Wt_magdb,'k',w,T_Az_act,'g', 'LineWidth',2);grid
legend(['Hinf ||T||inf = ' num2str(Tmax) ' (dB)'],'invWt',...
    ['RSLQR ||T||inf = ' num2str(Tmax_act) ' (dB)'],'Location','Best');
title('Comp Sens T');
xlabel('Freq (rps)');ylabel('Mag (dB)');
if(save_plots == 1) saveppt2(plot_file_name); end  

figure('Name','Sens S');
semilogx(w,S_Az,'b',w,Ws_magdb,'k', w,S_Az_act,'g', 'LineWidth',2);grid
legend(['Hinf ||S||inf = ' num2str(Smax) ' (dB)'],'invWs', ...
    ['RSLQR ||S||inf = ' num2str(Smax_act) ' (dB)'],'Location','Best');
title('Sens S');
xlabel('Freq (rps)');ylabel('Mag (dB)');
if(save_plots == 1) saveppt2(plot_file_name); end   


return     


