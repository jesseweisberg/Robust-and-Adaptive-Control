% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Purpose:
%   Take-home Midterm: RSLQR controller commanding AOA with state feedback
%
% Author: Jesse Weisberg
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% ESE 547 Take-home Midterm
% Author: Jesse Weisberg 
clear all
close all

%% Define Plant Matrices & Wiggle Matrix, Constants %%

% a - alpha, d - delta, w - omega, z - zeta
Za_V = -1.21;
Ma = 44.2506;
Zd_V = -.1987;
Md = -97.3213;
V = 886.78; % (fps)
Za = V*Za_V;
Zd = V*Zd_V;
w_act = pi*14; % actuator natural frequency rps
z_act = 0.6; % actuator damping
grav = 32.174; % gravity in (fps2)

%Define plant matrices without actuator (Used for RSLQR design)
Ap1 = [Za_V 1; Ma 0];
Bp1 = [Zd_V; Md];
Cp1 = [1 0; 0 1];
Dp1 = [0; 0];

%Define wiggle matrices
o1 = [0; 0];
Aw = [0, Cp1(1,:); o1, Ap1];
Bw = [Dp1(1); Bp1];

%Define plant matrices with actuator
Ap2 = [Za_V, 1, Zd_V, 0; 
    Ma, 0, Md, 0; 
    0, 0, 0, 1; 
    0, 0, -w_act*w_act, -2*z_act*w_act];
Bp2 = [0; 0; 0; w_act*w_act];
Cp2 = eye(4);
Dp2 = 0.*Cp2*Bp2;


%% Begin RSLQR Design %%

% Setup range of penalties for the LQR
Q=0.*Aw; % sets up Q as a matrix of 0's that is the size of Aw
R=1; 
xeig=[];
%qq is a vector which each element scales the LQR Q matrix
qq = logspace(-3,2, 500); %these are the varying values of q11, evenly spaced on log scale
w = logspace(-3, 3, 500); %varying frequency values
t = linspace(0, 2, 1000); %time scale  (t = 0:.002:2;)
rtd = 180/pi; %radian to degree scale
% for plotting unit circle in nyquist later
dd = 0.:.001:2*pi;
xx1 = cos(dd)-1; 
yy1 = sin(dd);

xopenloop = eig(Aw); % ??
% Preallocate matrices for speed 
% Stores the step responses of each state for each value of qq
a_st = 0.*ones(numel(qq),numel(t));
q_st = 0.*ones(numel(qq),numel(t));
del_st = 0.*ones(numel(qq),numel(t));
deldot_st = 0.*ones(numel(qq),numel(t));
T_st = 0.*ones(numel(qq),numel(w));
S_st = 0.*ones(numel(qq),numel(w));
sigDif_st = 0.*(numel(qq));

% Loop LQR design increasing the LQR penealty (qq)
% and compute the time domain and frequency domain metrics as the penalty
% is varied
npts = numel(qq);
ii = 372;
while ii < 373
    % The only penalty is the 1,1 element on the error state
    % Design the gains (1 integral error gain and 2 gains for a, q)
    Q(1,1)=qq(ii);
    [Kc,~,~] = lqr(Aw,Bw,Q,R); % Kc - RSLQR gain matrix produced from ARE solution 

    % Define Controller Matrices
    Ac = 0.;
    Bc1 = [1. 0. 0. 0.];
    Bc2 = -1;
    Cc = -Kc(1);
    Dc1 = [-Kc(2:3) 0. 0.];
    Dc2 = 0.;
    
    % Form the closed loop system of RSLQR designed controller with full plant (w/ actuator model)  
    Z = inv(eye(size(Dc1*Dp2))-Dc1*Dp2);
    Acl = [ (Ap2+Bp2*Z*Dc1*Cp2) (Bp2*Z*Cc);
        (Bc1*(Cp2+Dp2*Z*Dc1*Cp2)) (Ac+Bc1*Dp2*Z*Cc)];
    Bcl = [ Bp2*Z*Dc2;
        (Bc2+Bc1*Dp2*Z*Dc2)];
    Ccl = [(Cp2+Dp2*Z*Dc1*Cp2) (Dp2*Z*Cc)];
    Dcl = (Dp2*Z*Dc2);
    sys_cl = ss(Acl,Bcl,Ccl,Dcl);
    
    %Closed Loop Eigenvalues for root locus plot
    disp('Eigenvalues of Closed Loop System')
    xx = eig(Acl)
    xeig=[xeig;xx]; %store eigenvalues for root locus plot later on
    disp('Eigenvectors of Closed Loop System')
    [eigVecs, ~] = eigs(Acl)
    
%     disp('Eigenvalues of Closed Loop System: Check')
%     xx_2 = eig(Aw-Bw*Kc)
    %xeig = [xeig;xx];
    
    %% State-Space model of loop gain at the plant input Lu
    A_Lu = [ Ap2 0.*Bp2*Cc; Bc1*Cp2 Ac];
    B_Lu = [ Bp2; Bc1*Dp2];
    C_Lu = -[ Dc1*Cp2 Cc]; %change sign for loop gain
    D_Lu = -[ Dc1*Dp2];
    sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);
    magdb = 20*log10(abs(squeeze(freqresp(sys_Lu,w))));
    wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
%     T  = freqresp(sys_cl,w); % transfer function of closed loop
%     S = 1 - T;
%     T_st(ii,:) = 20*log10(abs(squeeze(T(1,1,:))));
%     S_st(ii,:) = 20*log10(abs(squeeze(S(1,1,:))));
    
%%%%%%%%%%%%
%     sigma(SYS,W,TYPE) or sigma(SYS,[],TYPE) draws the following modified 
%     SV plots depending on the value of TYPE:
%            TYPE = 1     -->     SV of  inv(SYS)
%            TYPE = 2     -->     SV of  I + SYS
%            TYPE = 3     -->     SV of  I + inv(SYS) 
%%%%%%%%%%%%
    
    %% Tune the RSLQR so that srmin == rdmin (almost)
    sr = sigma(sys_Lu,w,3); % singular values of stability robustness
    srmin = min(abs(sr));   % minimium singular value of stability robustness
    rd = sigma(sys_Lu,w,2); % singular values of return difference
    rdmin = min(abs(rd));   % minimum singular value of return difference
    sigDif = abs(rdmin-srmin);  
%     if abs(rdmin-srmin) < .0005
    if rdmin<.1
        disp('Exit loop: rdmin sufficiently close to srmin')
        ii
        ii = npts;
    end
    sigDif_st(ii) = sigDif; %sigDif value 
    
    %% Time Domain Plots to Form Time Histories for a, q, del, deldot
    y = step(sys_cl,t);
    a = y(:,1); % AOA (radians)
    ae = abs(ones(size(a))-a); % error for a (i.e. 1-a for all a)
    taur = 0.; taus= 0.; % rise time and settling time
    fv = ae(numel(ae)); % final value of the error
    e_n = ae - fv*ones(size(ae)) - 0.36*ones(size(ae)); 
    e_n1 = abs(e_n) + e_n;
    taur = crosst(e_n1,t); % rise time
    e_n = ae - fv*ones(size(ae)) - 0.05*ones(size(ae));
    e_n1 = abs(e_n) + e_n;
    taus = crosst(e_n1,t); % settling time
%     amin = abs(min(a))*100; % undershoot
%     amax = (abs(max(a))-1)*100; % overshoot
%     dmax = max(abs(y(:,4)))*rtd*grav; % compute in per g commanded
%     ddmax = max(abs(y(:,5)))*rtd*grav;
%     metric=[qq(ii) rdmin srmin wc taur taus amin amax dmax ddmax];
%     data(ii,:) = metric;
    a_st(ii,:) = a;
    q_st(ii,:) = y(:,2);
    del_st(ii,:) = y(:,3);
    deldot_st(ii,:) = y(:,4);

    ii = ii+1;
end

%% Find Optimal index, ii, (and thus gain for Q) for RSLQR
[minSigDif, optimal_ii] = min(sigDif_st)

%% Step Responses of Tuned RSLQR
% 
tunedStep = step(sys_cl,t); %step response of tuned RSLQR
a = tunedStep(:,1);
q = tunedStep(:,2);
del = tunedStep(:,3);
deldot = tunedStep(:,4);
aStepInfo = stepinfo(a,t,'RiseTimeLimits',[0.00,0.63], 'SettlingTimeThreshold',.05)
aRT = aStepInfo.RiseTime;
aST = aStepInfo.SettlingTime;

figure('Name','AOA Tuned (Step Response)')
plot(t,a);
legend(['63% Tr =  ' num2str(aRT), ',  95% Ts =  ' num2str(aST)],'Location','Best');
grid on;
xlabel('time (s)');
ylabel('AOA (radians)');
title('AOA Tuned (Step Response)'); 
hold on

figure('Name','Pitch Tuned (Step Response)')
plot(t,q);
grid on;
xlabel('time');
ylabel('Pitch (radians)');
title('Pitch Tuned (Step Response)');
hold on

figure('Name','Elevon Deflection Tuned (Step Response)')
plot(t,del);
grid on;
xlabel('time (s)');
ylabel('Elevon Deflection (radians)');
title('Elevon Deflection Tuned (Step Response)');
hold on

figure('Name','Elevon Deflection Rate Tuned (Step Response)')
plot(t,deldot);
grid on;
xlabel('time (s)');
ylabel('Elevon Deflection Rate (radians/s)');
title('Elevon Deflection Rate Tuned (Step Response)');
hold on

%% Plot Time Histories for a, q, del, deldot
%Stores the step responses of each state for each value of qq

x=2;
figure('Name','AOA Time History (Step Response)')
while x<ii  
    plot(t,a_st(x,:));
    %plot(t,y(:,1),'b', 'LineWidth',2);
    grid on;
    title('AOA Time History (Step Response)');
    xlabel('time (s)');
    ylabel('AOA (radians)');
    hold on
    x = x+1;
end

x=2;
figure('Name','Pitch Time History (Step Response)')
while x<ii  
    plot(t,q_st(x,:));
    %plot(t,y(:,1),'b', 'LineWidth',2);
    grid on;
    title('Pitch Time History (Step Response)');
    xlabel('time (s)');
    ylabel('Pitch (radians)');
    hold on
    x = x+1;
end

x=2;
figure('Name','Elevon Time History (Step Response)')
while x<ii  
    plot(t,del_st(x,:));
    %plot(t,y(:,1),'b', 'LineWidth',2);
    grid on;
    title('Elevon Deflection Time History (Step Response)');
    xlabel('Time (s)');
    ylabel('Elevon Deflection (radians)');
    hold on
    x = x+1;
end

x=2;
figure('Name','Elevon Rate Time History (Step Response)')
while x<ii  
    plot(t,deldot_st(x,:));
    %plot(t,y(:,1),'b', 'LineWidth',2);
    grid on;
    title('Elevon Deflection Rate Time History (Step Response)');
    xlabel('Time (s)');
    ylabel('Elevon Deflection Rate (radians/s)');
    hold on
    x = x+1;
end

%% Plant Input Freqeuncy Domain Analysis (check)
Ain = A_Lu;
Bin = B_Lu;
Cin = C_Lu;
Din = D_Lu;

Lu = freqresp(sys_Lu,w);
T  = freqresp(sys_cl,w);
S = 1 - T;

[nCp, nAp] = size(Cp2);
[~, nBp] = size(Bp2);

for i=1:numel(w)
    s = sqrt(-1)*w(i); % s=jw, iterating through w's
    GG = Cp2*inv(s*eye(size(Ap2))-Ap2)*Bp2+Dp2;
    KK = Cc*inv(s*eye(size(Ac))-Ac)*Bc1+Dc1;
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

%% Frequency Domain Analysis at the Plant Input

figure('Name','Nyquist Plot at Plant Input'),
plot(xx1,yy1,'r',real(squeeze(Lu)),imag(squeeze(Lu)),'k','LineWidth',2);grid
axis([-2 2 -2 2]);
legend('Unit Circle','Loop','Short Form','Location','Best');
xlabel('Re(L)')
ylabel('Im(L)')
title('Nyquist Plot at Plant Input')

figure('Name','Bode Magnitude at Plant Input'),
semilogx(w,20*log10(abs(squeeze(Lu))),'k','LineWidth',2);grid
legend('Loop','Short Form','Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag')
title('Bode at Input')

figure('Name', 'Bode from Matlab'),
margin(sys_Lu)

figure('Name','Return Difference at Plant Input'),
semilogx(w,20*log10(abs(rd)),'k','LineWidth',2);grid
legend([' min(I+Lu) = ' num2str(rdmin)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input')

figure('Name','Stability Robustness at Plant Input'),
semilogx(w,20*log10(abs(sr)),'k','LineWidth',2);grid
legend([' min(I+invLu) = ' num2str(srmin)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')


% Loop Gain/Phase Crossover Frequency  
disp('Classical Margins')
allmargin(sys_Lu) 

%Singular Value Gain/Phase Margins
disp('\\\ Singular Value Margins \\\')
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
   
    
    
    
    
    
    
    
    

