%% Initializing the process
clear;
load MSD2024_P2_Plant.mat % Load the file
s = tf('s');
% Define options for Bode plot
opts = bodeoptions;
%opts.FreqUnits = 'Hz'; % Change frequency units to Hertz, commented out
%due to 

% figure;
% bode(G); grid on;
% title('System Response');

%% Feedback Controller tuning
% Shapeit parameters
Kp = 0.5;
fi = 230; %Hz
fd =  130 ; %Hz
ft = 90000; %Hz
% Convert to rad/s
wi = fi*2*pi; %Integrator zero, rad/s
wd = fd*2*pi; %Differentiator pole, rad/s
wt = ft*2*pi; %(Tamed)Differentiator zero, rad/s
%Create PID transfer function
C_PID = (Kp + wi/s + ((s/wd) + 1) / ((s/wt) + 1));      %proportional, integral, lead/lag

% Regular notch filter implementation
f1 = 738;           % Notch frequency (Hz)
zeta_1 = 0.01;      % Damping for numerator (sharp)
zeta_2 = 0.9;       % Damping for denominator (wider)
% Convert to rad/s
w1 = f1*2*pi; %rad/s
% Create transfer function
num_notch = [1, 2*zeta_1*w1, w1^2];  % Numerator coefficients
den_notch = [1, 2*zeta_2*w1, w1^2];  % Denominator coefficients
notch = tf(num_notch, den_notch);

% Skewed notch implementation
f21 = 1009;          % Frequency for numerator (Hz)
zeta_21 = 0.01;      % Damping ratio for numerator
f22 = 971.7;         % Frequency for denominator (Hz)
zeta_22 = 0.005;     % Damping ratio for denominator
% Convert to rad/s
w21 = f21*2*pi; %Notch zero, rad/s
w22 = f22*2*pi; %Notch pole, rad/s
% Create transfer function
num_skewed = [1, 2*zeta_21*w21, w21^2];  % Numerator coefficients
den_skewed = [1, 2*zeta_22*w22, w22^2];  % Denominator coefficients
skewed_notch = tf(num_skewed, den_skewed);

%construct controller
C = C_PID * notch * skewed_notch;

figure;
bode(C); grid on;
title('Controller Response');

%% Open loop 

disp(opts);

L = G*C;
poles = pole(L);
disp(poles);
figure;
bodeplot(L, opts); 
grid on;
margin(L);
title('Open Loop Response');

%% Closed loop transfer functions
T = L / (1 + L);
S = 1 / (1 + L);


figure;
bode(T); grid on;
margin(T);
title('Complimentary Sensitivity');

figure;
bode(S); grid on;
margin(S);
title('Sensitivity');
% Might have to add nyquist plots to evaluate modulus margin here!!!

%% Feedback step response
G.u='u'; G.y='x';
C.u='e'; C.y='v'; 

S1 = sumblk("e = r - y");
S2 = sumblk("u = d + v");
S3 = sumblk("y = x + n");

Tsuma = connect(G,C,S1,S2,S3,"r","y");


figure;
step(Tsuma); grid on;
title('Closed Loop Step Response (Feedback only)');
info_step_feedback = stepinfo(Tsuma);

%% Feedforward Controller tuning
% New controller parameters for ff
% Gain
Kp = 1.91;
% Skewed Notch
fz1 = 739; %Skewed notch zero, Hz
zz1 = 0.009; %Skewed notch zero damping coefficient
fp1 = 972; %Skewed notch pole, Hz
zp1 = 0.005; %Skewed notch pole damping coefficient
% Notch
f2 = 1008; %Notch zero, Hz
zz2 = 0.014; %Notch zero damping coefficient
zp2 = 0.41; %Notch pole damping coefficient
% Lowpass
flp = 1003; %Low pass zero, Hz
% Convert shapeit values to angular frequencies and quality factors
wz1 = fz1*2*pi; % rad/s
wp1 = fp1*2*pi; % rad/s
w2 = f2*2*pi; % rad/s
wlp = flp*2*pi; % rad/s

Qz1 = 1/(2*zz1); %Quality factor (Q1)
Qp1 = 1/(2*zp1); %Quality factor (Q2)
Qz2 = 1/(2*zz2); %Quality factor (Q1)
Qp2 = 1/(2*zp2); %Quality factor (Q2)

% Controller parameters
C_skewednotch = ((s/wz1)^2+s/(Qz1*wz1)+1)/((s/wp1)^2+s/(Qp1*wp1)+1);
C_notch = ((s/w2)^2+s/(Qz2*w2)+1)/((s/w2)^2+s/(Qp2*w2)+1);
C_lowpass = 1/(s/wlp+1);
Cff = C_skewednotch*C_notch*C_lowpass;

figure;
bode(Cff); grid on;
title('Feedforward Controller Response');

if isproper(Cff)
    disp('Feedforward controller is proper! Check response for strictly proper')
end
%IS PROPER VERY NICE! (Borat voice)

%% Combined step response
Cff.u='r'; Cff.y='f';

S1 = sumblk("e = r - y");
S2 = sumblk("u = d + v + f");
S3 = sumblk("y = x + n");

Tsumb = connect(G,C,Cff,S1,S2,S3,"r","y");

figure;
step(Tsumb); grid on;
title('Closed Loop Step Response (Feedback & Feedforward)');
info_step_feedfoward_feedback = stepinfo(Tsumb);

%% Architecture comparison
% Start with reference tracking
figure;
bode(Tsuma); grid on;
title('Reference Tracking (Feedback only)');

figure;
bode(Tsumb); grid on;
title('Reference Tracking (Feedback & Feedforward)');

% Disturbance rejection y/d
GSa = connect(G,C,S1,S2,S3,"d","y");
GSb = connect(G,C,Cff,S1,S2,S3,"d","y");

figure;
bode(GSa); grid on;
title('Disturbance Rejection (Feedback only)');

figure;
bode(GSb); grid on;
title('Disturbance Rejection (Feedback & Feedforward)');

% Noise Attenuation y/n
Sa = connect(G,C,S1,S2,S3,"n","y");
Sb = connect(G,C,Cff,S1,S2,S3,"n","y");

figure;
bode(Sa); grid on;
title('Noise Attenuation (Feedback only)');

figure;
bode(Sb); grid on;
title('Noise Attenuation (Feedback & Feedforward)');