
clear;

load MSD2024_P2_Plant.mat % Load the file

s = tf('s');
% Define options for Bode plot
opts = bodeoptions;
opts.FreqUnits = 'Hz'; % Change frequency units to Hertz

% figure;
% bode(G); grid on;
% title('System Response');

%% tune controller variables
wc = 16.07; %radians

Kp = 0.5;

wi = 230; %radians

wd =  130 ; %radians

wt = 90000; %radians

C_PID = (Kp + wi/s + ((s/wd) + 1) / ((s/wt) + 1))      %proportional, integral, lead/lag


% w11 = 738;
% zeta_11 = 0.01 ;
% w12 = 738;
% zeta_12 = 0.9 ;
% Q11 = 1/(2*zeta_11);
% Q12 = 1/(2*zeta_12);
% notch = ((s/w11)^2 + (s/(Q11*w11)) + 1) / ((s/w12)^2 + (s/(Q12*w12)) + 1);

w1 = 738;           % Notch frequency (rad/s)
zeta_1 = 0.01;      % Damping for numerator (sharp)
zeta_2 = 0.9;       % Damping for denominator (wider)
num_notch = [1, 2*zeta_1*w1, w1^2];  % Numerator coefficients
den_notch = [1, 2*zeta_2*w1, w1^2];  % Denominator coefficients
notch = tf(num_notch, den_notch);




% w21 = 1009;
% zeta_21 = 0.01;
% w22 = 971.7;
% zeta_22 = 0.005;
% Q21 = 1/(2*zeta_21);
% Q22 = 1/(2*zeta_22);
% skewed_notch = ((2/w21)^2 + (s/(Q21*w21)) + 1) / ((s/22)^2 + (s/(Q22*w22)) + 1);

w21 = 1009;          % Frequency for numerator (rad/s)
zeta_21 = 0.01;      % Damping ratio for numerator
w22 = 971.7;         % Frequency for denominator (rad/s)
zeta_22 = 0.005;     % Damping ratio for denominator
num_skewed = [1, 2*zeta_21*w21, w21^2];  % Numerator coefficients
den_skewed = [1, 2*zeta_22*w22, w22^2];  % Denominator coefficients

skewed_notch = tf(num_skewed, den_skewed);


%construct controller
C = C_PID * notch * skewed_notch;

bode(C);

% figure;
% bode(C); grid on;
% title('Controller Response');

%% Open loop 
% 
% disp(opts);
% 
% L = G*C;
% poles = pole(L);
% disp(poles);
% figure;
% bodeplot(L, opts); 
% grid on;
% margin(L);
% title('Open Loop Response');

% %% Sensitivity stuff
% T = L / (1 + L);
% S = 1 / (1 + L);

% 
% figure;
% bode(T); grid on;
% margin(T);
% title('Complimentary Sensitivity');
% 
% figure;
% bode(S); grid on;
% margin(S);
% title('Sensitivity');


%% Step response
G.u='u'; G.y='x';
C.u='e'; C.y='v';

S1 = sumblk("e = r - y");
S2 = sumblk("u = d + v");
S3 = sumblk("y = x + n");

Tsum = connect(G,C,S1,S2,S3,"r","y");

% 
% figure;
% step(Tsum); grid on;
% title('Closed Loop Step Response');

% info = stepinfo(Tsum);
% 
% info

%End of Closed Loop Controller


%% Feedforward Definition
F = inv(G);

% figure;
% bode(F); grid on;
% title('Feedforward Controller Response');


disp(isproper(F));
%its not proper..
%make proper using low/pass filter
% wclp = 2 * wc;
% LowPass = 1 / (((s / wclp) + 1)^2)    ;

% disp(LowPass);
% LowPass
% PropaF = LowPass * F;
% figure;
% bode(PropaF); grid on;
% title('Proper Feedforward Controller Response');

% [num, den] = tfdata(F, 'v'); % Replace F with your feedforward transfer function
% 
% % Compare degrees
% deg_num = length(num) - 1; % Degree of numerator
% deg_den = length(den) - 1; % Degree of denominator

% deg_num
% deg_den