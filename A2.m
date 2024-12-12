load MSD2024_P2_Plant.mat % Load the file

s = tf('s');

%Bode Plant

% figure;
% bode(G); grid on;
% title('System Response');

%% tune controller variables
wc = 16.07; %radians

Kp = 0.5;

wi = 230 * 2 * pi; %radians

wd =  130 * 2 * pi; %radians

wt = 90000 * 2 * pi; %radians


%construct controller
C = Kp * (1 + wi/s) * (s/wd + 1)/(s/wt + 1);

% bode(C)

% figure;
% bode(C); grid on;
% title('Controller Response');

%% Open loop 
L = G*C
poles = pole(L);
disp(poles);
% figure;
% bode(L); grid on;
% margin(L);
% title('Open Loop Response');

%% Sensitivity stuff
T = L / (1 + L);
S = 1 / (1 + L);

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

info = stepinfo(Tsum);

info

%End of Closed Loop Controller


%% Feedforward Definition
F = inv(G);

% figure;
% bode(F); grid on;
% title('Feedforward Controller Response');


disp(isproper(F));
%its not proper..
%make proper using low/pass filter
wclp = 2 * wc;
LowPass = 1 / (((s / wclp) + 1)^2)    ;

disp(LowPass);
LowPass
PropaF = LowPass * F;
% figure;
% bode(PropaF); grid on;
% title('Proper Feedforward Controller Response');

[num, den] = tfdata(F, 'v'); % Replace F with your feedforward transfer function

% Compare degrees
deg_num = length(num) - 1; % Degree of numerator
deg_den = length(den) - 1; % Degree of denominator

deg_num
deg_den