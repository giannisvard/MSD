%% Initializing the process
clf;
clear;
close all;

load MSD2024_P2_Plant.mat % Load the file
s = tf('s');


%% Construct signals
t = linspace(0,10,100);

f_considerd = 10;
w_considered = 2 * pi * f_considerd;

% Input signal r: pure sinusoid
mag_r = 10e-6;            % m
f_r   = f_considerd;               % Hz
w_r   = f_r * 2 * pi;     % rad/s

r = mag_r * sin(w_r * t); % time domain
R = mag_r * 1j * w_r;     % Laplace domain: pure sinusoid


% Process disturbance d: pure sinusoid
mag_d = 20;               % V
f_d   = f_considerd;               % Hz
w_d   = f_d * 2 * pi;     % rad/s

d = mag_d * sin(w_d * t); % time domain
D = mag_d * 1j * w_d;     % Laplace domain: pure sinusoid

% Output disturbance n: pure sinusoid
mag_n = 10e-9;            % m
f_n   = f_considerd;               % Hz
w_n = f_n * 2 * pi;       % rad/s

n = mag_n * sin(w_n * t); % time domain
N = mag_n * 1j * w_n;     % Laplace domain: pure sinusoid



%% Copied controller from Assingment 2


%construct controller
C = controller_assignment_2();
%% Construct system transfer functions

L = C * G;               % Open loop transfer function

S = 1 / (1 + L);         % Sensitivity function

PS = G * S;              % Process * sensitivity function

T = L / (1 + L);         % Complimentary sensitivity function

%% Compute error

S_response = evalfr(S, 1j * w_considered);
PS_response = evalfr(PS, 1j * w_considered);
T_response = evalfr(T, 1j * w_considered);


e_real = sqrt(abs(S_response * R)^2 + abs(PS_response * D)^2 + abs(T_response * N)^2)






