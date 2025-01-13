%% Setup parameters and inputs
load MSD2024_P3_Signals.mat
ts = 30e-6; %Sampling time, s

%% B1
% Create time vector
N = size(d,2);
t = linspace(0,ts*(N-1),N); 

% Plot signals in time domain
figure;
plot(t,n,'b-','LineWidth',2);
hold on
plot(t,d,'r-','LineWidth',2);
title('Provided time signals')
legend('Process Disturbance', 'Output Disturbance')
hold off
grid on

% Plot signals in frequency domain
fft_d = fft(d); %Obtain the magnitudes
figure;
plot(fft_d)
