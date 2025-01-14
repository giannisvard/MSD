%% Setup parameters and inputs
load MSD2024_P3_Signals.mat
ts = 30e-6; %Sampling time, s

%% B1
% Create time vector
N = size(d,2);                  % Number of samples
t = linspace(0,(N-1),N)*ts;     % Time vector 

% Plot signals in time domain
figure;
plot(t,n,'b-','LineWidth',2);
hold on
plot(t,d,'r-','LineWidth',2);
title('Provided Signals in the Time Domain')
legend('Output Disturbance','Process Disturbance')
xlabel("t (s)")
ylabel("Signal")
hold off
grid on

% Plot signals in frequency domain
fft_d = fft(d);                 % Obtain the amplitude of d
fft_n = fft(n);                 % Obtain the amplitude of n

figure;
semilogx((0:N-1)/(N*ts),abs(fft_n),'b',"LineWidth",2)
title("Discrete Fourier Transform of Output Disturbance")
xlabel("f (Hz)")
ylabel("|n| (\mum)")
grid on

figure;
semilogx((0:N-1)/(N*ts),abs(fft_d),'r',"LineWidth",2)
title("Discrete Fourier Transform of Process Disturbance")
xlabel("f (Hz)")
ylabel("|d| (V)")
grid on

%% B2
% Input 
