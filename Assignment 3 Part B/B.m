close all
clear
%% Setup parameters and inputs
load MSD2024_P3_Signals.mat
load MSD2024_P2_Plant.mat
ts = 30e-6;                     % Sampling time, s
fs = 1/ts;                      % Sampling frequency, Hz
s = tf('s');                    % Laplace variable
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

% Transform signals to the frequency domain
fft_d_double = fft(d)/N;            % Obtain the amplitude of d, scale by N
fft_n_double = fft(n)/N;            % Obtain the amplitude of n, scale by N
f_vec = (0:N-1)/(N*ts);             % Obtain frequency vector (unshifted)
f_vec_p = f_vec(1:N/2);             % Frequency vector for positive frequencies only
% Obtain single-sided spectrum
fft_d = zeros(1,N/2);
fft_n = zeros(1,N/2);
% Accumulate to positive frequencies
fft_d = 2*fft_d_double(1:N/2);
fft_n = 2*fft_n_double(1:N/2);
fft_d(1) = fft_d_double(1);         % DC is excluded
fft_n(1) = fft_n_double(1);         % DC is excluded

% Plot signals in frequency domain
figure;
semilogx(f_vec_p,abs(fft_n),'b',"LineWidth",2)
title("Discrete Fourier Transform of Output Disturbance")
xlabel("f (Hz)")
ylabel("|n| (\mum)")
grid on

figure;
semilogx(f_vec_p,abs(fft_d),'r',"LineWidth",2)
title("Discrete Fourier Transform of Process Disturbance")
xlabel("f (Hz)")
ylabel("|d| (V)")
grid on

%% B2
% Compute single-sided PSD
psd_d = (abs(fft_d).^2)/(N*fs/2);    % Single-sided PSD
psd_n = (abs(fft_n).^2)/(N*fs/2);    % Single-sided PSD
% Plot the results
figure;
semilogx(f_vec_p,psd_d,'r',"LineWidth",2);
title("PSD of Process Disturbance")
xlabel("f (Hz)")
ylabel("PSD (V^2/Hz)")
grid on

figure;
semilogx(f_vec_p,psd_n,'b',"LineWidth",2);
title("PSD of Output Disturbance")
xlabel("f (Hz)")
ylabel("PSD (\mum^2/Hz)")
grid on

%% B3
% Obtain system's transfer functions
C = controller_assignment_2;            % Controller from assignment 2
P = G;                                  % Define plant tf as P
GS = P/(1+P*C);                         % Process sensitivity
S = 1/(1+P*C);                          % Output sensitivity
% Calculate function frequency response
[A_GS,f_GS] = freqresp(GS,f_vec_p,'Hz');    % Freqresponse for GS
A_GS = squeeze(A_GS)';                      % Remove empty dimensions
[A_S,f_S] = freqresp(S,f_vec_p,'Hz');       % Freqresponse for S
A_S = squeeze(A_S)';                        % Remove empty dimensions
% Calculate PSD of y
psd_y = ((abs(A_GS).^2).*psd_d) + ((abs(A_S).^2).*psd_n);   % PSD of y, the result of this does not make sense
% Plot PSD of y
figure;
semilogx(f_vec_p,psd_y,'k',"LineWidth",2);
title("PSD of Output y")
xlabel("f (Hz)")
ylabel("PSD (\mum^2/Hz)")
grid on
figure
loglog(f_vec_p,abs(A_S))
title("Response of S")
figure
loglog(f_vec_p,abs(A_GS))
title("Response of GS")
%% B4
% Calculate Cumulative Power Spectrum
cps_y = cumsum(psd_y.*f_vec_p(2));
cps_d = cumsum(psd_d.*f_vec_p(2));
cps_n = cumsum(psd_n.*f_vec_p(2));
% Plot PSD of y
figure;
loglog(f_vec_p,cps_y,'k',"LineWidth",2);
hold on
loglog(f_vec_p,cps_d,'r',"LineWidth",2);
hold on
loglog(f_vec_p,cps_n,'b',"LineWidth",2);
hold off
title("CPS")
xlabel("f (Hz)")
ylabel("CPS (\mum^2)")
legend('Output','Process Disturbance','Output Disturbance')
grid on