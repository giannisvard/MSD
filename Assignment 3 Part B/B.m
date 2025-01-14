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

% Plot signals in frequency domain
fft_d = fft(d);                 % Obtain the amplitude of d
fft_n = fft(n);                 % Obtain the amplitude of n
f_vec = (0:N-1)/(N*ts);         % Obtain frequency vector

figure;
semilogx(f_vec,abs(fft_n),'b',"LineWidth",2)
title("Discrete Fourier Transform of Output Disturbance")
xlabel("f (Hz)")
ylabel("|n| (\mum)")
grid on

figure;
semilogx(f_vec,abs(fft_d),'r',"LineWidth",2)
title("Discrete Fourier Transform of Process Disturbance")
xlabel("f (Hz)")
ylabel("|d| (V)")
grid on

%% B2
% Compute single-sided PSD
psd_d = (abs(fft_d).^2)/(N*fs);         % Single-sided PSD
psd_n = (abs(fft_n).^2)/(N*fs);    % Single-sided PSD
% Plot the results
figure;
semilogx(f_vec,psd_d,'r',"LineWidth",2);
title("PSD of Process Disturbance")
xlabel("f (Hz)")
ylabel("PSD (V^2/Hz)")
grid on

figure;
semilogx(f_vec,psd_n,'b',"LineWidth",2);
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
[A_GS,f_GS] = freqresp(GS,f_vec,'Hz'); % Freqresponse for GS
A_GS = squeeze(A_GS)';                          % Remove empty dimensions
[A_S,f_S] = freqresp(S,f_vec,'Hz');    % Freqresponse for S
A_S = squeeze(A_S)';                            % Remove empty dimensions
% Calculate PSD of y
psd_y = abs(A_GS).^2.*psd_d + abs(A_S).^2.*psd_n;   % PSD of y
% Plot PSD of y
figure;
loglog(f_vec,psd_y,'k',"LineWidth",2);     % Use log on both axes due to large variations
title("PSD of Output y")
xlabel("f (Hz)")
ylabel("PSD (\mum^2/Hz)")
grid on

%% B4
% Calculate Cumulative Power Spectrum
cps_y = cumsum(psd_y.*f_vec(2));
% Plot PSD of y
figure;
loglog(f_vec,cps_y,'k',"LineWidth",2);  % Use log on both axes due to large variations
title("CPS of Output y")
xlabel("f (Hz)")
ylabel("CPS (\mum^2)")
grid on
% Calculate complimentary sensitivity
T = P*C/(1+P*C);                        % tf format
[A_T,f_T] = freqresp(T,f_vec,'Hz');     % Freqresponse for T
A_T = squeeze(A_T)';                    % Remove empty dimensions
% Error contributions
e_d = abs(A_GS).*d;                     % Process disturbance error
e_n = abs(A_T).*n;                      % Output disturbance error
e_total = sqrt(e_d.^2+e_n.^2);            % Total error
% Convert error performance from time to frequency domain
E_d = fft(e_d);                         % Discrete Fourier Tranform of e_d
E_n = fft(e_n);                         % Discrete Fourier Tranform of e_n
E_total = fft(e_total);                 % Discrete Fourier Tranform of e_total
% Plot error contributions
figure;
plot(f_vec,E_total,'k',"LineWidth",2);
hold on
plot(f_vec,E_d,'r',"LineWidth",2);
hold on
plot(f_vec,E_n,'b',"LineWidth",2);
hold off
title("Error of Output y")
xlabel("f (Hz)")
ylabel("e (\mum^2)")
legend('Total error','Process Disturbance error','Output Disturbance error')