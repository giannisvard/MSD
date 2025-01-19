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
delta_f = 1/(N*ts);             % Frequency resolution

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
f_vec = (0:N/2)/(N*ts);             % Obtain frequency vector (unshifted)

%f_vec_p = f_vec(1:N/2);             % Frequency vector for positive frequencies only
% Obtain single-sided spectrum
fft_d = fft_d_double(1:N/2+1);
fft_n = fft_n_double(1:N/2+1);
fft_d(2:end-1) = 2*fft_d(2:end-1);         % DC and nyq is excluded
fft_n(2:end-1) = 2*fft_n(2:end-1);         % DC and nyq is excluded

% Plot signals in frequency domain
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
% Compute double-sided PSD
psd_d_double = (abs(fft(d)).^2)/(N*fs);    % Double-sided PSD
psd_n_double = (abs(fft(n)).^2)/(N*fs);    % Double-sided PSD
% Compute single-sided PSD
% Compute double-sided PSD
psd_d = psd_d_double(1:N/2+1);    % Single-sided PSD
psd_d(2:end-1) = 2*psd_d(2:end-1);
psd_n = psd_n_double(1:N/2+1);    % Single-sided PSD
psd_n(2:end-1) = 2*psd_n(2:end-1);
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
[A_GS,~] = freqresp(GS,2*pi*f_vec);    % Freqresponse for GS
A_GS = squeeze(A_GS)';                      % Remove empty dimensions
[A_S,~] = freqresp(S,2*pi*f_vec);       % Freqresponse for S
A_S = squeeze(A_S)';                        % Remove empty dimensions
% Calculate PSD components
H2_d = abs(A_GS).^2;
H2_n = abs(A_S).^2;
% Calculate PSD of y
psd_y = (H2_d.*psd_d) + (H2_n.*psd_n);   % PSD of y, the result of this does not make sense
% Plot PSD of y
figure;
semilogx(f_vec,psd_y,'k',"LineWidth",2);
title("PSD of Output y")
xlabel("f (Hz)")
ylabel("PSD (\mum^2/Hz)")
grid on
%% B4
% Calculate Cumulative Power Spectrum
cps_y = cumsum(psd_y)*delta_f;
cps_d = cumsum(H2_d.*psd_d)*delta_f;
cps_n = cumsum(H2_n.*psd_n)*delta_f;
% Plot PSD of y
figure;
loglog(f_vec,cps_y,'k',"LineWidth",2);
hold on
loglog(f_vec,cps_d,'r',"LineWidth",2);
hold on
loglog(f_vec,cps_n,'b',"LineWidth",2);
hold off
title("CPS")
xlabel("f (Hz)")
ylabel("CPS (\mum^2)")
legend('Output','Process Disturbance','Output Disturbance')
grid on

% Define the frequency range
f_min = 1e2; % 100 Hz
f_max = 1e4; % 10,000 Hz

% Filter the data to include only the specified frequency range
valid_indices = (f_vec >= f_min) & (f_vec <= f_max);
f_filtered = f_vec(valid_indices);
cps_y_filtered = cps_y(valid_indices);
cps_d_filtered = cps_d(valid_indices);
cps_n_filtered = cps_n(valid_indices);

% Plot PSD within the specified frequency range
figure;
loglog(f_filtered, cps_y_filtered, 'k', "LineWidth", 2); % Output
hold on;
loglog(f_filtered, cps_d_filtered, 'r', "LineWidth", 2); % Process Disturbance
loglog(f_filtered, cps_n_filtered, 'b', "LineWidth", 2); % Output Disturbance
hold off;

% Add title, labels, legend, and grid
title("CPS (Filtered Frequency Range)");
xlabel("f (Hz)");
ylabel("CPS (\mum^2)");
legend('Output', 'Process Disturbance', 'Output Disturbance');
grid on;

%% B5
% Import updated controller
C_new = revised_controller(C,G);

% Update transfer functions
GS = P/(1+P*C);                         % Process sensitivity
S = 1/(1+P*C);                          % Output sensitivity
% Calculate function frequency response
[A_GS,~] = freqresp(GS,2*pi*f_vec);    % Freqresponse for GS
A_GS = squeeze(A_GS)';                      % Remove empty dimensions
[A_S,~] = freqresp(S,2*pi*f_vec);       % Freqresponse for S
A_S = squeeze(A_S)';                        % Remove empty dimensions
% Calculate PSD components
H2_d = abs(A_GS).^2;
H2_n = abs(A_S).^2;
% Calculate PSD of y
psd_y = (H2_d.*psd_d) + (H2_n.*psd_n);   % PSD of y, the result of this does not make sense
% Plot PSD of y
figure;
semilogx(f_vec,psd_y,'k',"LineWidth",2);
title("PSD of Output y")
xlabel("f (Hz)")
ylabel("PSD (\mum^2/Hz)")
grid on

% Calculate Cumulative Power Spectrum
cps_y_updated = cumsum(psd_y)*delta_f;
cps_d_updated = cumsum(H2_d.*psd_d)*delta_f;
cps_n_updated = cumsum(H2_n.*psd_n)*delta_f;
% Plot PSD of y
figure;
loglog(f_vec,cps_y_updated,'k',"LineWidth",2);
hold on
loglog(f_vec,cps_d_updated,'r',"LineWidth",2);
hold on
loglog(f_vec,cps_n_updated,'b',"LineWidth",2);
hold off
title("CPS (Updated)")
xlabel("f (Hz)")
ylabel("CPS (\mum^2)")
legend('Output','Process Disturbance','Output Disturbance')
grid on

% Filter the data to include only the specified frequency range
cps_y_filtered = cps_y_updated(valid_indices);
cps_d_filtered = cps_d_updated(valid_indices);
cps_n_filtered = cps_n_updated(valid_indices);

% Plot PSD within the specified frequency range
figure;
loglog(f_filtered, cps_y_filtered, 'k', "LineWidth", 2); % Output
hold on;
loglog(f_filtered, cps_d_filtered, 'r', "LineWidth", 2); % Process Disturbance
loglog(f_filtered, cps_n_filtered, 'b', "LineWidth", 2); % Output Disturbance
hold off;

% Add title, labels, legend, and grid
title("CPS (Updated, Filtered Frequency Range)");
xlabel("f (Hz)");
ylabel("CPS (\mum^2)");
legend('Output', 'Process Disturbance', 'Output Disturbance');
grid on;

%% Debugging for B5

Delta_cps = cps_y_updated-cps_y;

figure;
plot(f_vec,Delta_cps)
title('CPS diff')
