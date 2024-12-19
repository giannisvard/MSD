clf;
clear;
close all;

load MSD2024_P2_Signals.mat % Load the file

%% Problem B.2 plot the signals in time domain
Ts = 30e-6;  %seconds
Fs = 1/Ts ;
samples = length(e); % signals are equally long


time = (0:samples-1) * Ts; % time vector is sample amount times sample time


figure;

subplot(2,2,1);
plot(time, r, 'LineWidth', 1.5);
grid on;
title('Signal r in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');    

subplot(2,2,2);
plot(time, e, 'LineWidth', 1.5);
grid on;
title('Signal e in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,2,3);
plot(time, u, 'LineWidth', 1.5);
grid on;
title('Signal u in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,2,4);
plot(time, y, 'LineWidth', 1.5);
grid on;
title('Signal y in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');


%% Problem B.3 Calculate and present freq response

% Define the window and overlap
window_length = round(length(e) / 9);    % 9 chirp signals, so the windows will only cover the chirp signals
window = rectwin(window_length);         % Use rectwin window
overlap = round(0.5 * window_length); % 50% overlap

%THERE ARE 10 SEGMENTS IN THE SIGNALS (1 STARTER AND 9 CHRIP RESPONSES) So
%this line divides the total input signal
ft = max(round(samples / 10), window_length);


input = r;
output1 = e;
output2 = u;
output3 = y;

%CHECK HERTZ OR RAD/S OUTPUT
Fs_rads = 2 * pi / Ts; % Sampling frequency in rad/s

% e/r
[T1, f1] = tfestimate(input, output1, window, overlap, ft, Fs); % Pass overlap and window explicitly
[C1, f1] = mscohere(input, output1, window, overlap, ft, Fs)

% u/r
[T2, f2] = tfestimate(input, output2, window, overlap, ft, Fs); % Pass overlap and window explicitly
[C2, f2] = mscohere(input, output2, window, overlap, ft, Fs)


% y/r   
[T3, f3] = tfestimate(input, output3, window, overlap, ft, Fs); % Pass overlap and window explicitly
[C3, f3] = mscohere(input, output3, window, overlap, ft, Fs)


figure;
%e/r TF and coherence
subplot(3,2,1);
semilogx(f1,mag2db(abs(T1))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('TF e/r');

subplot(3,2,2);
semilogx(f1,mag2db(abs(C1))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Coherence e/r');


%u/r TF and coherence
subplot(3,2,3);
semilogx(f2,mag2db(abs(T2))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('TF u/r');

subplot(3,2,4);
semilogx(f2,mag2db(abs(C2))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Coherence u/r');


%y/r TF and coherence
subplot(3,2,5);
semilogx(f3,mag2db(abs(T3))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('TF y/r');

subplot(3,2,6);
semilogx(f3,mag2db(abs(C3))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Coherence y/r');

%% Problem B4: calculate and present freq response of plant and controller
  







%% Debug
% disp(['Sampling Frequency (Fs): ', num2str(Fs), ' Hz']);
% disp(['Nyquist Frequency: ', num2str(Fs/2), ' Hz']);
% 
% disp(['Frequency Vector (f): ', num2str(min(f1)), ' to ', num2str(max(f1)), ' Hz']);
% 
% 
% % Convert frequency vector from radians/sample to Hz
% f1_Hz = f1 * (Fs / (2 * pi));
% disp(['Frequency Vector (Hz): ', num2str(min(f1_Hz)), ' to ', num2str(max(f1_Hz)), ' Hz']);










