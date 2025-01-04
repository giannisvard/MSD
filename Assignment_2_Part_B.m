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

% e/r   Sensitivity
[T1, f1] = tfestimate(input, output1, window, overlap, ft, Fs); % Pass overlap and window explicitly
[C1, f1] = mscohere(input, output1, window, overlap, ft, Fs)

% u/r   Noise/Control Sensitivity
[T2, f2] = tfestimate(input, output2, window, overlap, ft, Fs); % Pass overlap and window explicitly
[C2, f2] = mscohere(input, output2, window, overlap, ft, Fs)


% y/r   Complimentary Sensitivity
[T3, f3] = tfestimate(input, output3, window, overlap, ft, Fs); % Pass overlap and window explicitly
[C3, f3] = mscohere(input, output3, window, overlap, ft, Fs)


figure;
%e/r TF and coherence
subplot(3,1,1);
semilogx(f1,mag2db(abs(T1))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Magnitude e/r: Sensitivity Function');

subplot(3,1,2);
semilogx(f1,rad2deg(angle(T1))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on;
title('Phase e/r: Sensitivity Function');

subplot(3,1,3);
semilogx(f1,C1); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;
title('Coherence e/r: Sensitivity Function');

figure;
%u/r TF and coherence
subplot(3,1,1);
semilogx(f2,mag2db(abs(T2))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Magnitude u/r: Control Sensitivity Function');

subplot(3,1,2);
semilogx(f2,rad2deg(angle(T2))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Phase (deg');
grid on;
title('Phase u/r: Control Sensitivity Function');

subplot(3,1,3);
semilogx(f2,C2); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;
title('Coherence u/r: Control Sensitivity Function');

figure;
%y/r TF and coherence
subplot(3,1,1);
semilogx(f3,mag2db(abs(T3))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Magnitude y/r: Complimentary Sensitivity Function');

subplot(3,1,2);
semilogx(f3,rad2deg(angle(T3))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on;
title('Phase y/r: Complimentary Sensitivity Function');

subplot(3,1,3);
semilogx(f3,C3); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;
title('Coherence y/r: Complimentary Sensitivity Function');

%% Problem B4: calculate and present freq response of plant and controller
%CT TF for S

C = T2 .* (1 + T3 ./ (1 - T3));
P = (1 - T1) ./ (T1 .* C);

figure;
subplot(2,1,1);
semilogx(f1,mag2db(abs(C))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Identified Controller Magnitude');

subplot(2,1,2);
semilogx(f1,rad2deg(angle(C))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on;
title('Identified Controller Phase');



figure;
semilogx(f3,mag2db(abs(P))); hold on;
xlim([1,max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('test P');





% figure(2);clf(2);
% subplot(3,1,1);semilogx(f,mag2db(abs(T))); grid on; hold on;
% title('Transfer Function from tfestimate (Magnitude)');
% subplot(3,1,2);semilogx(f,rad2deg(angle(T))); grid on; hold on;
% title('Transfer Function from tfestimate (Phase');
% subplot(3,1,3);semilogx(f,C); grid on; hold on;
% title('Coherence function from mscohere');





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










