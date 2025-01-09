clf;
clear;
close all;

load MSD2024_P2_Signals.mat % Load the file
load MSD2024_P2_Plant.mat


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
window_length = round(length(e) / 11);    % 9 chirp signals, so the windows will only cover the chirp signals
window = hann(window_length);         % Use rectwin window
overlap = round(0.8 * window_length); % 50% overlap


%THERE ARE 10 SEGMENTS IN THE SIGNALS (1 STARTER AND 9 CHRIP RESPONSES) So
%this line divides the total input signal
% ft = max(round(samples / 10), window_length);
ft = 10000;

input = r;
output1 = e;
output2 = u;
output3 = y;


% window = hann(10000);
% overlap = round(10000/2);
% ft = [];
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
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Magnitude e/r: Sensitivity Function');

subplot(3,1,2);
semilogx(f1,rad2deg(angle(T1))); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on;
title('Phase e/r: Sensitivity Function');

subplot(3,1,3);
semilogx(f1,C1); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;
title('Coherence e/r: Sensitivity Function');

figure;
%u/r TF and coherence
subplot(3,1,1);
semilogx(f2,mag2db(abs(T2))); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Magnitude u/r: Control Sensitivity Function');

subplot(3,1,2);
semilogx(f2,rad2deg(angle(T2))); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Phase (deg');
grid on;
title('Phase u/r: Control Sensitivity Function');

subplot(3,1,3);
semilogx(f2,C2); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;
title('Coherence u/r: Control Sensitivity Function');

figure;
%y/r TF and coherence
subplot(3,1,1);
semilogx(f3,mag2db(abs(T3))); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Magnitude y/r: Complimentary Sensitivity Function');

subplot(3,1,2);
semilogx(f3,rad2deg(angle(T3))); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on;
title('Phase y/r: Complimentary Sensitivity Function');

subplot(3,1,3);
semilogx(f3,C3); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;
title('Coherence y/r: Complimentary Sensitivity Function');

%% Problem B4: calculate and present freq response of plant and controller
%CT TF for S

C = T2 .* (1 + T3 ./ (1 - T3));
% C_coh = C2 .* (1 + C3 ./ (1 - C3));
% C = T2 ./T1;
C_coh = C2./C1;
% C_coh = 1 - c_coh;


P = (1 - T1) ./ (T1 .* C);
% P_coh = (1 - C1) ./ (C1 .* C_coh);
P_coh = C3./C1 .* (1./(C2./C1));

figure;
subplot(3,1,1);
semilogx(f1,mag2db(abs(C))); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Identified Controller Magnitude');

subplot(3,1,2);
semilogx(f1,rad2deg(angle(C))); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on;
title('Identified Controller Phase');

subplot(3,1,3);
semilogx(f1,C_coh); hold on;
xlim([min(f1),10000]);
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;
title('Identified Controller Coherence');




figure;
subplot(3,1,1);
semilogx(f3,mag2db(abs(P))); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Identified Plant Magnitude');


subplot(3,1,2);
semilogx(f3,rad2deg(angle(P))); hold on;
xlim([min(f1),max(f1)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Identified Plant Magnitude');

subplot(3,1,3);
semilogx(f1,P_coh); hold on;
xlim([min(f1),10000]);
ylim([0, 1]);
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;
title('Identified Plant Coherence');




%% Problem B7: Recreate controller
% %Behaviour identified on OneNote
s = tf('s');

K1 = 1600;       % Proportional gain: moves entire plot up or down
z = 733;        % Zero of 2nd order PD: break frequency
Q = 6;          % Damping of 2nd order PD: 10 
p1 = 5000;      % Pole 1 
p2 = 6000;     % Pole 2 

% Low freq integrator (pole)
I_low = 1 / s;

% 2nd Order PD Filter
PD = K1 * ((s/(2*pi*z))^2 + (s/(Q*(2*pi*z))) + 1);

% High freq integrator (pole)
I_high = 1 / (1 + s/(2*pi*p1));

% Higher freq integrator (pole)
I_higher = 1 / (1.20 + s/(2*pi*p2));  %gain to lower plot

C_an = I_low * PD * I_high * I_higher;



% [T_tf, f_tf] = tfestimate(e, u, window, overlap, ft, Fs);
% [C_tf, f_tf] = mscohere(e, u, window, overlap, ft, Fs)




f = logspace(0,5,1000);
w = 2 * pi * f;

[mag, phase] = bode(C_an, w);
mag_rec = squeeze(mag);
phase_rec = squeeze(phase);

%Plot the magnitude responses
figure;
subplot(2,1,1);
semilogx(f1, mag2db(abs(C))); % Recreated controller in red
hold on;
semilogx(f, 20*log10(mag_rec), 'r--', 'LineWidth', 1.5); % Identified controller in blue dashed line
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Response: Recreated vs Identified Controllers');
legend('Identified Controller', 'Recreated Controller');

% Plot the phase responses
subplot(2,1,2);
semilogx(f1, rad2deg(angle(C))); % Recreated controller in red
hold on;
semilogx(f, phase_rec, 'r--', 'LineWidth', 1.5); % Identified controller in blue dashed line
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Phase Response: Recreated vs Identified Controllers');
legend('Identified Controller', 'Recreated Controller');

%% Tfestimate test

% Define the window and overlap
window_length = round(length(e) / 11);    % 9 chirp signals, so the windows will only cover the chirp signals
window = hann(window_length);         % Use rectwin window
overlap = round(0.8 * window_length); % 50% overlap

%THERE ARE 10 SEGMENTS IN THE SIGNALS (1 STARTER AND 9 CHRIP RESPONSES) So
%this line divides the total input signal
ft = 10000;
% ft = logspace(0,5,10000);

[T_tf, f_tf] = tfestimate(e, u, window, overlap, ft, Fs);
[C_tf, f_tf] = mscohere(e, u, window, overlap, ft, Fs)

[CP_tf, f_tf] = mscohere(u, y, window, overlap, ft, Fs);


% figure;
% %e/r TF and coherence
% subplot(3,1,1);
% semilogx(f_tf,mag2db(abs(T_tf))); hold on;
% xlim([1,max(f_tf)]);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% grid on;
% title('Estimated C');
% 
% subplot(3,1,2);
% semilogx(f_tf,rad2deg(angle(T_tf))); hold on;
% xlim([1,max(f_tf)]);
% xlabel('Frequency (Hz)');
% ylabel('Phase (deg)');
% grid on;
% title('Ectimated C');
% 
% subplot(3,1,3);
% semilogx(f_tf,C_tf); hold on;
% xlim([1,max(f1)]);
% xlabel('Frequency (Hz)');
% ylabel('Coherence');
% grid on;
% title('Coherence estimated C');
% 
% 
% 
% 
% %Prank plot controller
% subplot(3,1,1);
% semilogx(f1,mag2db(abs(C))); hold on;
% xlim([min(f1),max(f1)]);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% grid on;
% title('Identified Controller Magnitude');
% 
% subplot(3,1,2);
% semilogx(f1,rad2deg(angle(C))); hold on;
% xlim([min(f1),max(f1)]);
% xlabel('Frequency (Hz)');
% ylabel('Phase (deg)');
% grid on;
% title('Identified Controller Phase');
% 
% subplot(3,1,3);
% semilogx(f1,C_tf); hold on;
% xlim([min(f1),max(f1)]);
% xlabel('Frequency (Hz)');
% ylabel('Coherence');
% grid on;
% title('Identified Controller Coherence');
% 
% 
% figure;
% subplot(3,1,1);
% semilogx(f3,mag2db(abs(P))); hold on;
% xlim([min(f1),max(f1)]);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% grid on;
% title('Identified Plant Magnitude');
% 
% 
% subplot(3,1,2);
% semilogx(f3,rad2deg(angle(P))); hold on;
% xlim([min(f1),max(f1)]);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% grid on;
% title('Identified Plant Magnitude');
% 
% subplot(3,1,3);
% semilogx(f1,CP_tf); hold on;
% xlim([min(f1),max(f1)]);
% xlabel('Frequency (Hz)');
% ylabel('Coherence');
% grid on;
% title('Identified Plant Coherence');


%% Combined plot
% 
% figure;
% %e/r TF and coherence
% subplot(2,1,1);
% semilogx(f_tf,mag2db(abs(T_tf))); hold on;
% semilogx(f1, mag2db(abs(C)), 'r--');
% xlim([1,max(f_tf)]);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% grid on;
% title('Compared analytical and estimated C');
% legend('Estimated Controller', 'Analytical Controller');
% 
% 
% subplot(2,1,2);
% semilogx(f_tf,rad2deg(angle(T_tf))); hold on;
% semilogx(f1, rad2deg(angle(C)), 'r--');
% xlim([1,max(f_tf)]);
% xlabel('Frequency (Hz)');
% ylabel('Phase (deg)');
% grid on;
% title('Compared analytical and estimated C');
% legend('Estimated Controller', 'Analytical Controller');




%% Try with tfest()

% data = frd(C,2*pi*f1,Ts); %Make FRD Data
% 
% % Use tfest() function to obtain transfer function from data
% % General Configuration: sys = tfest(data,np,nz,iodelay);
% % Important to model delay for controller design
% 
% np = 3; % Tune number of poles
% nz = 2; % Tune number of zeros
% iodelay = 1; % Tune delay 
% sys = tfest(data,np,nz,iodelay);
% Pnump = sys.Numerator;
% Pdenp = sys.Denominator;
% Ptf = tf(Pnump,Pdenp);
% figure;
% bode(Ptf); %Plant Transfer Function from Identification
% title('Plant Transfer Function from Identification');


% Define options for Bode plot
% opts = bodeoptions;
% opts.FreqUnits = 'Hz'; % Change frequency units to Hertz
% % Plot the Bode diagram using bodeplot
% figure;
% bodeplot(G, opts); % Use bodeplot instead of bode
% grid on;
% title('Open Loop Response (Frequency in Hz)');


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



disp(window_length);

disp(C1);





