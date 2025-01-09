s = tf('s');

% Define the plant
P = 1 / (s^2 + 3*s + 2);
figure;
bode(P);

% Design a lead compensator to improve phase margin
C_lead = (1 + 0.1*s) / (1 + 0.01*s);

% Add a low-pass filter to suppress high-frequency noise
C_lpf = 1 / (1 + s/100);

% Combine components to form the controller
C = C_lead * C_lpf;

% Open-loop transfer function
L = C * P;

% Plot Bode diagram
figure;
bode(L);
grid on;
title('Open-Loop Transfer Function');

% Analyze margins
margin(L);