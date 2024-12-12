PM = 60; % degrees
GM = 6; % dB
MM = 6; % dB

% Assume PID phase is max 50 deg
phase_PID = 50; % deg
phase_system = PM - (phase_PID+180);

% Rule of thumb
G_wc = 1.585; % linear units
kp_pid = 0.33/G_wc;