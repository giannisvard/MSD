function C_real = controller_assignment_2() 
    s = tf('s');
    % Shapeit parameters
    Kp = 1.86;
    fi = 230; %Hz
    fd =  137 ; %Hz
    ft = 90000; %Hz
    % Convert to rad/s
    wi = fi*2*pi; %Integrator zero, rad/s
    wd = fd*2*pi; %Differentiator pole, rad/s
    wt = ft*2*pi; %(Tamed)Differentiator zero, rad/s
    %Create PID transfer function
    C_PID_real = Kp*(1+wi/s)*((s/wd+1)/(s/wt + 1));  %proportional, integral, lead/lag

    % Regular notch filter implementation
    f1 = 738;           % Notch frequency (Hz)
    zeta_1 = 0.01;      % Damping for numerator (sharp)
    zeta_2 = 0.9;       % Damping for denominator (wider)
    % Convert to rad/s
    w1 = f1*2*pi; %rad/s
    % Create transfer function
    num_notch = [1, 2*zeta_1*w1, w1^2];  % Numerator coefficients
    den_notch = [1, 2*zeta_2*w1, w1^2];  % Denominator coefficients
    notch = tf(num_notch, den_notch);

    % Skewed notch implementation
    f21 = 1009;          % Frequency for numerator (Hz)
    zeta_21 = 0.01;      % Damping ratio for numerator
    f22 = 971.7;         % Frequency for denominator (Hz)
    zeta_22 = 0.005;     % Damping ratio for denominator
    % Convert to rad/s
    w21 = f21*2*pi; %Notch zero, rad/s
    w22 = f22*2*pi; %Notch pole, rad/s
    % Create transfer function
    num_skewed = [1, 2*zeta_21*w21, w21^2];  % Numerator coefficients
    den_skewed = [1, 2*zeta_22*w22, w22^2];  % Denominator coefficients
    skewed_notch = tf(num_skewed, den_skewed);


    %construct controller
    C_real = C_PID_real * notch * skewed_notch;
end

