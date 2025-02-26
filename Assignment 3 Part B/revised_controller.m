function C_updated = revised_controller(C_old,G)
    % Goal of update
    f_considerd = 10;
    w_considered = 2 * pi * f_considerd;
    
    % Construct signals
    t = linspace(0,10,100);

    % Input signal r: pure sinusoid
    mag_r = 10e-6;            % micrometers   10e-6 meters
    f_r   = f_considerd;               % Hz
    w_r   = f_r * 2 * pi;     % rad/s
    
    r = mag_r * sin(w_r * t); % time domain
    R = mag_r * 1j * w_r;     % Laplace domain: pure sinusoid
    
    
    % Process disturbance d: pure sinusoid
    mag_d = 20;               % V
    f_d   = f_considerd;               % Hz
    w_d   = f_d * 2 * pi;     % rad/s
    
    d = mag_d * sin(w_d * t); % time domain
    D = mag_d * 1j * w_d;     % Laplace domain: pure sinusoid
    
    % Output disturbance n: pure sinusoid
    mag_n = 10e-9;            % micrometers    10e-9 meters
    f_n   = f_considerd;               % Hz
    w_n = f_n * 2 * pi;       % rad/s
    
    n = mag_n * sin(w_n * t); % time domain
    N = mag_n * 1j * w_n;     % Laplace domain: pure sinusoid

    % Construct system transfer functions
    L = C_old * G;               % Open loop transfer function
    
    S = 1 / (1 + L);         % Sensitivity function
    
    PS = G * S;              % Process * sensitivity function
    
    T = L / (1 + L);         % Complimentary sensitivity function
    
    %% Compute error
    C_response = abs(evalfr(C_old, 1j * 20*pi*2));
    G_response = evalfr(G, 1j * w_considered);
    
    S_response = evalfr(S, 1j * w_considered);
    PS_response = evalfr(PS, 1j * w_considered);
    T_response = evalfr(T, 1j * w_considered);
    
    r_contribution = abs(S_response * mag_r);
    d_contribution = abs(PS_response * mag_d);
    n_contribution = abs(T_response * mag_n);
    
    e_real = sqrt(r_contribution^2 + d_contribution^2 + n_contribution^2);  

    % Anti notch implementation
    
    f_anti = 10;                    % Activation frequency for anti notch action
    w_anti = f_anti * 2 * pi;       % Activation frequency in rad/s
    % Q2 / Q1 = extra gain, 
    new_gain = abs((5 * (1 + C_response * G_response ) - 1) * (1 / G_response));
    Q_2_den_anti = 34.195;
    Q_1_num_anti = new_gain / Q_2_den_anti;
    % purpously built Q1 like this to keep is really thin
    
    
    num_anti = [1, 1/(Q_1_num_anti*w_anti), w_anti^2 ];
    den_anti = [1, 1/(Q_2_den_anti*w_anti), w_anti^2 ];
    
    anti_notch = tf(num_anti, den_anti);
    
    %construct controller
    C_updated = C_old * anti_notch;
end