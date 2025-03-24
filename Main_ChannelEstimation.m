% Compressed Sensing Based Channel Estimation 
% for Movable Antenna Communications
% Author for this code: Jiannan Wang
clear; clc; close all;

%% Parameter Settings
Lt = 3;                 % Number of multipaths at the transmitter
Lr = 3;                 % Number of multipaths at the receiver
lambda = 0.01;          % Carrier wavelength
A = 4*lambda;           % Area size
M = 256;                % Number of T-MA movement positions
N = 256;                % Number of R-MA movement positions
G = 200;                % Angular grid resolution
SNR_dB = 10;            % Signal-to-noise ratio (dB)
SNR_linear = 10^(SNR_dB/10);     
sigma2 = 1;             % Noise variance
P = sigma2*SNR_linear;  % Transmit power

% Generate true AoD/AoA and PRM randomly
[theta_t, phi_t] = RandomAngle_Generator([0,pi],[0,pi],Lt);
[theta_r, phi_r] = RandomAngle_Generator([0,pi],[0,pi],Lr);
true_Angle = struct('theta_t', theta_t, ...
                    'phi_t', phi_t, ...
                    'theta_r', theta_r, ...
                    'phi_r', phi_r);

% PRM = (randn(Lr, Lt) + 1j*randn(Lr, Lt))/sqrt(2);

% assuming that Lr = Lt
eta = 1;  % the ratio of the average power for diagnonal elements to that for non-diasonal elements
diag_var = eta/(eta+1)*Lr;
nondiag_var = 1/((eta+1)*(Lr-1)*Lr);
PRM = zeros(Lr,Lr);
for p=1:Lr
    for q=1:Lr
        if p==q
            PRM(p,p) = sqrt(diag_var/2)*(randn+1j*randn);
        else
            PRM(p,q) = sqrt(nondiag_var/2)*(randn+1j*randn);
        end
    end
end
%% Steering matrix for AoD/AoA estimation
% pos = struct('tm_x', (2*rand(M,1)-1)*A/2, ...
%              'tm_y', (2*rand(M,1)-1)*A/2, ...
%              'rn_x', (2*rand(N,1)-1)*A/2, ...
%              'rn_y', (2*rand(N,1)-1)*A/2);
% square shape
space = 4*A/M;
upper_edge = [-A/2:space:A/2-space;A/2*ones(1,M/4)];  % start point is closed, end point is open
right_edge = [A/2*ones(1,M/4);A/2:-space:-A/2+space];
lower_edge = [A/2:-space:-A/2+space;-A/2*ones(1,M/4)];
left_edge = [-A/2*ones(1,M/4);-A/2:space:A/2-space];

pos.tm_x = [upper_edge(1,:),lower_edge(1,:),right_edge(1,:),left_edge(1,:)].';
pos.tm_y = [upper_edge(2,:),lower_edge(2,:),right_edge(2,:),left_edge(2,:)].';
pos.rn_x = pos.tm_x;
pos.rn_y = pos.tm_y;
figure,scatter(pos.tm_x,pos.tm_y,'filled');
title('square shape MA');grid on;

[A_alphat, Psi_yt, B_alphat, Psi_yr] = Receiver(pos, true_Angle, lambda, P);
zt = sqrt(sigma2/2)*(randn(M,1)+1j*randn(M,1));
yt = conj(Psi_yt*PRM(:)) + zt;  % equ 10

zr = sqrt(sigma2/2)*(randn(N,1)+1j*randn(N,1));
yr = Psi_yr*PRM(:) + zr;  % equ 10
%% Construct Overcomplete Dictionary
theta_grid = linspace(-1, 1, G);
phi_grid = linspace(-1, 1, G);
[theta_mesh, phi_mesh] = meshgrid(theta_grid, phi_grid);
theta_pairs = theta_mesh(:);
phi_pairs = phi_mesh(:);
A_bar = exp(-1j*2*pi/lambda*(pos.tm_x*theta_pairs.'+pos.tm_y*phi_pairs.'));
B_bar = exp(-1j*2*pi/lambda*(pos.rn_x*theta_pairs.'+pos.rn_y*phi_pairs.'));
%% Sparse Recovery Using OMP Algorithm
% AoD and AoA estimation
Lt_esti = Lt;
Lr_esti = Lr;
[est_theta_t, est_phi_t] = angles_estimator(yt, A_bar, Lt_esti, theta_pairs, phi_pairs);
[est_theta_r, est_phi_r] = angles_estimator(yr, B_bar, Lr_esti, theta_pairs, phi_pairs);
esti_Angle = struct('theta_t', est_theta_t, ...
                    'phi_t', est_phi_t, ...
                    'theta_r', est_theta_r, ...
                    'phi_r', est_phi_r);
%% Display Results
figure;
scatter(theta_t, phi_t, 100, 'r', 'filled'); hold on;
scatter(est_theta_t, est_phi_t, 200, 'r', '*');hold on;
scatter(theta_r, phi_r, 100, 'b', 'filled'); hold on;
scatter(est_theta_r, est_phi_r, 200, 'b', 'x');
legend('True AoD', 'Estimated AoD', 'True AoA', 'Estimated AoA');
xlabel('theta'); ylabel('phi');
title('AoD and AoA Estimation Results');
grid on;
%% Estimate PRM with Additional Optimized Positions
K = max(Lt_esti * Lr_esti - Lt_esti - Lr_esti, 0)+1;  % number of additional measurement
if K > 0
    %% Define optimization problem for additional positions
    pos_init = struct('ta_x', (2*rand(K,1)-1)*A/2, ...
                      'ta_y', (2*rand(K,1)-1)*A/2, ...
                      'ra_x', (2*rand(K,1)-1)*A/2, ...
                      'ra_y', (2*rand(K,1)-1)*A/2);
    
    % Objective function handle
    obj_fun = @(x) cond_obj_function(x, K, lambda, true_Angle, P, [Psi_yt; Psi_yr]);
    
    % Optimize using fmincon
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
        'Display', 'off', 'MaxFunctionEvaluations', 10000);
    % Bounds for positions
    lb = -A/2*ones(4*K, 1);
    ub = A/2*ones(4*K, 1);
    x0 = [pos_init.ta_x;pos_init.ta_y;pos_init.ra_x;pos_init.ra_y];
    [pos_opt, ~] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);
    
    % Extract optimized positions
    ta_x_opt = pos_opt(1:K);
    ta_y_opt = pos_opt(K+1:2*K);
    ra_x_opt = pos_opt(2*K+1:3*K);
    ra_y_opt = pos_opt(3*K+1:4*K);
    
    %% Generate additional channel measurements using TRUE channel parameters
    y_additional = zeros(K, 1);
    Psi_additional = zeros(K, Lt_esti*Lr_esti);
    for k = 1:K
        % True steering vectors
        g_tka_true = exp(-1j*2*pi/lambda*(ta_x_opt(k)*true_Angle.theta_t+ta_y_opt(k)*true_Angle.phi_t));
        f_rka_true = exp(-1j*2*pi/lambda*(ra_x_opt(k)*true_Angle.theta_r+ra_y_opt(k)*true_Angle.phi_r));
        % Generate noise
        z_ka = sqrt(sigma2/2) * (randn + 1j*randn);
        y_additional(k) = sqrt(P)*kron(g_tka_true.',f_rka_true')*PRM(:) + z_ka;

        % Estimated steering vectors
        g_tka_esti = exp(-1j*2*pi/lambda*(ta_x_opt(k)*esti_Angle.theta_t+ta_y_opt(k)*esti_Angle.phi_t));
        f_rka_esti = exp(-1j*2*pi/lambda*(ra_x_opt(k)*esti_Angle.theta_r+ra_y_opt(k)*esti_Angle.phi_r));
        Psi_additional(k, :) = sqrt(P)*kron(g_tka_esti.', f_rka_esti');
    end
    
    Psi = [Psi_yt; Psi_yr; Psi_additional];
    rank([Psi_yt;Psi_yr])
    rank(Psi_additional)
    rank(Psi)
    y = [conj(yt); yr; y_additional];

    if rank(Psi)<size(Psi,2)
        disp('Attention! Observed matrix is not full-column rank');
    else
        disp('Good! Observed matrix is full-column rank');
    end

    epsilon_hat = Psi \ y;  % least square estimator
    PRM_hat = reshape(epsilon_hat, Lr_esti, Lt_esti);
    
    mean(abs(epsilon_hat-PRM(:)).^2)
    %% Display PRM estimation results
    figure;
    subplot(1,2,1);
    imagesc(abs(PRM));
    title('True PRM');
    colorbar;
    subplot(1,2,2);
    imagesc(abs(PRM_hat));
    title('Estimated PRM');
    colorbar;
    % PRM
    % PRM_hat
    
else
    disp('No additional measurements needed.');
end
%% Channel Reconstruction
delta = lambda/20;       % distance of adjecent grids
D = (A/delta)^2;         % number of total grids   

fixed_Tx_MA = [-A/2+delta/2,A/2-delta/2];
% fixed_Tx_MA = [0,0];

row = (-A/2+delta/2):delta:(A/2-delta/2);
varied_Rx_MA.x = repmat(row,[sqrt(D),1]);
col = ((A/2-delta/2):-delta:(-A/2+delta/2))';
varied_Rx_MA.y = repmat(col,[1,sqrt(D)]);

H_hat = zeros(sqrt(D),sqrt(D));
for i=1:sqrt(D)
    for j=1:sqrt(D)
        g_t_hat = exp(-1j*2*pi/lambda*(fixed_Tx_MA(1)*esti_Angle.theta_t+fixed_Tx_MA(2)*esti_Angle.phi_t));
        f_r_hat = exp(-1j*2*pi/lambda*(varied_Rx_MA.x(i,j)*esti_Angle.theta_r+varied_Rx_MA.y(i,j)*esti_Angle.phi_r));
        H_hat(i,j) = f_r_hat'*PRM_hat*g_t_hat;
    end
end
figure,imagesc(row/lambda,flip(col)/lambda,abs(H_hat)./max(abs(H_hat(:))));colorbar;
xlabel('Normalized x/{/lambda}'),ylabel('Normalized y/{/lambda}');
title('Tx-MA is fixed while Rx-MA is moving along all the grids');