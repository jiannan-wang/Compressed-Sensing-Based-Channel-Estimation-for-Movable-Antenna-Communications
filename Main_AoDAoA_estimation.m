% Compressed Sensing Based Channel Estimation 
% for Movable Antenna Communications
clear; clc; close all;

%% Parameter Settings
Lt = 3;                 % Number of multipaths at the transmitter
Lr = 3;                 % Number of multipaths at the receiver
lambda = 0.01;          % Carrier wavelength
A = 4*lambda;           % Area size
delta = lambda/5;       % distance of adjecent grids
M = 256;                % Number of T-MA movement positions
N = 256;                % Number of R-MA movement positions
G = 200;                % Angular grid resolution
SNR_dB = 10;            % Signal-to-noise ratio (dB)
SNR_linear = 10^(SNR_dB/10);     
sigma2 = 1;             % Noise variance
P = sigma2*SNR_linear;  % Transmit power

% Generate true AoD/AoA and PRM 
[theta_t, phi_t] = RandomAngle_Generator([0,pi],[0,pi],Lt);
[theta_r, phi_r] = RandomAngle_Generator([0,pi],[0,pi],Lr);
true_Angle = struct('theta_t', theta_t, ...
                    'phi_t', phi_t, ...
                    'theta_r', theta_r, ...
                    'phi_r', phi_r);
PRM = (randn(Lr, Lt) + 1j*randn(Lr, Lt))/sqrt(2);
%% Steering matrix for AoD/AoA estimation
pos = struct('tm_x', (2*rand(M,1)-1)*A/2, ...
             'tm_y', (2*rand(M,1)-1)*A/2, ...
             'rn_x', (2*rand(N,1)-1)*A/2, ...
             'rn_y', (2*rand(N,1)-1)*A/2);
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
A_bar = exp(-1j*2*pi/lambda * (pos.tm_x*theta_pairs.' + pos.tm_y*phi_pairs.'));
B_bar = exp(-1j*2*pi/lambda * (pos.rn_x*theta_pairs.' + pos.rn_y*phi_pairs.'));
%% Sparse Recovery Using OMP Algorithm
% AoD and AoA estimation
Lt_esti = Lt;
Lr_esti = Lr;
[est_theta_t, est_phi_t] = angles_estimator(yt, A_bar, Lt_esti, theta_pairs, phi_pairs);
[est_theta_r, est_phi_r] = angles_estimator(yr, B_bar, Lr_esti, theta_pairs, phi_pairs);
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