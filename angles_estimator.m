function [est_theta, est_phi] = angles_estimator(y, dict, L, theta_pairs, phi_pairs)
%   This function implements a OMP-based angle estimation pipeline
%
%   Inputs:
%   y            - Received signal vector (column vector)
%   dict         - Dictionary matrix (N x G^2)
%   L            - Number of multipath components to estimate
%   theta_pairs  - Azimuth angle mapping table (G^2 x 1)
%   phi_pairs    - Elevation angle mapping table (G^2 x 1)
%
%   Outputs:
%   est_theta    - Estimated azimuth angles (radians)
%   est_phi      - Estimated elevation angles (radians)
%

    % OMP-based support recovery
    support = OMP(y, dict, L);
    
    % Sparse solution estimation
    sub_dict = dict(:, support);
    x_bar = zeros(size(dict, 2), 1);
    x_bar(support) = sub_dict \ y;  % Matrix left division for min L2-norm solution
    
    % Peak identification and angle mapping
    [~, est_idx] = maxk(abs(x_bar), L);
    est_theta = theta_pairs(est_idx);
    est_phi = phi_pairs(est_idx);
end
