function [theta, phi] = RandomAngle_Generator(elevation_bound, azimuth_bound, num_paths)
%   This function generates virtual angles (theta, phi) through nonlinear 
%   transformations of uniformly distributed random angles within specified ranges.
%
%   INPUTS:
%   - elevation_bound: [1x2] vector specifying [min, max] elevation angle range (radians)
%   - azimuth_bound:   [1x2] vector specifying [min, max] azimuth angle range (radians)
%   - num_paths:       Number of propagation paths to generate
%
%   OUTPUTS:
%   - theta:          [num_paths x 1] transformed elevation angles
%   - phi:            [num_paths x 1] transformed azimuth angles

    theta_tilde = (elevation_bound(2)-elevation_bound(1))*rand(num_paths,1) + elevation_bound(1);     
    phi_tilde = (azimuth_bound(2)-azimuth_bound(1))*rand(num_paths,1) + azimuth_bound(1);    
    
    % virtual AoDs/AoAs
    theta = sin(theta_tilde) .* cos(phi_tilde);
    phi = cos(theta_tilde);
end