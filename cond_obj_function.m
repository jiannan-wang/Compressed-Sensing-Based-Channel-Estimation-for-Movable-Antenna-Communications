function [cond_value, Psi, Psi_additional] = cond_obj_function(pos, K, lambda, Angle, P, Psi)
    % Extract positions from optimization variable
    ta_x = pos(1:K);
    ta_y = pos(K+1:2*K);
    ra_x = pos(2*K+1:3*K);
    ra_y = pos(3*K+1:4*K);
    
    Psi_additional = [];
    % Add new measurements based on current positions
    for k = 1:K
        g_tka = exp(-1j*2*pi/lambda * (ta_x(k)*Angle.theta_t + ta_y(k)*Angle.phi_t));
        f_rka = exp(-1j*2*pi/lambda * (ra_x(k)*Angle.theta_r + ra_y(k)*Angle.phi_r));
        new_row = sqrt(P) * kron(g_tka.', f_rka');
        Psi_additional = [Psi_additional; new_row];
    end
    
    Psi = [Psi; Psi_additional];

    % Compute condition number
    cond_value = cond(Psi);
end