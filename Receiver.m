function [A, Psi_yt, B, Psi_yr] = Receiver(pos, Angle, lambda, P)
    % Observed matrix for AoD estimation
    A = exp(-1j*2*pi/lambda*(pos.tm_x*Angle.theta_t.'+ ...
                                 pos.tm_y*Angle.phi_t.'));  % M x Lt_esti
    f_r1 = exp(-1j * 2*pi/lambda * (pos.rn_x(1)*Angle.theta_r + ...
                                    pos.rn_y(1)*Angle.phi_r));
    Psi_yt = sqrt(P) * kron(conj(A), f_r1');  % equ 10
    
    % Observed matrix for AoA estimation
    B = exp(-1j * 2*pi/lambda * (pos.rn_x * Angle.theta_r.' + ...
                                     pos.rn_y * Angle.phi_r.'));
    g_tM = exp(-1j * 2*pi/lambda * (pos.tm_x(end)*Angle.theta_t + ...
                                    pos.tm_y(end)*Angle.phi_t));  % N x Lr_esti
    Psi_yr = sqrt(P) * kron(g_tM.', B);  % equ 10

    
end