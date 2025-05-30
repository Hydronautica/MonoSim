function Bplus = computeBplus(m, n, k_m_vec, omegaVec, h)
% COMPUTEBPLUS Computes B+ coefficient for second-order wave theory
%
% Inputs:
%   m        - First frequency index
%   n        - Second frequency index
%   k_m_vec  - Wave number vector
%   omegaVec - Angular frequency vector
%   h        - Water depth
%
% Outputs:
%   Bplus    - B+ coefficient

% Extract wave numbers and frequencies
k_m = k_m_vec(m);
k_n = k_m_vec(n);
omega_m = omegaVec(m);
omega_n = omegaVec(n);

% Sum of wave numbers and frequencies
k_p = k_m + k_n;
omega_p = omega_m + omega_n;

% Calculate D+ denominator
D_plus = (omega_p^2 - 9.81*k_p*tanh(k_p*h));

% Calculate B+ coefficient
Bplus = (omega_m*omega_n/(9.81)) * ...
        (1/(k_m*k_n)) * ...
        ((k_m*k_n)/(k_p) - 0.5*(k_m^2 + k_n^2)) * ...
        (1/D_plus) * ...
        (omega_p^2 + 9.81*k_p*tanh(k_p*h));

end