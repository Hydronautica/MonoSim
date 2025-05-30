function Bminus = computeBminus(m, n, k_m_vec, omegaVec, h)
% COMPUTEBMINUS Computes B- coefficient for second-order wave theory
%
% Inputs:
%   m        - First frequency index
%   n        - Second frequency index
%   k_m_vec  - Wave number vector
%   omegaVec - Angular frequency vector
%   h        - Water depth
%
% Outputs:
%   Bminus   - B- coefficient

% Extract wave numbers and frequencies
k_m = k_m_vec(m);
k_n = k_m_vec(n);
omega_m = omegaVec(m);
omega_n = omegaVec(n);

% Difference of wave numbers and frequencies
k_d = abs(k_m - k_n);
omega_d = abs(omega_m - omega_n);

% Calculate D- denominator
D_minus = (omega_d^2 - 9.81*k_d*tanh(k_d*h));

% Calculate B- coefficient
Bminus = (omega_m*omega_n/(9.81)) * ...
         (1/(k_m*k_n)) * ...
         ((k_m*k_n)/(k_d) + 0.5*(k_m^2 + k_n^2)) * ...
         (1/D_minus) * ...
         (omega_d^2 + 9.81*k_d*tanh(k_d*h));

end