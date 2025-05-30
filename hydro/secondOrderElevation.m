function eta2 = secondOrderElevation(time, params)
% SECONDORDERELEVATION Computes second-order surface elevation
%
% Inputs:
%   time   - Current time
%   params - Structure containing simulation parameters
%
% Outputs:
%   eta2   - Second-order surface elevation component

% Extract parameters
Nfreq = params.Nfreq;
omegaVec = params.omegaVec;
A_m = params.A_m;
phiRand = params.phiRand;
k_m_vec = params.k_m_vec;
h = params.h;

% Initialize second-order elevation
eta2 = 0;

% Double summation for second-order components
for m = 1:Nfreq
    for n = 1:Nfreq
        % Calculate B+ and B- coefficients
        Bplus = computeBplus(m, n, k_m_vec, omegaVec, h);
        Bminus = computeBminus(m, n, k_m_vec, omegaVec, h);
        
        % Sum and difference of frequencies
        omega_p = omegaVec(m) + omegaVec(n);
        omega_m = omegaVec(m) - omegaVec(n);
        
        % Phase angles
        phi_m = phiRand(m);
        phi_n = phiRand(n);
        
        % B+ term (sum frequency)
        eta2 = eta2 + 0.5 * A_m(m) * A_m(n) * Bplus * ...
               cos(-omega_p*time + phi_m + phi_n);
        
        % B- term (difference frequency)
        eta2 = eta2 + 0.5 * A_m(m) * A_m(n) * Bminus * ...
               cos(-omega_m*time + phi_m - phi_n);
    end
end

end