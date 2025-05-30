function [waveVel, waveAcc, eta2_history] = precomputeWaveKinematics(time, elemZ, params)
% PRECOMPUTEWAVEKINEMATICS Pre-computes wave kinematics for all time steps
%
% Inputs:
%   time   - Time vector for simulation
%   elemZ  - Z-coordinates of nodes
%   params - Structure containing simulation parameters
%
% Outputs:
%   waveVel     - Pre-computed wave velocities [nDOF x nTimeSteps]
%   waveAcc     - Pre-computed wave accelerations [nDOF x nTimeSteps]
%   eta2_history - Second-order surface elevation history

% Extract parameters
g = 9.81;
rho_w = params.rho_w;
Cd = params.Cd;
Cm = params.Cm;
D_outer1 = params.D_outer1;
D_outer2 = params.D_outer2;
sectionCutOff = params.sectionCutOff;
h = params.h;
Le = params.Le;
Nfreq = params.Nfreq;
omegaVec = params.omegaVec;
A_m = params.A_m;
phiRand = params.phiRand;
k_m_vec = params.k_m_vec;
secondOrder = params.secondOrder;
irregular = params.irregular;
Hs = params.H_wave;
Tp = params.T_wave;

% Add thetaVec if not present (direction vector, assume 0 for 1D case)
if ~isfield(params, 'thetaVec')
    thetaVec = zeros(size(omegaVec));
else
    thetaVec = params.thetaVec;
end

% Initialize arrays
nTimeSteps = length(time);
nNodes = length(elemZ);
waveVel = zeros(2*nNodes, nTimeSteps);
waveAcc = zeros(2*nNodes, nTimeSteps);
eta2_history = zeros(1, nTimeSteps);

% First-order surface elevation
eta1 = zeros(1, nTimeSteps);
for m = 1:Nfreq
    for t = 1:nTimeSteps
        eta1(t) = eta1(t) + A_m(m) * cos(omegaVec(m)*time(t) + phiRand(m));
    end
end

% Calculate second-order surface elevation if needed
if secondOrder && irregular
    for t = 1:nTimeSteps
        eta2_history(t) = secondOrderElevation(time(t), params);
    end
end

% Loop through each node
for j = 1:nNodes
    z_b = elemZ(j);
    
    % Skip nodes outside water column
    if z_b < 0 || z_b > h
        continue;
    end
    
    z_w = z_b - h;  % Depth relative to water surface (negative below surface)
    
    % Loop through each time step
    for t = 1:nTimeSteps
        current_time = time(t);
        
        % Calculate wave kinematics at node position
        if ~irregular
            % Regular wave kinematics
            omega0 = 2*pi/Tp;
            k0 = fzero(@(k) omega0^2 - g*k*tanh(k*h), omega0^2/g);
            a = Hs/2;
            C1 = cosh(k0*(z_w + h))/sinh(k0*h);
            u = a*omega0 * C1 * cos(omega0*current_time);
            udot = -a*omega0^2 * C1 * sin(omega0*current_time);
        else
            % First-order irregular wave kinematics
            u = 0;
            udot = 0;
            for m = 1:Nfreq
                k_m = k_m_vec(m);
                omega_m = omegaVec(m);
                A = A_m(m);
                phi = phiRand(m);
                
                % Calculate wave number and frequency components
                decay = cosh(k_m*(z_w + h))/sinh(k_m*h);
                u = u + A * omega_m * decay * cos(omega_m*current_time + phi);
                udot = udot - A * omega_m^2 * decay * sin(omega_m*current_time + phi);
            end
            
            % Add second-order components if enabled
            if secondOrder
                try
                    % Use external functions for second-order components if available
                    [u2x, ~, ~] = secondOrderVelocity(A_m, phiRand, omegaVec, k_m_vec, ...
                                                    thetaVec, h, 0, z_b, current_time);
                    [a2x, ~, ~] = secondOrderAcceleration(A_m, phiRand, omegaVec, k_m_vec, ...
                                                    thetaVec, h, 0, z_b, current_time);
                    u = u + u2x;
                    udot = udot + a2x;
                catch
                    % Fallback to inline calculation if external functions not available
                    u2 = 0;
                    udot2 = 0;
                    
                    for m = 1:Nfreq
                        for n = 1:Nfreq
                            % Calculate B+ and B- coefficients with epsilon to avoid division by zero
                            Bplus = computeBplus(m, n, k_m_vec, omegaVec, h);
                            Bminus = computeBminus(m, n, k_m_vec, omegaVec, h);
                            
                            % Sum and difference of frequencies and wave numbers
                            omega_p = omegaVec(m) + omegaVec(n);
                            omega_m = omegaVec(m) - omegaVec(n);
                            k_p = k_m_vec(m) + k_m_vec(n);
                            k_m = abs(k_m_vec(m) - k_m_vec(n));
                            
                            % Phase angles
                            phi_m = phiRand(m);
                            phi_n = phiRand(n);
                            
                            % B+ term (sum frequency)
                            u2 = u2 + 0.5 * A_m(m) * A_m(n) * Bplus * omega_p * ...
                                 cosh(k_p*(z_w + h))/sinh(k_p*h) * ...
                                 cos(omega_p*current_time + phi_m + phi_n);
                            
                            % B- term (difference frequency)
                            u2 = u2 + 0.5 * A_m(m) * A_m(n) * Bminus * omega_m * ...
                                 cosh(k_m*(z_w + h))/sinh(k_m*h) * ...
                                 cos(omega_m*current_time + phi_m - phi_n);
                            
                            % B+ term (sum frequency) - acceleration
                            udot2 = udot2 - 0.5 * A_m(m) * A_m(n) * Bplus * omega_p^2 * ...
                                    cosh(k_p*(z_w + h))/sinh(k_p*h) * ...
                                    sin(omega_p*current_time + phi_m + phi_n);
                            
                            % B- term (difference frequency) - acceleration
                            udot2 = udot2 - 0.5 * A_m(m) * A_m(n) * Bminus * omega_m^2 * ...
                                    cosh(k_m*(z_w + h))/sinh(k_m*h) * ...
                                    sin(omega_m*current_time + phi_m - phi_n);
                        end
                    end
                    
                    % Add to first-order components
                    u = u + u2;
                    udot = udot + udot2;
                end
            end
        end
        
        % Store in arrays (only horizontal DOF)
        dof_j = 2*(j-1) + 1;  % Horizontal DOF
        waveVel(dof_j, t) = u;
        waveAcc(dof_j, t) = udot;
    end
end

end