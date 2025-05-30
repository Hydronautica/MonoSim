function [FHyd, FHyd_history] = computeMorisonForces(U, Udot, time, elemZ, params, FHyd_history)
% COMPUTEMORISONFORCES Computes hydrodynamic forces using Morison equation
%
% Inputs:
%   U            - Current displacement vector
%   Udot         - Current velocity vector
%   time         - Current time
%   elemZ        - Z-coordinates of nodes
%   params       - Structure containing simulation parameters
%   FHyd_history - History of hydrodynamic forces (for updating)
%
% Outputs:
%   FHyd         - Current hydrodynamic force vector
%   FHyd_history - Updated history of hydrodynamic forces

% Extract parameters
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

% Check if eta2 exists in params, otherwise set to 0
if isfield(params, 'eta2')
    eta2 = params.eta2;
else
    eta2 = 0;
end

secondOrder = params.secondOrder;
irregular = params.irregular;
step = params.step;

% Initialize force vector
FHyd = zeros(size(U));

% Loop through each node
for j = 1:length(elemZ)
    z = elemZ(j);
    
    % Skip nodes above water level
    if z > h
        continue;
    end
    
    % Determine diameter based on section
    if z <= sectionCutOff
        D = D_outer1;
    else
        D = D_outer2;
    end
    
    % Calculate wave kinematics at node position
    if irregular
        % First-order irregular wave kinematics
        u = 0;
        udot = 0;
        for m = 1:Nfreq
            k_m = k_m_vec(m);
            omega_m = omegaVec(m);
            A = A_m(m);
            phi = phiRand(m);
            
            % Calculate wave number and frequency components
            u = u + A * omega_m * cosh(k_m*(z+h))/sinh(k_m*h) * cos(-omega_m*time + phi);
            udot = udot - A * omega_m^2 * cosh(k_m*(z+h))/sinh(k_m*h) * sin(-omega_m*time + phi);
        end
        
        % Add second-order components if enabled
        if secondOrder
            % Get second-order velocity and acceleration
            u2 = secondOrderVelocity(time, z, params);
            udot2 = secondOrderAcceleration(time, z, params);
            
            % Add to first-order components
            u = u + u2;
            udot = udot + udot2;
        end
    else
        % Regular wave kinematics (not implemented in this version)
        u = 0;
        udot = 0;
    end
    
    % Get structural velocity at node
    dof_j = 2*(j-1) + 1;  % Horizontal DOF
    Udot_j = Udot(dof_j);
    
    % Calculate relative velocity
    u_rel = u - Udot_j;
    
    % Compute Morison force components
    F_drag = 0.5 * rho_w * Cd * D * Le * abs(u_rel) * u_rel;
    F_inertia = rho_w * Cm * pi * D^2/4 * Le * udot;
    
    % Total force at node
    F_j = F_drag + F_inertia;
    
    % Add to global force vector
    FHyd(dof_j) = F_j;
    
    % Update history if provided
    if nargin > 5 && ~isempty(FHyd_history)
        FHyd_history(dof_j, step) = F_j;
    end
end

end

% Helper function for second-order velocity
function u2 = secondOrderVelocity(time, z, params)
    % Extract parameters
    h = params.h;
    Nfreq = params.Nfreq;
    omegaVec = params.omegaVec;
    A_m = params.A_m;
    phiRand = params.phiRand;
    k_m_vec = params.k_m_vec;
    
    u2 = 0;
    for m = 1:Nfreq
        for n = 1:Nfreq
            % Calculate B+ and B- coefficients
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
                 cosh(k_p*(z+h))/sinh(k_p*h) * ...
                 cos(-omega_p*time + phi_m + phi_n);
            
            % B- term (difference frequency)
            u2 = u2 + 0.5 * A_m(m) * A_m(n) * Bminus * omega_m * ...
                 cosh(k_m*(z+h))/sinh(k_m*h) * ...
                 cos(-omega_m*time + phi_m - phi_n);
        end
    end
end

% Helper function for second-order acceleration
function udot2 = secondOrderAcceleration(time, z, params)
    % Extract parameters
    h = params.h;
    Nfreq = params.Nfreq;
    omegaVec = params.omegaVec;
    A_m = params.A_m;
    phiRand = params.phiRand;
    k_m_vec = params.k_m_vec;
    
    udot2 = 0;
    for m = 1:Nfreq
        for n = 1:Nfreq
            % Calculate B+ and B- coefficients
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
            udot2 = udot2 + 0.5 * A_m(m) * A_m(n) * Bplus * omega_p^2 * ...
                    cosh(k_p*(z+h))/sinh(k_p*h) * ...
                    sin(-omega_p*time + phi_m + phi_n);
            
            % B- term (difference frequency)
            udot2 = udot2 + 0.5 * A_m(m) * A_m(n) * Bminus * omega_m^2 * ...
                    cosh(k_m*(z+h))/sinh(k_m*h) * ...
                    sin(-omega_m*time + phi_m - phi_n);
        end
    end
end