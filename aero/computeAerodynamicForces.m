function [FAero, FAero_history] = computeAerodynamicForces(time, elemZ, params, FAero_history)
% COMPUTEAERODYNAMICFORCES Computes aerodynamic forces on the structure
%
% Inputs:
%   time          - Current time
%   elemZ         - Z-coordinates of nodes
%   params        - Structure containing simulation parameters
%   FAero_history - History of aerodynamic forces (for updating)
%
% Outputs:
%   FAero         - Current aerodynamic force vector
%   FAero_history - Updated history of aerodynamic forces

% Extract parameters
rho_a = params.rho_a;
Cd_a = params.Cd_a;
D_outer2 = params.D_outer2;
Le = params.Le;
h = params.h;
U_ref = params.U_ref;
z_ref = params.z_ref;
alpha_wind = params.alpha_wind;
step = params.step;

% Initialize force vector
FAero = zeros(length(elemZ)*2, 1);

% Loop through each node
for j = 1:length(elemZ)
    z = elemZ(j);
    
    % Skip nodes below water level
    if z <= h
        continue;
    end
    
    % Calculate wind velocity at node height using power law profile
    if z < z_ref
        V_z = U_ref * (z/z_ref)^alpha_wind;
    else
        V_z = U_ref;
    end
    
    % Calculate aerodynamic force
    F_aero_j = 0.5 * rho_a * Cd_a * D_outer2 * Le * V_z^2;
    
    % Add to global force vector (horizontal DOF)
    dof_j = 2*(j-1) + 1;
    FAero(dof_j) = F_aero_j;
    
    % Update history if provided
    if nargin > 3 && ~isempty(FAero_history)
        %FAero_history(dof_j, step) = F_aero_j;
    end
end

end