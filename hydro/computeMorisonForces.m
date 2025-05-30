function [FHyd, FHyd_history] = computeMorisonForces(U, Udot, time, elemZ, params, FHyd_history, waveVel, waveAcc)
% COMPUTEMORISONFORCES Computes hydrodynamic forces using Morison equation
%
% Inputs:
%   U            - Current displacement vector
%   Udot         - Current velocity vector
%   time         - Current time
%   elemZ        - Z-coordinates of nodes
%   params       - Structure containing simulation parameters
%   FHyd_history - History of hydrodynamic forces (for updating)
%   waveVel      - Pre-computed wave velocities at current time step
%   waveAcc      - Pre-computed wave accelerations at current time step
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
    
    % Get horizontal DOF index
    dof_j = 2*(j-1) + 1;  % Horizontal DOF
    
    % Get wave kinematics from pre-computed arrays
    u = waveVel(dof_j);
    udot = waveAcc(dof_j);
    
    % Get structural velocity at node
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