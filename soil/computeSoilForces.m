function [F_py, soilStates, thetaStates] = computeSoilForces(U, V, elemZ, params, soilStates, thetaStates)
% COMPUTESOILFORCES Computes soil reaction forces using p-y and m-theta curves
%
% Inputs:
%   U           - Current displacement vector
%   V           - Current velocity vector
%   elemZ       - Z-coordinates of nodes
%   params      - Structure containing simulation parameters
%   soilStates  - Array of structures for p-y curve state
%   thetaStates - Array of structures for m-theta curve state
%
% Outputs:
%   F_py        - Vector of soil reaction forces
%   soilStates  - Updated array of structures for p-y curve state
%   thetaStates - Updated array of structures for m-theta curve state

% Extract parameters
nNode = params.nNode;
Le = params.Le;
useLinearPY = params.linearPY;

% Initialize soil force vector
F_py = zeros(2*nNode, 1);

% Loop through all nodes
for j = 1:nNode
    z = elemZ(j);
    
    % Only apply soil forces below mudline (z < 0)
    if z < 0
        % Get DOFs for this node
        dof = 2*(j-1) + 1;       % Lateral displacement DOF
        dof_rot = 2*(j-1) + 2;   % Rotational DOF
        
        % Get current displacement and velocity
        y_disp = U(dof);
        y_vel = V(dof);
        depth = -z;  % Positive depth below mudline
        
        % Compute lateral spring force (p-y curve)
        [p_spring, soilStates(j)] = pySimpleHysteretic(y_disp, y_vel, soilStates(j), params, depth);
        
        % Compute lateral damping force
        c_ref = 0;               % N·s/m³, standard for sand
        c_damp = c_ref * depth;  % [N·s/m²]
        f_damp = -c_damp * y_vel;
        
        % Total lateral soil force
        p_total = p_spring + f_damp;
        
        % Apply lateral force to global vector
        F_py(dof) = F_py(dof) - p_total * Le;
        
        % Get current rotation and rotational velocity
        theta = U(dof_rot);
        theta_dot = V(dof_rot);
        
        % Compute rotational spring moment (m-theta curve)
        [m_spring, thetaStates(j)] = mrSimpleHysteretic(theta, theta_dot, thetaStates(j), params, depth);
        
        % Compute rotational damping
        D_pile = params.D_outer_soil;  % pile diameter at this depth
        c_ref = 0;                     % N·s/m³ (translational base)
        c_rot_ref = c_ref * (D_pile/2)^2;  % ≃ N·m·s/rad per m depth
        c_rot = c_rot_ref * depth;        % [N·m·s/rad]
        m_damp = -c_rot * theta_dot;
        
        % Total rotational moment
        m_total = m_spring + m_damp;
        
        % Apply rotational moment to global vector
        F_py(dof_rot) = F_py(dof_rot) - m_total * Le;
    end
end

end