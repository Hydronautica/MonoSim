function F_guy = computeGuyWireForces(U, V, params)
% COMPUTEGUYWIREFORCES Computes forces from guy wires
%
% Inputs:
%   U      - Current displacement vector
%   V      - Current velocity vector
%   params - Structure containing simulation parameters
%
% Outputs:
%   F_guy  - Vector of guy wire forces

% Extract parameters
nDOF = length(U);
dof_guy = params.dof_guy;
z_nodeG1 = params.z_nodeG1;

% Right guy wire parameters
anchorR_x = params.anchorR_x;
anchorR_z = params.anchorR_z;
k_guyR = params.k_guyR;
L0_guyR = params.L0_guyR;

% Left guy wire parameters
anchorL_x = params.anchorL_x;
anchorL_z = params.anchorL_z;
k_guyL = params.k_guyL;
L0_guyL = params.L0_guyL;

% Initialize guy wire force vector
F_guy = zeros(nDOF, 1);

% Current lateral displacement and velocity at the attachment node
x_node = U(dof_guy);
x_dot = V(dof_guy);

% Damping coefficients
c_guyR = 1e6;  % Damping coefficient for right guy wire
c_guyL = 1e6;  % Damping coefficient for left guy wire

% Right guy wire
dxR = anchorR_x - x_node;
dzR = anchorR_z - z_nodeG1;
L_R = sqrt(dxR^2 + dzR^2);

if L_R > L0_guyR
    % Smoothing parameter for tanh function
    eps_smooth = 0.01;
    
    % Calculate elongation
    delta_L_R = max(0, L_R - L0_guyR);
    
    % Calculate tension using smoothed function
    T_R = k_guyR * tanh(delta_L_R / eps_smooth) * eps_smooth;
    
    % Direction cosine in x-direction
    dirR = dxR / L_R;
    
    % Damping force (opposes velocity)
    F_damp_R = -c_guyR * x_dot * dirR;
    
    % Total horizontal force
    Fx_R = T_R * dirR + F_damp_R;
    
    % Apply to global force vector
    F_guy(dof_guy) = F_guy(dof_guy) - Fx_R;
end

% Left guy wire
dxL = anchorL_x - x_node;
dzL = anchorL_z - z_nodeG1;
L_L = sqrt(dxL^2 + dzL^2);

if L_L > L0_guyL
    % Smoothing parameter for tanh function
    eps_smooth = 0.01;
    
    % Calculate elongation
    delta_L_L = max(0, L_L - L0_guyL);
    
    % Calculate tension using smoothed function
    T_L = k_guyL * tanh(delta_L_L / eps_smooth) * eps_smooth;
    
    % Direction cosine in x-direction
    dirL = dxL / L_L;
    
    % Damping force (opposes velocity)
    F_damp_L = -c_guyL * x_dot * dirL;
    
    % Total horizontal force
    Fx_L = T_L * dirL + F_damp_L;
    
    % Apply to global force vector
    F_guy(dof_guy) = F_guy(dof_guy) - Fx_L;
end

end