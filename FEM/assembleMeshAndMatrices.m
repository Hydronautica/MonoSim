function [K, M, C, elemZ, nNode, nDOF, free, fixed] = assembleMeshAndMatrices(params)
% ASSEMBLEMESHANDMATRICES Creates the mesh and assembles global stiffness and mass matrices
%
% Inputs:
%   params - Structure containing simulation parameters
%
% Outputs:
%   K      - Global stiffness matrix
%   M      - Global mass matrix
%   C      - Global damping matrix
%   elemZ  - Z-coordinates of nodes
%   nNode  - Number of nodes
%   nDOF   - Number of degrees of freedom
%   free   - Indices of free DOFs
%   fixed  - Indices of fixed DOFs

% Extract parameters
nElem = params.nElem;
L_pile = params.L_pile;
L_above = params.L_above;
h = params.h;
L_total = params.L_total;
sectionCutOff = params.sectionCutOff;
D_outer1 = params.D_outer1;
thick1 = params.thick1;
E1 = params.E1;
rho1 = params.rho1;
D_outer2 = params.D_outer2;
thick2 = params.thick2;
E2 = params.E2;
rho2 = params.rho2;
D_outer_soil = params.D_outer_soil;
thick_soil = params.thick_soil;
E_soil = params.E_soil;
rho_soil = params.rho_soil;
alpha_ray = params.alpha_ray;
beta_ray = params.beta_ray;
m_hub = params.m_hub;
I_hub = params.I_hub;

% Create mesh
nNode = nElem + 1;
Le = L_total / nElem;
elemZ = linspace(-L_pile, L_above + h, nNode);  % from -L_pile (soil) to +L_above (tip)

% Initialize global matrices
nDOF = 2 * nNode;
K = zeros(nDOF);
M = zeros(nDOF);

% Assemble element matrices
for e = 1:nElem
    z1 = elemZ(e);
    z2 = elemZ(e+1);
    zmid = 0.5 * (z1 + z2);

    if zmid < 0
        % Below mudline (embedded in soil)
        D_outer = D_outer_soil;
        thick   = thick_soil;
        E       = E_soil;
        rho     = rho_soil;
    elseif zmid <= sectionCutOff
        % Above-ground Section 1
        D_outer = D_outer1;
        thick   = thick1;
        E       = E1;
        rho     = rho1;
    else
        % Above-ground Section 2
        D_outer = D_outer2;
        thick   = thick2;
        E       = E2;
        rho     = rho2;
    end

    D_inner = D_outer - 2*thick;
    A_cs = pi/4 * (D_outer^2 - D_inner^2);
    Iz = pi/64 * (D_outer^4 - D_inner^4);

    % Element stiffness matrix
    Ke = E*Iz/Le^3 * [12,6*Le,-12,6*Le;
                      6*Le,4*Le^2,-6*Le,2*Le^2;
                      -12,-6*Le,12,-6*Le;
                      6*Le,2*Le^2,-6*Le,4*Le^2];
    
    % Element mass matrix
    Me = rho*A_cs*Le/420 * [156,22*Le,54,-13*Le;
                           22*Le,4*Le^2,13*Le,-3*Le^2;
                           54,13*Le,156,-22*Le;
                           -13*Le,-3*Le^2,-22*Le,4*Le^2];

    % Assemble into global matrices
    dofs = [2*(e-1)+1, 2*(e-1)+2, 2*e+1, 2*e+2];
    K(dofs, dofs) = K(dofs, dofs) + Ke;
    M(dofs, dofs) = M(dofs, dofs) + Me;
end

% Add tip mass and rotational inertia
topDOF = 2*(nNode-1) + 1;
M(topDOF, topDOF) = M(topDOF, topDOF) + m_hub;
topRotDOF = topDOF + 1;    % rotational DOF of the top node
M(topRotDOF, topRotDOF) = M(topRotDOF, topRotDOF) + I_hub;

% Compute damping matrix using Rayleigh damping
C = alpha_ray*M + beta_ray*K;

% Define boundary conditions
fixed = [1, 2];  % Fixed DOFs at the base (displacement and rotation)
all = 1:nDOF;
free = setdiff(all, fixed);  % Free DOFs

end