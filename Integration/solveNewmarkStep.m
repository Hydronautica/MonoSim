function [U_next, V_next, Aacc_next] = solveNewmarkStep(U, V, Aacc, F, K, M, C, free, fixed, a0c, a1c, a2c, a3c, a4c, a5c, Keff, dt, gammaN)
% SOLVENEWMARKSTEP Solves one step of Newmark-beta integration
%
% Inputs:
%   U      - Current displacement vector
%   V      - Current velocity vector
%   Aacc   - Current acceleration vector
%   F      - Current external force vector
%   K      - Global stiffness matrix
%   M      - Global mass matrix
%   C      - Global damping matrix
%   free   - Indices of free DOFs
%   fixed  - Indices of fixed DOFs
%   a0c    - Newmark constant
%   a1c    - Newmark constant
%   a2c    - Newmark constant
%   a3c    - Newmark constant
%   a4c    - Newmark constant
%   a5c    - Newmark constant
%   Keff   - Effective stiffness matrix
%   dt     - Time step
%   gammaN - Newmark gamma parameter
%
% Outputs:
%   U_next    - Displacement vector at next time step
%   V_next    - Velocity vector at next time step
%   Aacc_next - Acceleration vector at next time step

% Initialize next step vectors
U_next = U;
V_next = V;
Aacc_next = Aacc;

% Extract free DOF components of the force vector
Ff = F(free);

% Compute effective force vector
Feff = Ff + M(free, free) * (a0c*U(free) + a2c*V(free) + a3c*Aacc(free)) + ...
       C(free, free) * (a1c*U(free) + a4c*V(free) + a5c*Aacc(free));

% Solve for displacement at next time step
U_next(free) = Keff \ Feff;

% Update acceleration at next time step
Aacc_next(free) = a0c*(U_next(free) - U(free)) - a2c*V(free) - a3c*Aacc(free);

% Update velocity at next time step
V_next(free) = V(free) + dt*((1-gammaN)*Aacc(free) + gammaN*Aacc_next(free));

% Fixed DOFs remain at zero
U_next(fixed) = 0;
V_next(fixed) = 0;
Aacc_next(fixed) = 0;

end