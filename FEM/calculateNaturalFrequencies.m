function [freq_hz, alpha_ray, beta_ray] = calculateNaturalFrequencies(K, M, free)
% CALCULATENATURALFREQUENCIES Calculates natural frequencies and recommended damping parameters
%
% Inputs:
%   K     - Global stiffness matrix
%   M     - Global mass matrix
%   free  - Indices of free DOFs
%
% Outputs:
%   freq_hz   - Natural frequencies in Hz
%   alpha_ray - Recommended Rayleigh damping alpha parameter
%   beta_ray  - Recommended Rayleigh damping beta parameter

% Use reduced matrices to avoid singularity
K_red = K(free, free);
M_red = M(free, free);

% Number of modes to extract
eigModes = 6;

% Solve eigenvalue problem (shift-invert around zero to get smallest nonzero frequencies)
[~, Omega2_red] = eigs(K_red, M_red, eigModes, 'smallestabs');

% Extract natural frequencies
omega_n = sqrt(diag(Omega2_red));  % rad/s
freq_hz = omega_n / (2*pi);        % Hz

% Display results
fprintf('First natural frequency: %.3f Hz\n', freq_hz(1));
fprintf('Second natural frequency: %.3f Hz\n', freq_hz(2));

% Calculate recommended Rayleigh damping parameters for 3% damping ratio
omega1 = 2*pi*freq_hz(1);    % First mode in rad/s
omega2 = 2*pi*freq_hz(2);    % Second mode in rad/s
zeta = 0.03;                  % Target damping ratio (3%)

% Calculate Rayleigh damping parameters
alpha_ray = 2*zeta*omega1*omega2/(omega1+omega2);
beta_ray = 2*zeta/(omega1+omega2);

% Display recommended damping parameters
fprintf('Recommended Alpha: %.3f\n', alpha_ray);
fprintf('Recommended Beta: %.3f\n', beta_ray);

end