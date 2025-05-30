function params = initializeParameters()
% INITIALIZEPARAMETERS Sets up all simulation parameters and flags
%
% Returns:
%   params - Structure containing all simulation parameters

%% Flags
params.video = false;
params.guywires = false;
params.secondOrder = true;
params.linearPY = true;
params.irregular = true;

%% Main Parameters
params.g = 9.81;            % gravity [m/s^2]
params.ramp_duration = 20;  % Duration of force linear ramp in seconds
params.dt = 0.02;           % time step [s]
params.t_total = 100;       % total simulation time [s]
params.nSteps = round(params.t_total / params.dt);

%% Section & Tip Properties
params.m_hub = 1;           % tip mass [kg]
params.I_hub = 1e1;         % tip rotational inertia [kg·m^2]
params.F_hub = 1e5;         % constant horizontal tip load [N]
params.F_hub_oscill = 0;    % oscillatory tip load amplitude [N]
params.F_hub_freq = 0.8;    % oscillatory tip load frequency [rad/s]

%% Wave Kinematics
params.gammaJ = 3.3;        % JONSWAP peak enhancement factor
params.H_wave = 4;       % wave height [m]
params.T_wave = 8;         % wave period [s]
params.h = 23;              % [m], water depth
params.waterDepth = params.h; % alias for water depth [m]
params.rho_w = 1025;    % water density [kg/m^3]
params.Cd = 1.0;           % Morison drag coefficient
params.Cm = 2.0;           % Morison inertia coefficient
params.Nfreq = 16;         % number of spectral components
params.T_min = 4;           % Low Period Cutoff
params.T_max = 30;          % High Period Cutoff

%% Aerodynamic (Air) Properties
params.rho_a = 1.225;     % air density [kg/m^3]
params.Cd_a = 1.2;        % aerodynamic drag coefficient
params.U_ref = 0;           % wind speed at reference height [m/s]
params.z_ref = 10;          % reference height for wind profile [m]
params.alpha_wind = 0.14;        % wind power-law exponent

%% Section Properties
params.D_outer1 = 7;        % Section 1 outer diameter [m]
params.thick1 = 0.07;       % Section 1 wall thickness [m]
params.E1 = 210e9;          % Section 1 Young's modulus [Pa]
params.rho1 = 7850;         % Section 1 density [kg/m^3]
params.D_outer2 = 7;        % Section 2 outer diameter [m]
params.thick2 = 0.07;       % Section 2 wall thickness [m]
params.E2 = 210e9;          % Section 2 Young's modulus [Pa]
params.rho2 = 7850;         % Section 2 density [kg/m^3]
params.D_outer_soil = 7;    % Outer diameter in soil [m]
params.thick_soil = 0.07;   % Wall thickness in soil [m]
params.E_soil = 210e9;      % Young's modulus in soil
params.rho_soil = 7850;     % Density in soil

%% Damping & Integration
params.alpha_ray = 0.077;   % Rayleigh damping alpha
params.beta_ray = 0.006;    % Rayleigh damping beta
params.gammaN = 0.5;        % Newmark-beta gamma
params.betaN = 0.25;        % Newmark-beta beta

%% Soil Properties
params.k_soil = 1e6;        % Lateral soil stiffness [N/m³]
params.k_soil_rot = 1e7;    % Rotational soil stiffness [N·m/rad/m²]

%% Geometry & Mesh
params.L_pile = 30;         % Embedment depth [m]
params.L_above = 120;       % Original above-ground length
params.TowerHeight = 90;
params.sectionCutOff = params.L_above - params.TowerHeight; % [m], section change coordinate

params.nElem = 100;         % total number of elements
params.d_factor = 1;        % deflection scale factor for visualization
params.z_stress = [1, 130]; % Z Location of Stress Measurement [m]
params.z_acc = [3, 13, 20]; % Z Location of Acceleration Measurement [m]
params.z_disp = 0;          % Z Location of Displacement Measurement [m]
params.L_total = params.L_pile + params.L_above + params.h; % Total length from -L_pile to +L_above
params.nNode = params.nElem + 1;
params.Le = params.L_total / params.nElem;

%% Guy Wire Properties
params.z_attach = 120;

% Find attachment node
[~, params.j_guy] = min(abs(linspace(-params.L_pile, params.L_above + params.h, params.nNode) - params.z_attach));
params.dof_guy = 2*(params.j_guy-1) + 1; % Horizontal DOF index
params.z_nodeG1 = (params.j_guy); % Z location of the attachment

% Right-side guy wire
params.anchorR_x = 120;                   % Anchor position in +X
params.anchorR_z = 0;                     % Ground level
params.k_guyR = 2e6;                      % [N/m]
params.L0_guyR = sqrt((params.anchorR_x)^2 + (params.anchorR_z - params.z_nodeG1)^2);  % unstretched length

% Left-side guy wire
params.anchorL_x = -120;                  % Anchor in −X direction
params.anchorL_z = 0;
params.k_guyL = 2e6;                      % [N/m]
params.L0_guyL = sqrt((params.anchorL_x)^2 + (params.anchorL_z - params.z_nodeG1)^2);
params.L0_guyR = params.L0_guyR * (1 - 0.1);  % 10% pre-tension
params.L0_guyL = params.L0_guyL * (1 - 0.1);  % 10% pre-tension

%% Creation of Wave Spectrum
params.omega_p = 2*pi/params.T_wave;
params.omega_min = 2*pi ./ params.T_max;   % corresponds to the longest period
params.omega_max = 2*pi ./ params.T_min;   % corresponds to the shortest period

params.omegaVec = linspace(params.omega_min, params.omega_max, params.Nfreq);

% JONSWAP spectrum S(ω)
sigma1 = 0.07; sigma2 = 0.09;
alphaJS = 0.076 * params.H_wave^2 * params.omega_p^4 / params.g^2 * (1 - 0.287*log(params.gammaJ));
S = alphaJS * params.g^2 ./ params.omegaVec.^5 ...
    .* exp(-1.25*(params.omega_p./params.omegaVec).^4) ...
    .* params.gammaJ.^exp(- (params.omegaVec - params.omega_p).^2 ...
                  ./ (2*( (params.omegaVec<params.omega_p)*sigma1^2 ...
                        + (params.omegaVec>=params.omega_p)*sigma2^2 ) ...
                     * params.omega_p^2) );

domega = params.omegaVec(2)-params.omegaVec(1);
m0 = trapz(params.omegaVec, S);              % discrete integral ≈ ∫S dω
S = S * (params.H_wave^2/16) / m0;           % rescale so that ∫S dω = H^2/16

params.A_m = sqrt(2 * S * domega);
params.phiRand = 2*pi*rand(params.Nfreq,1);      % ONE set of random phases
params.k_m_vec = arrayfun(@(w) fzero(@(kk) w^2 - params.g*kk*tanh(kk*params.h), w^2/params.g), params.omegaVec);

% Pre-compute wave number k from dispersion (linear wave theory)
omega = 2*pi / params.T_wave;
func = @(k) omega^2 - params.g * k * tanh(k * params.waterDepth);
params.k = fzero(func, omega^2/params.g);

% Time discretization
params.time = (0:params.nSteps) * params.dt;
params.ramp = min(params.time / params.ramp_duration, 1);  % Linear ramp from 0 to 1 over ramp_duration seconds

end