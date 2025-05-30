% Euler–Bernoulli beam dynamic analysis with two material sections
% using Newmark-beta integration, Rayleigh damping, Morison loads from wave kinematics,
% natural frequency extraction, and snapshot of Morison loads

clear; clc; close all;
% For reproducibility 
rng(30)
set(groot, 'DefaultAxesFontSize', 24); % Set default axes font size
set(groot, 'DefaultTextFontSize', 24); % Set default text font size
%% Flags
video = false ;
guywires = false;
secondOrder = false;
useLinearPY = true ;
irregular = true;
%% Main Parameters
g = 9.81;            % gravity [m/s^2]
ramp_duration   = 20; % Duration of force linear ramp in seconds
dt              = 0.02;     % time step [s]
t_total         = 100;       % total simulation time [s]

% Section & Tip Properties
m_hub           = 1;      % tip mass [kg]
I_hub           = 1e1;      % tip rotational inertia [kg·m^2]
F_hub           = 1e5;      % constant horizontal tip load [N]
F_hub_oscill    =  0;      % oscillatory tip load amplitude [N]
F_hub_freq      = 0.8;      % oscillatory tip load frequency [rad/s]

% Wave Kinematics
gammaJ          = 3.3;      % JONSWAP peak enhancement factor
H_wave          = 1e-6;        % wave height [m]
T_wave          = 16;    % wave period [s]
h               = 23;       % [m], water depth
waterDepth      = h;        % alias for water depth [m]
rho_water       = 1025;     % water density [kg/m^3]
Cd              = 1.0;     % Morison drag coefficient
Cm              = 2.0;     % Morison inertia coefficient
Nfreq           = 128;       % number of spectral components
T_min           = 4;          % Low Period Cutoff
T_max           = 30;          % High Period Cutoff

% Aerodynamic (Air) Properties
rho_air         = 1.225;    % air density [kg/m^3]
Cd_air          = 1.2;      % aerodynamic drag coefficient
U_ref           = 0;       % wind speed at reference height [m/s]
z_ref           = 10;       % reference height for wind profile [m]
alpha           = 0.14;     % wind power-law exponent

% Section Properties
D_outer1        = 7;        % Section 1 outer diameter [m]
thick1          = 0.07;      % Section 1 wall thickness [m]
E1              = 210e9;    % Section 1 Young's modulus [Pa]
rho1            = 7850;     % Section 1 density [kg/m^3]
D_outer2        = 7;        % Section 2 outer diameter [m]
thick2          = 0.07;      % Section 2 wall thickness [m]
E2              = 210e9;    % Section 2 Young's modulus [Pa]
rho2            = 7850;     % Section 2 density [kg/m^3]
D_outer_soil = 7;       % Outer diameter [m]
thick_soil   = 0.07;    % Wall thickness [m]
E_soil       = 210e9;   % Young0s modulus
rho_soil     = 7850;    % Density

% Damping & Integration
alpha_ray       = 0.077;    % Rayleigh damping alpha
beta_ray        = 0.006;    % Rayleigh damping beta
gammaN          = 0.5;      % Newmark-beta gamma
betaN           = 0.25;     % Newmark-beta beta

% Geometry & Mesh
L_pile = 50;                % Embedment depth [m]
L_above = 120;                % Original above-ground length
TowerHeight = 90; 
sectionCutOff   = L_above-TowerHeight;       % [m], section change coordinate

nElem = 300;  % total number of elements
d_factor        = 1;                % deflection scale factor for visualization
z_stress        = [1,130];               % Z Location of Stress Measurement [m]
z_acc           = [3, 13, 20 ] ; 
z_disp          = 0;               % Z Location of Displacement Measurement [m]
L_total = L_pile + L_above + h;   % Total length from -L_pile to +L_above
nNode = nElem + 1;
Le = L_total / nElem;
elemZ = linspace(-L_pile, L_above + h, nNode);  % from -L_pile (soil) to +L_above (tip)

%% GUY WIRE PROPERTIES
z_attach = [120] ;
% Node to attach both wires
j_guy = find(abs(elemZ - z_attach) < 1);  % Attachment node index
[~, j_guy] = min(abs(elemZ - z_attach));
dof_guy = 2*(j_guy-1) + 1;                   % Horizontal DOF index

z_nodeG1 = elemZ(j_guy);  % Z location of the attachment

% --- Right-side guy wire ---
anchorR_x = 120;                   % Anchor position in +X
anchorR_z = 0;                     % Ground level
k_guyR = 2e6;                      % [N/m]
L0_guyR = sqrt((anchorR_x)^2 + (anchorR_z - z_nodeG1)^2);  % unstretched length

% --- Left-side guy wire ---
anchorL_x = -120;                  % Anchor in −X direction
anchorL_z = 0;
k_guyL = 2e6;                   % [N/m]
L0_guyL = sqrt((anchorL_x)^2 + (anchorL_z - z_nodeG1)^2);
L0_guyR = L0_guyR * (1 - 0.1);  % 1% pre-tension
L0_guyL = L0_guyL * (1 - 0.1);  % 2% pre-tension

%% Creation of Wave Spectrum
omega_p  = 2*pi/T_wave;
omega_min = 2*pi ./ T_max;   % corresponds to the longest period
omega_max = 2*pi ./ T_min;   % corresponds to the shortest period

omegaVec = linspace(omega_min, omega_max, Nfreq);

% JONSWAP spectrum S(ω)
sigma1 = 0.07; sigma2 = 0.09;
alphaJS = 0.076 * H_wave^2 * omega_p^4 / g^2 * (1 - 0.287*log(gammaJ));
S = alphaJS * g^2 ./ omegaVec.^5 ...
    .* exp(-1.25*(omega_p./omegaVec).^4) ...
    .* gammaJ.^exp(- (omegaVec - omega_p).^2 ...
                  ./ (2*( (omegaVec<omega_p)*sigma1^2 ...
                        + (omegaVec>=omega_p)*sigma2^2 ) ...
                     * omega_p^2) );

domega   = omegaVec(2)-omegaVec(1);
m0    = trapz(omegaVec, S);              % discrete integral ≈ ∫S dω
S     = S * (H_wave^2/16) / m0;           % rescale so that ∫S dω = H^2/16

A_m   = sqrt(2 * S * domega);
phiRand  = 2*pi*rand(Nfreq,1);      % ONE set of random phases
k_m_vec  = arrayfun(@(w) fzero(@(kk) w^2 - g*kk*tanh(kk*h), w^2/g), omegaVec);

% Pre-compute wave number k from dispersion (linear wave theory)
omega = 2*pi / T_wave;
func = @(k) omega^2 - g * k * tanh(k * waterDepth);
k = fzero(func, omega^2/g);


% Time discretization
nSteps = round(t_total / dt);
time = (0:nSteps) * dt;
ramp = min(time / ramp_duration, 1);  % Linear ramp from 0 to 1 over ramp_duration seconds

%% -------------------- MESH & GLOBAL MATRICES --------------------
nDOF = 2 * nNode;
K = zeros(nDOF);
M = zeros(nDOF);

% assemble element matrices
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

    Ke = E*Iz/Le^3 * [12,6*Le,-12,6*Le;6*Le,4*Le^2,-6*Le,2*Le^2;-12,-6*Le,12,-6*Le;6*Le,2*Le^2,-6*Le,4*Le^2];
    Me = rho*A_cs*Le/420 * [156,22*Le,54,-13*Le;22*Le,4*Le^2,13*Le,-3*Le^2;54,13*Le,156,-22*Le;-13*Le,-3*Le^2,-22*Le,4*Le^2];

    dofs = [2*(e-1)+1, 2*(e-1)+2, 2*e+1, 2*e+2];
    K(dofs, dofs) = K(dofs, dofs) + Ke;
    M(dofs, dofs) = M(dofs, dofs) + Me;
end


% lump tip mass into M & damping
topDOF = 2*(nNode-1) + 1;
M(topDOF, topDOF) = M(topDOF, topDOF) + m_hub;
topRotDOF = topDOF + 1;    % rotational DOF of the top node
M(topRotDOF, topRotDOF) = M(topRotDOF, topRotDOF) + I_hub;

C = alpha_ray*M + beta_ray*K;

%% AERODYNAMIC LOAD
% index of your tip‐bending DOF
tipDispDOF = 2*(nNode-1) + 1;

% define your constant horizontal tip load [N]

% make a constant‐load vector
F_const = zeros(nDOF,1);
% Cos Unsteady (TestO=)
%F_const(tipDispDOF) = F_hub.*cos(0.4.*time);
% Steady
F_const(tipDispDOF) = F_hub ;

%% -------------------- TIME INTEGRATION (Newmark-Beta) --------------------
% boundary and initial DOFs
fixed = [1,2];
all    = 1:nDOF;
free   = setdiff(all, fixed);

U = zeros(nDOF, nSteps+1);
V = zeros(nDOF, nSteps+1);
Aacc = zeros(nDOF, nSteps+1);

% initial acceleration at t = 0
F0 = computeMorisonForces(elemZ, Nfreq, omegaVec, A_m, phiRand, k_m_vec, ...
                             sectionCutOff, D_outer1, D_outer2, ...
                             Cd, Cm, rho_water, time(1),h,H_wave,T_wave,Le,secondOrder,irregular ) ;
Fm_air  = computeAerodynamicForces(elemZ, sectionCutOff, D_outer1, D_outer2, h, Le,U_ref,z_ref,alpha,Cd_air,rho_air);
F_py = zeros(nDOF, 1);

ramp_i = ramp(1);
F0 = (F0 + F_const + Fm_air) * ramp_i;
Aacc(free, 1) = M(free, free) \ (F0(free) - C(free, free) * V(free, 1) - K(free, free) * U(free, 1));

% Newmark constants
a0c = 1/(betaN * dt^2);
a1c = gammaN / (betaN * dt);
a2c = 1/(betaN * dt);
a3c = 1/(2*betaN) - 1;
a4c = gammaN / betaN - 1;
a5c = dt * (gammaN/(2*betaN) - 1);
Keff = K(free, free) + a0c*M(free, free) + a1c*C(free, free);
FHyd_history = zeros(nDOF, nSteps+1);
FAero_history = zeros(nDOF, nSteps+1);
% Initialize soil state for each node
emptyState = struct( ...
    'y_prev', 0, ...
    'p_prev', 0, ...
    'direction', 0, ...
    'y_rev', 0, ...
    'p_rev', 0, ...
    'k_unload', 0 ...
);

soilStates = repmat(emptyState, nNode, 1);
thetaStates = repmat(struct( ...
    'theta_prev', 0, ...
    'm_prev', 0, ...
    'direction', 0), nNode, 1);
j_watch = 28;  % node index you want to track
dof_watch = 2*(j_watch - 1) + 1;

p_record = zeros(1, nSteps+1);
y_record = zeros(1, nSteps+1);
Fsoil_history = zeros(nDOF, nSteps+1);  % Store F_py for all time steps
tic
[Fi,eta2] = computeMorisonForces( elemZ, Nfreq, omegaVec, A_m, phiRand, k_m_vec, ...
                             sectionCutOff, D_outer1, D_outer2, ...
                             Cd, Cm, rho_water, time,h,H_wave,T_wave,Le,secondOrder,irregular );
toc


tic;
%%
for i = 1:nSteps
    % print progress
    fprintf('Step %d / %d   t = %.2f s\n', i, nSteps, time(i));

    ramp_i = ramp(i);
    F_hydro = Fi(:, i);  % ← hydro force at current time step
    Fm_air  = computeAerodynamicForces(elemZ, sectionCutOff, D_outer1, D_outer2, h, Le,U_ref,z_ref,alpha,Cd_air,rho_air);
    F_const(tipDispDOF) = F_hub + F_hub_oscill * cos(F_hub_freq * time(i));
    F_py = zeros(nDOF, 1);
    for j = 1:nNode
        z = elemZ(j);
        c_ref = 5e4;  % [N·s/m³]
    
        if z < 0
            dof = 2*(j-1) + 1;
            dof_rot = 2*(j-1)+2 ;
            y_disp = U(dof, i);
            y_vel  = V(dof, i);  % ← nodal velocity
            depth  = -z;
        
            % Spring force (nonlinear)
            [p_spring, soilStates(j)] = pySimpleHysteretic(y_disp, depth, soilStates(j),useLinearPY);
        
            % Damping coefficient (linear with depth)
            c_ref = 0;               % N·s/m³, standard for sand
            c_damp = c_ref * depth;    % [N·s/m²]
        
            % Damping force
            f_damp = -c_damp * y_vel;
        
            % Total soil force (spring + damping)
            p_total = p_spring + f_damp;
        
            % Optional: record values for diagnostics
            if j == j_watch
                y_record(i) = y_disp;
                p_record(i) = p_total;
            end
        
            % Apply total force to global vector
            F_py(dof) = F_py(dof) - p_total * Le;
            theta = U(dof_rot, i);        % current rotation
            theta_dot = V(dof_rot, i);    % current rotational velocity

            [m_spring, thetaStates(j)] = mrSimpleHysteretic(theta, depth, thetaStates(j),useLinearPY);
            % 2) rotational damping coefficient (per unit depth)
            %    here we start from the same c_ref [N·s/m³] and convert to N·m·s/rad
            D_pile = D_outer;                    % pile diameter at this depth
            c_ref  = 0;                        % N·s/m³ (translational base)
            c_rot_ref = c_ref * (D_pile/2)^2;    % ≃ N·m·s/rad per m depth
            c_rot    = c_rot_ref * depth;        % [N·m·s/rad]
        
            % 3) damping moment opposing motion
            m_damp = -c_rot * theta_dot;
        
            % 4) total moment
            m_total = m_spring + m_damp;
        
            % 5) apply to global vector (moment × element length)
            F_py(dof_rot) = F_py(dof_rot) - m_total * Le;
        end
        % Current lateral displacement at the node
        if guywires == true
            x_node = U(dof_guy, i);
            x_dot  = V(dof_guy, i);  % ← velocity at the node
            zeta_guy = 0.03;  % 3% damping
            m_eff = 1e5;      % [kg] effective lumped mass at node (estimate)
            
            % Right cable
            c_guyR = 1e6;
            
            % Left cable
            c_guyL = 1e6;
            dxR = anchorR_x - x_node;
            dzR = anchorR_z - z_nodeG1;
            L_R = sqrt(dxR^2 + dzR^2);
            
            if L_R > L0_guyR
                eps_smooth = 0.01;
                delta_L_R = max(0, L_R - L0_guyR);
                T_R = k_guyR * tanh(delta_L_R / eps_smooth) * eps_smooth;
            
                % Direction
                dirR = dxR / L_R;
            
                % Damping force (opposes velocity)
                F_damp_R = -c_guyR * x_dot * dirR;
            
                % Total horizontal force
                Fx_R = T_R * dirR + F_damp_R;
            
                F_py(dof_guy) = F_py(dof_guy) - Fx_R;
            end
            dxL = anchorL_x - x_node;
            dzL = anchorL_z - z_nodeG1;
            L_L = sqrt(dxL^2 + dzL^2);
            
            if L_L > L0_guyL
                eps_smooth = 0.01;
                delta_L_L = max(0, L_L - L0_guyL);
                T_L = k_guyL * tanh(delta_L_L / eps_smooth) * eps_smooth;
            
                dirL = dxL / L_L;
                F_damp_L = -c_guyL * x_dot * dirL;
            
                Fx_L = T_L * dirL + F_damp_L;
            
                F_py(dof_guy) = F_py(dof_guy) - Fx_L;
            end
    
        end
                Fsoil_history(:, i) = F_py;

    end
    
    Fi_t    = ramp_i * (F_hydro + F_const + Fm_air + F_py);  % full hydro + aero + soil
    Ff = Fi_t(free);
    Feff = Ff + M(free, free) * (a0c*U(free, i) + a2c*V(free, i) + a3c*Aacc(free, i)) + ...
           C(free, free) * (a1c*U(free, i) + a4c*V(free, i) + a5c*Aacc(free, i));
    U(free, i+1) = Keff \ Feff;
    Aacc(free, i+1) = a0c*(U(free, i+1) - U(free, i)) - a2c*V(free, i) - a3c*Aacc(free, i);
    V(free, i+1) = V(free, i) + dt*((1-gammaN)*Aacc(free, i) + gammaN*Aacc(free, i+1));
    FHyd_history(:,i) = F_hydro;
    FAero_history(:,i) = Fm_air + F_const;
    

end
toc;

%% -------------------- POST-PROCESSING: TIP RESPONSE & STRESS --------------------
figure;
subplot(2,1,1)
plot(time, U(topDOF, :), 'LineWidth', 4);
xlabel('Time [s]'); ylabel('Tip disp [m]'); grid on;
grid minor;
[max_disp, max_idx] = max(U(topDOF, :));
fprintf('Maximum tip displacement: %.4f m at t = %.2f s\n', max_disp, time(max_idx));

nLocs    = numel(z_stress);

% Find the nearest element index for each z_stress
stress_plot = arrayfun(@(z) ...
    find(abs(elemZ - z) == min(abs(elemZ - z)), 1), ...
    z_stress);

% Preallocate: each row is one location, each column one time step
stress_time = zeros(nLocs, nSteps+1);

% Loop over each location
for iLoc = 1:nLocs
    node_loc = stress_plot(iLoc);
    
    % pick material properties based on your cutoff
    if z_stress(iLoc) <= sectionCutOff
        D_loc = D_outer1; E_loc = E1;
    else
        D_loc = D_outer2; E_loc = E2;
    end
    c = D_loc/2;
    
    % compute stress time history at this node
    for j = 1:(nSteps+1)
        % degrees of freedom: i = this node, p = previous, n = next
        idx_i = 2*(node_loc-1) + 1;
        idx_p = 2*(node_loc-2) + 1;
        idx_n = 2*(node_loc)   + 1;
        
        curvature = ( U(idx_p, j) ...
                    - 2*U(idx_i, j) ...
                    +    U(idx_n, j) ) ...
                    / Le^2;
        stress_time(iLoc, j) = E_loc * curvature * c;
    end
end

% Plot all locations on one figure
subplot(2,1,2)
 hold on; grid on; grid minor;
colors = lines(nLocs);
for iLoc = 1:nLocs
    plot(time, stress_time(iLoc,:), ...
         'LineWidth',5, ...
         'Color', colors(iLoc,:));
end
xlabel('Time [s]');
ylabel('Bending stress [Pa]');
legend(arrayfun(@(z) sprintf('z = %.1f m', z), z_stress, ...
               'UniformOutput',false), ...
       'Location','Best');
% find the node nearest z=50 m
[~, node_plot] = min(abs(elemZ - z_disp));

% vertical‐disp. DOF for that node
dispDOF = 2*(node_plot-1) + 1;

% plot its time history
figure;
plot(time, U(dispDOF, :), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Displacement [m]');
title(sprintf('Beam displacement at z = %.1f m', elemZ(node_plot)));
grid on;
grid minor;



%


%% --- Nodal acceleration time series at specified z_stress locations ---
acc_plot = arrayfun(@(z) ...
    find(abs(elemZ - z) == min(abs(elemZ - z)), 1), ...
    z_acc);
nLocs    = numel(z_acc);

accel_time = zeros(nLocs, nSteps+1);

for iLoc = 1:nLocs
    node_loc = acc_plot(iLoc);   % same nearest node index as above

    % direct grab of acceleration time series for this node
    accel_time(iLoc, :) = Aacc(node_loc, :);
end

% Plot acceleration at each location
figure; hold on; grid on;
colors = lines(nLocs);
for iLoc = 1:nLocs
    plot(time, accel_time(iLoc,:), ...
         'LineWidth',1.5, ...
         'Color', colors(iLoc,:));
end
xlabel('Time [s]');
ylabel('Nodal acceleration [m/s^2]');
legend(arrayfun(@(z) sprintf('z = %.1f m', z), z_stress, ...
               'UniformOutput',false), ...
       'Location','Best');
title('Nodal Acceleration at Selected Locations');
%% -------------------- NATURAL FREQUENCIES --------------------
% Use reduced matrices to avoid singularity
K_red = K(free, free);
M_red = M(free, free);

eigModes = 6;  % number of modes to extract
% Shift-invert around zero to get smallest nonzero frequencies
[Phi_red, Omega2_red] = eigs(K_red, M_red, eigModes, 'smallestabs');
omega_n = sqrt(diag(Omega2_red));  % rad/s
freq_hz = omega_n / (2*pi);       % Hz
fprintf('First natural frequency: %.3f Hz\n', freq_hz(1));
fprintf('Second natural frequency: %.3f Hz\n', freq_hz(2));
omega1 = 2*pi*freq_hz(1);    % ≃1.602 rad/s
omega2 = 2*pi*freq_hz(2);    % ≃10.542 rad/s
zeta   = 0.03;

alpha_ray = 2*zeta*omega1*omega2/(omega1+omega2);   % ≃0.0835
beta_ray  = 2*zeta/(omega1+omega2);                 % ≃0.00494

fprintf('Recomended Alpha:%.3f \n', alpha_ray);
fprintf('Recomended Beta:%.3f \n', beta_ray);

% %% Surface Elevation]
nSteps = round(t_total/dt);
time   = (0:nSteps) * dt;
eta_irreg = zeros(1, nSteps+1);
omegaVec = omegaVec(:);   % Nfreq×1
A_m      = A_m(:);        % Nfreq×1
phiRand  = phiRand(:);    % Nfreq×1

for i = 1:nSteps+1
    t = time(i);
    eta_irreg(i) = sum( A_m .* cos(omegaVec*t + phiRand) );
end
%%
%%— Plot both on the same figure —%%
figure; hold on; grid on;
plot(time, eta_irreg,time,eta_irreg+eta2, 'LineWidth', 2.0, 'DisplayName', 'Irregular \eta');
xlabel('Time [s]');
ylabel('Surface elevation \eta [m]');
title('Wave surface elevation time–series');

% 
% 
% %%
% 
% F_vert_hist = FHyd_history(1:2:end, :);      % nNode × (nSteps+1)
% F_Hydro     = sum(F_vert_hist,1) ;    % 1×(nSteps+1)
% figure;
% plot(time, F_Hydro, 'LineWidth', 1.5);
% xlabel('Time [s]');
% ylabel('Total horizontal hydro-force [N]');
% grid on;
% 
% F_vert_hist_aero = FAero_history(1:2:end, :);      % nNode × (nSteps+1)
% F_Aero     = sum(F_vert_hist_aero,1) ;    % 1×(nSteps+1)
% figure;
% plot(time, F_Aero, 'LineWidth', 1.5);
% xlabel('Time [s]');
% ylabel('Total horizontal aero-force [N]');
% grid on;
% 
% 
% %% Second Order Regular Wave Plot
% % % Precompute k & Stokes coefficients
% % omega0 = 2*pi/T_wave;
% % % solve dispersion for k0 once
% % k0     = fzero(@(kk) omega0^2 - g*kk*tanh(kk*waterDepth), omega0^2/g);
% % a      = H_wave/2;                    % first-order amp
% % % 2nd-order coefficient (Stokes drift term)
% % B2     = (a^2 * k0) / (2 * sinh(k0*waterDepth)^2);
% % 
% % % build surface elevation time series
% % eta_stokes = zeros(1, nSteps+1);
% % for i = 1:nSteps+1
% %   t = time(i);
% %   eta1 = a * cos(omega0*t);
% %   eta2 = B2 * (cosh(2*k0*waterDepth) / sinh(k0*waterDepth)^2) * cos(2*omega0*t);
% %   eta_stokes(i) = eta1 + eta2;
% % end
% % figure;
% % plot(time, eta_stokes, 'LineWidth', 1.5);
% % xlabel('Time [s]');
% % ylabel('\eta(t) [m]');
% % title('2nd-Order Stokes Surface Elevation');
% % grid on;
% 
% %% Save relevant data
% 
% 
% 
% 
% 
% 
% 
secondOrder = [time;eta2;stress_time;U(topDOF, :)   ];
% 
% 
% 
% 
% 
% 
% 
% 
% 
%% Generate Figure
data = U(topDOF, :);        % full tip–displacement vector
N    = length(data);        % total number of samples
startIdx = floor(0.70*N) + 1;  % index where the last 70% begins
[umax_rel, relIdx] = max( data(startIdx:end) );

% convert back to the full‐series index
idx_max = startIdx - relIdx - 1;
umax = data(idx_max);
t_max = time(idx_max);

fprintf('Peak tip displacement %.3f m at t = %.3f s\n', umax, t_max);


% 2) Get the nodal deflections at that instant
U_nodes = U(1:2:end, idx_max);    % vertical DOFs only

% 3) Compute bending stress at each node via second‐difference curvature
stress_nodes = zeros(nNode,1);
for j = 2:nNode-1
    %  curvature at node j
    kappa = (U_nodes(j-1) - 2*U_nodes(j) + U_nodes(j+1)) / Le^2;
    %  section properties at this node
    if elemZ(j) <= sectionCutOff
        E_j    = E1;
        D_j    = D_outer1;
    else
        E_j    = E2;
        D_j    = D_outer2;
    end
    c_j        = D_j/2;
    stress_nodes(j) = E_j * c_j * kappa;
end
%  mirror end‐values
stress_nodes(1)     = stress_nodes(2);
stress_nodes(end)   = stress_nodes(end-1);

% 4) Plot deflected shape colored by stress
figure('Color','w'); 
scatter(U_nodes.*d_factor, elemZ, 50, stress_nodes, 'filled'); hold on;
plot(U_nodes.*d_factor, elemZ, 'k-', 'LineWidth',1.5);
cb = colorbar;
cb.Ticks = (min(stress_nodes)*0.8:1e7: max(stress_nodes)*1.2);
cb.Label.String = 'Bending Stress [Pa]';
colormap(jet);
xlabel('Deflection [m]');
ylabel('Elevation [m]');
title(sprintf('SF = %.1f', d_factor));
grid on; axis equal tight;
%% Beam Theory Comparison
% theoretical deflection at each z
w_th = F_hub .* elemZ.^2 .* (3*L_total - elemZ) ./ (6 * E1 * Iz);

% now plot side by side
figure('Color','w'); hold on; grid on;
% your numerics at t_max
plot(U_nodes, elemZ, 'r-','LineWidth',3,'DisplayName','MonoSim');
% beam‐theory curve
plot(w_th, elemZ, 'k--','LineWidth',3,'DisplayName','Beam Theory Eq.');
xlabel('Horizontal deflection [m]');
ylabel('Elevation [m]');
legend('Location','best');

% Stress Comparison
M_th      = F_hub * (L_total - elemZ);            
% section modulus distance (outer fiber)
c         = D_outer1/2;                            
% theoretical bending stress
sigma_th  = M_th .* c ./ Iz;                       

% --- Plot FE‐model stress vs. theoretical stress ---
figure('Color','w'); hold on; grid on;
% FE‐model stress (change 'sigma_model' to your actual variable)
plot(stress_nodes, elemZ, 'r-','LineWidth',3,'DisplayName','MonoSim');
% theoretical beam stress
plot(sigma_th,     elemZ, 'k--','LineWidth',3,'DisplayName','Beam Theory Eq.');
xlabel('Bending Stress [Pa]');
ylabel('Elevation z [m]');
legend('Location','best');


%%

%Parameters
data = U(topDOF, :);        % full tip–displacement vector
N    = length(data);        % total number of samples
startIdx = floor(0.70*N) + 1;  % index where the last 70% begins
[umax_rel, relIdx] = max( data(startIdx:end) );

% convert back to the full‐series index
idx_max = startIdx - relIdx - 1;
umax = data(idx_max);
t_max = time(idx_max);

fprintf('Peak tip displacement %.3f m at t = %.3f s\n', umax, t_max);


% 2) Get the nodal deflections at that instant
U_nodes = U(1:2:end, idx_max);    % vertical DOFs only

% 3) Compute bending stress at each node via second‐difference curvature
stress_nodes = zeros(nNode,1);
for j = 2:nNode-1
    %  curvature at node j
    kappa = (U_nodes(j-1) - 2*U_nodes(j) + U_nodes(j+1)) / Le^2;
    %  section properties at this node
    if elemZ(j) <= sectionCutOff
        E_j    = E1;
        D_j    = D_outer1;
    else
        E_j    = E2;
        D_j    = D_outer2;
    end
    c_j        = D_j/2;
    stress_nodes(j) = E_j * c_j * kappa;
end
%  mirror end‐values
stress_nodes(1)     = stress_nodes(2);
stress_nodes(end)   = stress_nodes(end-1);
nCirc   = 36;           % number of points around circumference
theta   = linspace(0,2*pi,nCirc);
[Z,Th]  = meshgrid(elemZ, theta);  % theta×z

% pick radius for each node
R = zeros(size(Z));
for j=1:numel(elemZ)
  if elemZ(j) <= sectionCutOff
    R(:,j) = D_outer1/2;
  else
    R(:,j) = D_outer2/2;
  end
end

% replicate deflection & stress around circumference
U_mat      = d_factor.*repmat(U_nodes(:).', size(theta,2), 1);     % theta×z
stress_mat = repmat(stress_nodes(:).', size(theta,2),1); % theta×z

% coordinates of the warped cylinder
X = R .* cos(Th) + U_mat;  % warp the X‐coordinate by deflection
Y = R .* sin(Th);
% Z stays the same

% Plot
figWidth  = 1080;   % Desired width in pixels
figHeight = 1080;   % Desired height in pixels
figure('Color','w', 'Position', [100, 100, figWidth, figHeight]);

hSurf = surf(X, Y, Z, stress_mat, ...
    'FaceColor','interp', 'EdgeColor','none');
hold on
rHub = 5;   % Example radius, set to your value
[xs, ys, zs] = sphere(40);  % creates a unit sphere with 40x40 resolution

% Position of hub: X = deflected, Y = 0, Z = elemZ(end)
Xhub = xs * rHub + d_factor * U_nodes(end);  % add deflection in X
Yhub = ys * rHub;
Zhub = zs * rHub + elemZ(end);

surf(Xhub, Yhub, Zhub, 'FaceColor',[0.4 0.7 1], 'EdgeColor','none', 'FaceAlpha', 0.8);

% (optional) improve lighting and plot appearance
camlight; lighting gouraud;
t0      = time(idx_max);
% now the water plane at t_max
x_range = linspace(-50, 50, 200);   % 200 points from -50 to 50 m
y_range = linspace(-50, 50, 200);   % ditto in y
[Xw, Yw] = meshgrid(x_range, y_range);
eta_xy  = zeros(size(Xw));

for m = 1:Nfreq
  km    = k_m_vec(m);              % wave number for mode m
  Am    = A_m(m);                  % amplitude
  phm   = phiRand(m);              % random phase
  om    = omegaVec(m);             % rad/s
  % phase field kx - omega t + phi
  phase = km*Xw - om*t0 + phm;     
  eta_xy = eta_xy + Am * cos(phase);
end
surf(Xw, Yw, eta_xy+h, 'FaceColor','cyan', 'EdgeColor','none', ...
     'FaceAlpha', 0.5);

% refine view
colormap(jet); colorbar;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
%title(sprintf('MonoSim @ t=%.3f s', t_max));
view(-90,0);
axis equal tight;
camlight headlight; lighting gouraud;
band_radius = D_outer1/2 * 1.05;   % slightly larger so bands are visible
band_width  = 2;                 % width of band in z-direction

z_band_locs = [3 13 20];           % z locations for bands

nCirc_band = 60;                   % smoothness of the band
theta_band = linspace(0, 2*pi, nCirc_band);

for k = 1:numel(z_band_locs)
    % Center z for this band (find closest in elemZ)
    [~, iz] = min(abs(elemZ - z_band_locs(k)));
    z_center = elemZ(iz);

    % Range for band: band_width/2 above and below center
    z_band = linspace(z_center-band_width/2, z_center+band_width/2, 6);

    [thetaGrid, zGrid] = meshgrid(theta_band, z_band);
    Xb = band_radius * cos(thetaGrid);
    Yb = band_radius * sin(thetaGrid);
    Zb = zGrid;

    surf(Xb, Yb, Zb, ...
        'FaceColor','k', 'EdgeColor','none', 'FaceAlpha', 0.7);
end

%% 2D Animation
% v = VideoWriter('morisonLoads_deflected.mp4','MPEG-4');
% v.FrameRate = 10;
% open(v);
% 
% figure('Color','w');
% for i = 1:10:nSteps+1
%     clf; hold on; grid on;
% 
%     % 1) plot deflected tower (scaled)
%     U_def = d_factor * U(1:2:end, i);  % U at each node, scaled
%     plot(U_def, elemZ, 'k-', 'LineWidth', 2);
% 
%     % 2) morison loads as arrows at the deflected tower
%     scale = Le*2 / max(abs(FHyd_history(:)));
%     for j = 1:nNode
%         Fji = FHyd_history(2*(j-1)+1, i)+FAero_history(2*(j-1)+1, i);
%         quiver(U_def(j), elemZ(j), Fji*scale, 0,0, 'r', ...
%                'MaxHeadSize', 3, 'LineWidth', 1);
%     end
% 
%     title(sprintf('Loads & Deflection (×%d) at t = %.2f s', d_factor, time(i)));
%     xlabel('Horizontal [m]'); ylabel('Elevation [m]');
%     axis tight; ylim([0 L_total+20]);xlim([-2 2])
%     
%     frame = getframe(gcf);
%     writeVideo(v, frame);
% end
% close(v);

% 3D Animation

%Pre‐compute geometry that doesn%t change
nCirc    = 36;
theta    = linspace(0,2*pi,nCirc)';
[Zb, Th] = meshgrid(elemZ, theta);
Rb       = zeros(size(Zb));
for j = 1:nNode
  Rb(:,j) = (elemZ(j)<=sectionCutOff) * (D_outer1/2) + (elemZ(j)>sectionCutOff) * (D_outer2/2);
end
F_py = Fsoil_history;
Fq = F_py(1:2:end);  % horizontal DOFs only

if video == true
    % Set up video
    v = VideoWriter('beam_wave_animation_guywire.mp4','MPEG-4');
    v.FrameRate = 30;
    open(v);
    smin = -5e8;
    smax = 5e8;
    scaleFactor = 1e-4;       % Adjust to get reasonable arrow length
    
    figWidth  = 720;   % Desired width in pixels
    figHeight = 1080;   % Desired height in pixels
    figure('Color','w', 'Position', [100, 100, figWidth, figHeight], ...
        'MenuBar', 'none', 'ToolBar', 'none');
    for k = floor(nSteps/2):10:nSteps+1
    
    %for k = 1:20:nSteps+1
      t0 = time(k);
      
      % 1) beam deflection & stress at time k
      Uk         = U(1:2:end, k);             % nodal vertical disp
      % recompute stress at each node
      stress_k   = zeros(nNode,1);
      for j = 2:nNode-1
        kappa = (Uk(j-1) - 2*Uk(j) + Uk(j+1)) / Le^2;
        if elemZ(j)<=sectionCutOff
          Ej = E1; Dj = D_outer1;
        else
          Ej = E2; Dj = D_outer2;
        end
        cj       = Dj/2;
        stress_k(j) = Ej * cj * kappa;
      end
      stress_k([1,end]) = stress_k(2)*0.5 + stress_k(end-1)*0.5;  % endpoints
    
      % expand around circumference
      Uk_mat     = d_factor.*repmat(Uk', nCirc, 1);
      stress_mat = repmat(stress_k', nCirc, 1);
      Xb = Rb .* cos(Th) + Uk_mat;
      Yb = Rb .* sin(Th);
    
      % 2) wave surface at time k
      xw = linspace(-50,50,200);
      yw = linspace(-50,50,200);
      [Xw,Yw] = meshgrid(xw,yw);
      eta_xy = zeros(size(Xw));
      for m = 1:Nfreq
        phase = k_m_vec(m)*Xw - omegaVec(m)*t0 + phiRand(m);
        eta_xy = eta_xy + A_m(m)*cos(phase);
      end
    
          % === 3) PLOT ===
        clf; hold on; grid on;
        
        % --- Plot monopile deformation + stress ---
        surf(Xb, Yb, Zb, stress_mat, 'FaceColor','interp', 'EdgeColor','none');
        
        % --- Transparent mudline plane at z = 0 ---
        [Xm, Ym] = meshgrid(linspace(-50, 50, 2), linspace(-50, 50, 2));
        Zm = zeros(size(Xm));
        surf(Xm, Ym, Zm, 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % --- Plot wave surface ---
        surf(Xw, Yw, eta_xy + h, 'FaceColor','cyan', 'EdgeColor','none', 'FaceAlpha',0.5);
        
        % --- Soil restoring force arrows below mudline ---
        Zq = elemZ(:);
        Xq = zeros(size(Zq));
        Yq = zeros(size(Zq));
        Fq = F_py(1:2:end, k);
        
        % Only nodes below mudline
        mask = Zq < 0;
        Zq = Zq(mask);
        Xq = Xq(mask);
        Yq = Yq(mask);
        Fq = Fq(mask);
        
        Uq = Fq * scaleFactor;  % x-direction arrows
        quiver3(Xq, Yq, Zq, Uq, zeros(size(Uq)), zeros(size(Uq)), 0, ...
                'k', 'LineWidth', 1.5, 'MaxHeadSize', 4.5);
        scatter3(Xq + Uq, Yq, Zq, 50, abs(Fq).*100, 'filled');
        
        % --- Plot guy wires as lines ---
        % Get current pile head location
        x_tip  = U(dof_guy, k);   % horizontal displacement of guy node
        z_tip  = elemZ(j_guy);    % vertical coordinate of guy node (fixed)
        if guywires == true

            % Right guy wire
            plot3([x_tip, anchorR_x], [0, 0], [z_tip, anchorR_z], ...
              'r-', 'LineWidth', 2, 'DisplayName','Right Guy Wire');
        
            % Left guy wire
            plot3([x_tip, anchorL_x], [0, 0], [z_tip, anchorL_z], ...
              'b-', 'LineWidth', 2, 'DisplayName','Left Guy Wire');
        end
        % --- Axis, colormap, lighting ---
        colormap(jet);
        cb = colorbar;
        caxis([smin, smax]);
        cb.Ticks = linspace(smin, smax, 5);
        cb.Label.String = 'Bending stress [Pa]';
        
        xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
        title(sprintf('Time = %.2f s', t0));
        view(-45, 15);
        axis equal tight;
        camlight headlight; lighting gouraud;
        
        % Optional: Add legend for guy wires
        %legend('Location', 'northeast');
        
        % Write frame to video
        frame = getframe(gcf);
        writeVideo(v, frame);
    
    end

    close(v);
    disp('Animation saved to beam_wave_animation.mp4');
else
end

% figure; hold on; grid on;
% plot(y_record, p_record, 'b-', 'LineWidth', 2);
% xlabel('Lateral displacement y [m]', 'FontSize', 14);
% ylabel('Soil reaction p [N/m]', 'FontSize', 14);
% title(sprintf('Hysteresis Loop at Node %d (z = %.2f m)', j_watch, elemZ(j_watch)), 'FontSize', 16);
% 
% % Assume: y_record, p_record, time are all (1×nSteps+1)
% t_vec = time(1:nSteps);         % Time vector
% y_vec = y_record(1:nSteps);     % Displacement
% p_vec = p_record(1:nSteps);     % Soil reaction
% 
% % Plot hysteresis with color by time
% figure; hold on; grid on;
% 
% scatter3(y_vec, p_vec, t_vec, 20, t_vec, 'filled');  % 3D scatter, color = time
% colormap(jet); colorbar;
% xlabel('Displacement y [m]', 'FontSize', 14);
% ylabel('Soil Resistance p [N/m]', 'FontSize', 14);
% zlabel('Time [s]', 'FontSize', 14);
% title(sprintf('Hysteresis at Node %d (z = %.2f m)', j_watch, elemZ(j_watch)), 'FontSize', 16);
% view(45, 25);
% 
% % --- Overlay Backbone Curve ---
% depth = -elemZ(j_watch);              % positive depth
% phi_eff = 35;                         % deg
% gamma   = 10e3;                       % N/m^3
% k_init  = 1e6;                        % N/m^3
% pult    = gamma * depth * 3 * tan(deg2rad(phi_eff));
% 
% y_back = linspace(min(x_vec), max(y_vec), 500);
% p_back = pult * tanh(k_init * depth * y_back / pult);
% 
% plot3(y_back, p_back, max(t_vec)*ones(size(y_back)), 'k--', 'LineWidth', 2, ...
%     'DisplayName', 'Backbone Curve');
% 
% legend('Hysteresis Colored by Time', 'Backbone Curve');

%% Pile Deflection Curve
% --- Parameters ---
E  = E_soil;           % Young's modulus [Pa]
I  = pi/64 * (D_outer_soil^4 - (D_outer_soil-2*thick_soil)^4);
EI = E * I;          % Flexural rigidity [Nm^2]

k  = 1e6;            % Modulus of subgrade reaction [N/m^2]
P  = F_hub  ;         % Lateral load at pile head [N]
Lm = L_above ;
L  = L_pile;             % Pile length [m]
M0 = P*Lm ;
lambda = (k / (4 * EI))^(1/4);

% Define symbolic variable and unknown coefficients
syms z A B C D real

% General solution (symbolic)
y_sym = A*cosh(lambda*z) + B*sinh(lambda*z) + C*cos(lambda*z) + D*sin(lambda*z);

% Derivatives
dy  = diff(y_sym, z);
d2y = diff(y_sym, z, 2);
d3y = diff(y_sym, z, 3);

% Apply boundary conditions
eqs = [ ...
    -EI * subs(d2y, z, 0) == M0;     % Moment at top
    -EI * subs(d3y, z, 0) == P;      % Shear at top
    subs(y_sym, z, L) == 0;         % Displacement at bottom
    subs(dy, z, L) == 0             % Rotation at bottom
];

% Solve the system
S = solve(eqs, [A, B, C, D]);

% Substitute solution into y(z)
y_full = simplify(subs(y_sym, [A, B, C, D], [S.A, S.B, S.C, S.D]));

% Convert to numeric function
y_func = matlabFunction(y_full, 'Vars', z);

% Plotting
z_vals = linspace(0, L, 300);
y_vals = y_func(z_vals);
%%
set(groot,'defaultAxesFontSize',32)
set(groot,'defaultTextFontSize',32)
figure('Color','w');
plot(-y_vals, -z_vals, 'b-', 'LineWidth', 4);
hold on
plot(U(1:2:end,end),elemZ, 'LineWidth', 4)
xlabel('y(z) [m]', 'Interpreter', 'latex');
ylabel('Depth [m]', 'Interpreter', 'latex');
grid on; box on;
axis([-0.001 0.035 -L 0])
xline(0,'k--','lineWidth',2)
legend 'Analytical' 'MonoSim'
% subplot(2,1,2)
% plot(M, -z, 'r-', 'LineWidth', 2);
% xlabel('Moment M(z) [Nm]', 'Interpreter', 'latex');
% ylabel('Depth [m]', 'Interpreter', 'latex');
% title('Bending Moment Along Pile', 'Interpreter', 'latex');
% grid on; box on;
% set(gca, 'YDir','reverse', 'FontSize', 14);


%% ------------------- HELPER FUNCTION --------------------
% 2) Modify computeMorisonForces to take those in:
function [Fm_wave, eta2] = computeMorisonForces(z, Nfreq, omegaVec, A_m, phiRand, k_m_vec, ...
                                                cutoff, D1, D2, Cd, Cm, rhoW, t, h, ...
                                                Hs, Tp, Le, secondOrder, irregular)
    g      = 9.81;
    nNode  = numel(z);
    Fm_wave = zeros(2*nNode, numel(t));
    % right after you unpack Nfreq, etc.
    thetaVec = zeros(1, Nfreq);   % all waves traveling in x–direction

    for i = 1:nNode
        z_b = z(i);
        if z_b < 0 || z_b > h, continue; end
        z_w = z_b - h;

        % ——— First‐order U & Udot (regular vs irregular) —————
        if ~irregular
            omega0 = 2*pi/Tp;
            k0     = omega0^2/g/tanh(omega0^2/g/tanh(omega0*h)*h);
            a      = Hs/2;
            C1     = cosh(k0*(z_w+h))/sinh(k0*h);
            U1     = a*omega0*C1.*cos(omega0*t);
            Udot1  = -a*omega0^2*C1.*sin(omega0*t);
            U2 = zeros(size(t)); Udot2 = zeros(size(t));

            if secondOrder
                % if you need 2nd-order velocity in Morison forces,
                % you could call a secondOrderVelocity() here instead.
                % For now, keep U2=0 so only η² is returned.
            end

            U    = U1 + U2;
            Udot = Udot1 + Udot2;
        else
            U    = zeros(size(t));
            Udot = zeros(size(t));
            for m = 1:Nfreq
                decay = cosh(k_m_vec(m)*(z_w + h)) / sinh(k_m_vec(m)*h);
                U    = U + A_m(m)*decay.*cos(omegaVec(m)*t + phiRand(m));
                Udot = Udot - A_m(m)*omegaVec(m)*decay.*sin(omegaVec(m)*t + phiRand(m));
            end

            % --- Add second-order irregular if requested ---
            if secondOrder
                % get the 2nd-order horizontal velocity & acceleration at x=0, z=z_b
                [u2x, ~, ~] = secondOrderVelocity( A_m, phiRand, omegaVec, k_m_vec, ...
                                                   thetaVec, h, 0, z_b, t );
                [a2x, ~, ~] = secondOrderAcceleration( A_m, phiRand, omegaVec, k_m_vec, ...
                                                       thetaVec, h, 0, z_b, t );
                U    = U    + u2x;
                Udot = Udot + a2x;
            end


        % ——— Morison drag + inertia ————————————————
        D    = (z_b <= cutoff)*D1 + (z_b > cutoff)*D2;
        A_cs = pi*(D^2)/4;
        Fd   = 0.5*rhoW*Cd*D.*abs(U).*U;
        Fi   = rhoW*Cm*A_cs.*Udot;

        Fm_wave(2*(i-1)+1, :) = (Fd + Fi).*Le;
    end

    % ——— Now get η² at z=0 using our helper —————————
    if secondOrder
        eta2 = secondOrderElevation(A_m, phiRand, omegaVec, k_m_vec, h, t);
    else
        eta2 = zeros(size(t));
    end
    end
end


function Fm_air = computeAerodynamicForces(z, cutoff, D1, D2, h, Le,U_ref,z_ref,alpha,Cd_air,rho_air)
  % Computes aerodynamic drag on above-water nodes using power-law wind
  nNode = numel(z);
  Fm_air = zeros(2*nNode,1);
  for i = 1:nNode
    z_b = z(i);
    if z_b <= h, continue; end
    % wind speed profile
    z_air = z_b - h;
    U_wind = U_ref * (z_air/z_ref)^alpha;
    % section properties
    D = (z_b<=cutoff)*D1 + (z_b>cutoff)*D2;
    % aerodynamic drag
    Fd_air = 0.5 * rho_air * Cd_air * D * abs(U_wind) * U_wind;
    Fm_air(2*(i-1)+1) = Fd_air * Le;
  end
end

function [p, state] = pySimpleHysteretic(y, depth, state, useLinear)
% pySimpleHysteretic: 1D lateral p–y spring, linear or nonlinear (API tanh)
%
% [p, state] = pySimpleHysteretic(y, depth, state, useLinear)
%   y         – lateral deflection [m]
%   depth     – depth below mudline [m] (positive)
%   state     – hysteretic state struct (y_prev, p_prev, direction, …)
%   useLinear – true → p = k*y; false → p = pult * tanh(k*y/pult)
%
% State fields updated even in linear mode for compatibility.

    % Soil params
    phi_eff = 35;                 % deg
    gamma   = 10000;              % N/m^3
    k_ref   = 1e6;                % N/m^3 per unit depth

    % Depth‐dependent stiffness & capacity
    k      = k_ref * depth;       
    pult   = gamma * depth * 3 * tan(deg2rad(phi_eff));

    if useLinear
        % --- Linear spring ---
        p = k * y;
    else
        % --- Nonlinear API‐style tanh ---
        p = pult * tanh(k * y / pult);
    end

    % Update state (for unloading logic later if desired)
    state.y_prev   = y;
    state.p_prev   = p;
    state.direction = sign(y);
end

function [m, state] = mrSimpleHysteretic(theta, depth, state, useLinear)
% mrSimpleHysteretic: 1D rotational m–θ spring, linear or nonlinear (API tanh)
%
% [m, state] = mrSimpleHysteretic(theta, depth, state, useLinear)
%   theta     – rotation [rad]
%   depth     – depth below mudline [m] (positive)
%   state     – hysteretic state struct (theta_prev, m_prev, direction, …)
%   useLinear – true → m = k_theta*theta; false → m_ult * tanh(k_theta*theta/m_ult)

    % Soil&pile params
    gamma        = 10000;         % N/m^3
    D_pile       = 7;             % m
    k_theta_ref  = 0;          % N·m/rad per meter depth

    % Depth‐dependent rotational stiffness & capacity
    k_theta = k_theta_ref * depth;
    m_ult   = 0.5 * gamma * depth * D_pile^2;

    if useLinear
        % --- Linear rotational spring ---
        m = k_theta * theta;
    else
        % --- Nonlinear API‐style tanh ---
        m = m_ult * tanh(k_theta * theta / m_ult);
    end

    % Update state
    state.theta_prev = theta;
    state.m_prev     = m;
    state.direction  = sign(theta);
end

% Second–Order Surface Elevation η^{(2)}
function eta2 = secondOrderElevation(A, phi, omega, k, h, t)
    % eta2(t) = sum_{n,m} A(n)*A(m)*[ L+_{nm} cos((ωn+ωm)t+φn+φm)
    %                               + L-_{nm} cos((ωn-ωm)t+φn-φm) ]
    % where L±_{nm} is given by Eq.3.26 in Kim et al. (2014). :contentReference[oaicite:0]{index=0}

    g    = 9.81;
    N    = numel(omega);
    eta2 = zeros(size(t));
    epsk = 1e-6;

    for n = 1:N
        Rn = k(n)*tanh(k(n)*h);               % Eq.3.28 :contentReference[oaicite:1]{index=1}
        for m = 1:N
            Rm = k(m)*tanh(k(m)*h);
            ksum  = k(n) + k(m);
            kdiff = abs(k(n) - k(m));

            % ——— sum-frequency transfer function L+_{nm} ———
            if ksum > epsk
                % numerator D+_{nm}, Eq.3.29 with “+”
                Dp = ( (sqrt(Rn)+sqrt(Rm)) * ...
                       ( sqrt(Rm)*(k(n)^2 - Rn^2) + sqrt(Rn)*(k(m)^2 - Rm^2) ) / ...
                       ( (sqrt(Rn)+sqrt(Rm))^2 - ksum*tanh(ksum*h) ) ) ...
                   + ( 2*(sqrt(Rn)+sqrt(Rm))^2 * ( k(n)*k(m) - Rn*Rm ) / ...
                       ( (sqrt(Rn)+sqrt(Rm))^2 - ksum*tanh(ksum*h) ) );
                % L+ per Eq.3.26
                Lp = 0.25 * ( Dp ...
                          - (k(n)*k(m) - Rn*Rm)/sqrt(Rn*Rm) ...
                          + (Rn + Rm) );
                wsum = omega(n) + omega(m);
                eta2 = eta2 + A(n)*A(m) * Lp .* cos(wsum*t + phi(n) + phi(m));
            end

            % ——— difference-frequency transfer function L-_{nm} ———
            if kdiff > epsk
                Dm = ( (sqrt(Rn)-sqrt(Rm)) * ...
                       ( sqrt(Rm)*(k(n)^2 - Rn^2) - sqrt(Rn)*(k(m)^2 - Rm^2) ) / ...
                       ( (sqrt(Rn)-sqrt(Rm))^2 - kdiff*tanh(kdiff*h) ) ) ...
                   + ( 2*(sqrt(Rn)-sqrt(Rm))^2 * ( k(n)*k(m) + Rn*Rm ) / ...
                       ( (sqrt(Rn)-sqrt(Rm))^2 - kdiff*tanh(kdiff*h) ) );
                Lm = 0.25 * ( Dm ...
                          - (k(n)*k(m) + Rn*Rm)/sqrt(Rn*Rm) ...
                          + (Rn + Rm) );
                wdiff = omega(n) - omega(m);
                eta2 = eta2 + A(n)*A(m) * Lm .* cos(wdiff*t + phi(n) - phi(m));
            end
        end
    end
end

function Bp = computeBplus(n,m,omega,k,theta,h)
    % compute the “+” quadratic transfer function B⁺_{nm} (Longuet–Higgins, Eq.4.38)
    g     = 9.81;
    wsum  = omega(n) + omega(m);
    kn    = k(n);    km   = k(m);
    thn   = theta(n); thm  = theta(m);
    % dot product kn·km = kn*km*cos(θn−θm)
    kn_km = kn*km*cos(thn-thm);
    num   = wsum * ( kn_km - kn*km );
    denom = wsum^2  - g*(kn + km);
    Bp    = (g^2/(omega(n)*omega(m))) * (num/denom);
end
function Bm = computeBminus(n,m,omega,k,theta,h)
    % compute the “–” quadratic transfer function B⁻_{nm} (Longuet–Higgins, Eq.4.38)
    g      = 9.81;
    wdiff  = omega(n) - omega(m);
    kn     = k(n);    km    = k(m);
    thn    = theta(n); thm   = theta(m);
    kn_km  = kn*km*cos(thn-thm);
    num    = wdiff * ( kn_km + kn*km );
    denom  = wdiff^2  - g*abs(kn - km);
    Bm     = (g^2/(omega(n)*omega(m))) * (num/denom);
end
function [u2x,u2y,u2z] = secondOrderVelocity(A,phi,omega,k,theta,h,x,z,t)
    % second‐order velocity u^{(2)} (Eqs.4.19–4.24), unidirectional
    N     = numel(omega);
    u2x   = zeros(size(t));  u2y = u2x;  u2z = u2x;
    epsw  = 1e-6;
    for n = 1:N
        for m = 1:N
            % sum‐frequency
            Bp   = computeBplus(n,m,omega,k,theta,h);
            wsum = omega(n)+omega(m);
            if abs(Bp)>epsw
                xUp = Bp*(k(n)*cos(theta(n))+k(m)*cos(theta(m)));
                yUp = Bp*(k(n)*sin(theta(n))+k(m)*sin(theta(m)));
                ksum= abs(k(n)*exp(1i*theta(n))+k(m)*exp(1i*theta(m)));
                zUp = 1i*Bp*ksum*tanh(ksum*(h+z));
                phaseP = cos(wsum*t + phi(n)+phi(m));
                u2x = u2x + A(n)*A(m)* xUp .* phaseP;
                u2y = u2y + A(n)*A(m)* yUp .* phaseP;
                u2z = u2z + real(A(n)*A(m)* zUp .* phaseP);
            end
            % diff‐frequency
            Bm    = computeBminus(n,m,omega,k,theta,h);
            wdiff = omega(n)-omega(m);
            if abs(Bm)>epsw
                xUm = Bm*(k(n)*cos(theta(n))-k(m)*cos(theta(m)));
                yUm = Bm*(k(n)*sin(theta(n))-k(m)*sin(theta(m)));
                kdiff = abs(k(n)*exp(1i*theta(n)) - k(m)*exp(1i*theta(m)));
                zUm = 1i*Bm*kdiff*tanh(kdiff*(h+z));
                phaseM = cos(wdiff*t + phi(n)-phi(m));
                u2x = u2x + A(n)*A(m)* xUm .* phaseM;
                u2y = u2y + A(n)*A(m)* yUm .* phaseM;
                u2z = u2z + real(A(n)*A(m)* zUm .* phaseM);
            end
        end
    end
end
function [a2x,a2y,a2z] = secondOrderAcceleration(A,phi,omega,k,theta,h,x,z,t)
    % second‐order acceleration a^{(2)} via analytic d/dt of u2
    N      = numel(omega);
    a2x    = zeros(size(t));  a2y = a2x;  a2z = a2x;
    epsw   = 1e-6;
    for n = 1:N
        for m = 1:N
            % sum‐frequency
            wsum = omega(n)+omega(m);
            if abs(wsum)>epsw
                Bp   = computeBplus(n,m,omega,k,theta,h);
                xUp = Bp*(k(n)*cos(theta(n))+k(m)*cos(theta(m)));
                yUp = Bp*(k(n)*sin(theta(n))+k(m)*sin(theta(m)));
                ksum= abs(k(n)*exp(1i*theta(n))+k(m)*exp(1i*theta(m)));
                zUp = 1i*Bp*ksum*tanh(ksum*(h+z));
                phaseP = sin(wsum*t + phi(n)+phi(m));
                a2x = a2x - A(n)*A(m)* xUp .* (wsum.*phaseP);
                a2y = a2y - A(n)*A(m)* yUp .* (wsum.*phaseP);
                a2z = a2z - real(A(n)*A(m)* zUp .* (wsum.*phaseP));
            end
            % diff‐frequency
            wdiff = omega(n)-omega(m);
            if abs(wdiff)>epsw
                Bm    = computeBminus(n,m,omega,k,theta,h);
                xUm = Bm*(k(n)*cos(theta(n))-k(m)*cos(theta(m)));
                yUm = Bm*(k(n)*sin(theta(n))-k(m)*sin(theta(m)));
                kdiff= abs(k(n)*exp(1i*theta(n)) - k(m)*exp(1i*theta(m)));
                zUm = 1i*Bm*kdiff*tanh(kdiff*(h+z));
                phaseM = sin(wdiff*t + phi(n)-phi(m));
                a2x = a2x - A(n)*A(m)* xUm .* (wdiff.*phaseM);
                a2y = a2y - A(n)*A(m)* yUm .* (wdiff.*phaseM);
                a2z = a2z - real(A(n)*A(m)* zUm .* (wdiff.*phaseM));
            end
        end
    end
end


