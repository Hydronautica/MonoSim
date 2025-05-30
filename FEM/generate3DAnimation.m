function generate3DAnimation(time, U, FHyd_history, FAero_history, Fsoil_history, elemZ, params)
% GENERATE3DANIMATION Creates a 3D animation of the monopile with wave surface
%
% Inputs:
%   time          - Time vector
%   U             - Displacement array
%   FHyd_history  - History of hydrodynamic forces
%   FAero_history - History of aerodynamic forces
%   Fsoil_history - History of soil reaction forces
%   elemZ         - Z-coordinates of nodes
%   params        - Structure containing simulation parameters

% Extract parameters
d_factor = params.d_factor;
D_outer1 = params.D_outer1;
D_outer2 = params.D_outer2;
sectionCutOff = params.sectionCutOff;
h = params.h;
nSteps = params.nSteps;
Nfreq = params.Nfreq;
omegaVec = params.omegaVec;
A_m = params.A_m;
phiRand = params.phiRand;
k_m_vec = params.k_m_vec;
guywires = params.guywires;
j_guy = params.j_guy;
dof_guy = params.dof_guy;
z_nodeG1 = params.z_nodeG1;
anchorR_x = params.anchorR_x;
anchorR_z = params.anchorR_z;
anchorL_x = params.anchorL_x;
anchorL_z = params.anchorL_z;

% Pre-compute geometry that doesn't change
nCirc = 36;
theta = linspace(0, 2*pi, nCirc)';
[Zb, Th] = meshgrid(elemZ, theta);
Rb = zeros(size(Zb));
for j = 1:length(elemZ)
    Rb(:,j) = (elemZ(j) <= sectionCutOff) * (D_outer1/2) + (elemZ(j) > sectionCutOff) * (D_outer2/2);
end

% Set up video
v = VideoWriter('beam_wave_animation.mp4', 'MPEG-4');
v.FrameRate = 30;
open(v);

% Set stress color limits
smin = -5e8;
smax = 5e8;

% Set force arrow scale factor
scaleFactor = 1e-4;  % Adjust to get reasonable arrow length

% Create figure
figWidth = 720;   % Desired width in pixels
figHeight = 1080;  % Desired height in pixels
figure('Color', 'w', 'Position', [100, 100, figWidth, figHeight], ...
    'MenuBar', 'none', 'ToolBar', 'none');

% Animation loop
for k = floor(nSteps/2):10:nSteps+1
    t0 = time(k);
    
    % 1) Beam deflection & stress at time k
    Uk = U(1:2:end, k);  % nodal vertical disp
    
    % Recompute stress at each node
    stress_k = zeros(length(elemZ), 1);
    for j = 2:length(elemZ)-1
        kappa = (Uk(j-1) - 2*Uk(j) + Uk(j+1)) / (elemZ(2)-elemZ(1))^2;
        if elemZ(j) <= sectionCutOff
            Ej = params.E1; Dj = D_outer1;
        else
            Ej = params.E2; Dj = D_outer2;
        end
        cj = Dj/2;
        stress_k(j) = Ej * cj * kappa;
    end
    stress_k([1,end]) = stress_k(2)*0.5 + stress_k(end-1)*0.5;  % endpoints
    
    % Expand around circumference
    Uk_mat = d_factor .* repmat(Uk', nCirc, 1);
    stress_mat = repmat(stress_k', nCirc, 1);
    Xb = Rb .* cos(Th) + Uk_mat;
    Yb = Rb .* sin(Th);
    
    % 2) Wave surface at time k
    xw = linspace(-50, 50, 200);
    yw = linspace(-50, 50, 200);
    [Xw, Yw] = meshgrid(xw, yw);
    eta_xy = zeros(size(Xw));
    for m = 1:Nfreq
        phase = k_m_vec(m)*Xw - omegaVec(m)*t0 + phiRand(m);
        eta_xy = eta_xy + A_m(m)*cos(phase);
    end
    
    % === 3) PLOT ===
    clf; hold on; grid on;
    
    % --- Plot monopile deformation + stress ---
    surf(Xb, Yb, Zb, stress_mat, 'FaceColor', 'interp', 'EdgeColor', 'none');
    
    % --- Transparent mudline plane at z = 0 ---
    [Xm, Ym] = meshgrid(linspace(-50, 50, 2), linspace(-50, 50, 2));
    Zm = zeros(size(Xm));
    surf(Xm, Ym, Zm, 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % --- Plot wave surface ---
    surf(Xw, Yw, eta_xy + h, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    % --- Soil restoring force arrows below mudline ---
    Zq = elemZ(:);
    Xq = zeros(size(Zq));
    Yq = zeros(size(Zq));
    Fq = Fsoil_history(1:2:end, k);
    
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
    if guywires
        % Get current pile head location
        x_tip = U(dof_guy, k);   % horizontal displacement of guy node
        z_tip = z_nodeG1;        % vertical coordinate of guy node (fixed)
        
        % Right guy wire
        plot3([x_tip, anchorR_x], [0, 0], [z_tip, anchorR_z], ...
              'r-', 'LineWidth', 2, 'DisplayName', 'Right Guy Wire');
        
        % Left guy wire
        plot3([x_tip, anchorL_x], [0, 0], [z_tip, anchorL_z], ...
              'b-', 'LineWidth', 2, 'DisplayName', 'Left Guy Wire');
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
    
    % Write frame to video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
disp('3D Animation saved to beam_wave_animation.mp4');

end