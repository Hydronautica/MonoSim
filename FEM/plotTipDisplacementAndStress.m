function plotTipDisplacementAndStress(time, U, elemZ, params)
% PLOTTIPDISPLACEMENTANDSTRESS Plots tip displacement and bending stress
%
% Inputs:
%   time   - Time vector
%   U      - Displacement array
%   elemZ  - Z-coordinates of nodes
%   params - Structure containing simulation parameters

% Extract parameters
sectionCutOff = params.sectionCutOff;
D_outer1 = params.D_outer1;
D_outer2 = params.D_outer2;
E1 = params.E1;
E2 = params.E2;
z_stress = params.z_stress;
topDOF = 2*(params.nNode-1) + 1;
Le = params.Le;

% Create figure
figure;

% Plot tip displacement
subplot(2,1,1);
plot(time, U(topDOF, :), 'LineWidth', 4);
xlabel('Time [s]'); ylabel('Tip disp [m]'); grid on;
grid minor;
[max_disp, max_idx] = max(U(topDOF, :));
fprintf('Maximum tip displacement: %.4f m at t = %.2f s\n', max_disp, time(max_idx));

% Find the nearest element index for each z_stress location
nLocs = numel(z_stress);
stress_plot = arrayfun(@(z) ...
    find(abs(elemZ - z) == min(abs(elemZ - z)), 1), ...
    z_stress);

% Preallocate: each row is one location, each column one time step
stress_time = zeros(nLocs, length(time));

% Loop over each location
for iLoc = 1:nLocs
    node_loc = stress_plot(iLoc);
    
    % Pick material properties based on cutoff
    if z_stress(iLoc) <= sectionCutOff
        D_loc = D_outer1; E_loc = E1;
    else
        D_loc = D_outer2; E_loc = E2;
    end
    c = D_loc/2;
    
    % Compute stress time history at this node
    for j = 1:length(time)
        % Degrees of freedom: i = this node, p = previous, n = next
        idx_i = 2*(node_loc-1) + 1;
        idx_p = 2*(node_loc-2) + 1;
        idx_n = 2*(node_loc) + 1;
        
        % Handle boundary cases
        if node_loc == 1
            idx_p = idx_i;
        end
        if node_loc == length(elemZ)
            idx_n = idx_i;
        end
        
        % Calculate curvature using central difference
        curvature = (U(idx_p, j) - 2*U(idx_i, j) + U(idx_n, j)) / Le^2;
        
        % Calculate bending stress
        stress_time(iLoc, j) = E_loc * curvature * c;
    end
end

% Plot all stress locations on one figure
subplot(2,1,2);
hold on; grid on; grid minor;
colors = lines(nLocs);
for iLoc = 1:nLocs
    plot(time, stress_time(iLoc,:), ...
         'LineWidth', 5, ...
         'Color', colors(iLoc,:));
end
xlabel('Time [s]');
ylabel('Bending stress [Pa]');
legend(arrayfun(@(z) sprintf('z = %.1f m', z), z_stress, ...
               'UniformOutput', false), ...
       'Location', 'Best');

end