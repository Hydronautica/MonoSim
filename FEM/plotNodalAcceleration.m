function plotNodalAcceleration(time, Aacc, elemZ, params)
% PLOTNODALACCELERATION Plots acceleration time series at specified locations
%
% Inputs:
%   time   - Time vector
%   Aacc   - Acceleration array
%   elemZ  - Z-coordinates of nodes
%   params - Structure containing simulation parameters

% Extract parameters
z_acc = params.z_acc;

% Find the nearest node index for each z_acc location
acc_plot = arrayfun(@(z) ...
    find(abs(elemZ - z) == min(abs(elemZ - z)), 1), ...
    z_acc);
nLocs = numel(z_acc);

% Preallocate acceleration time series array
accel_time = zeros(nLocs, length(time));

% Extract acceleration time series for each location
for iLoc = 1:nLocs
    node_loc = acc_plot(iLoc);
    dof_loc = 2*(node_loc-1) + 1;  % Horizontal DOF
    accel_time(iLoc, :) = Aacc(dof_loc, :);
end

% Plot acceleration at each location
figure; hold on; grid on;
colors = lines(nLocs);
for iLoc = 1:nLocs
    plot(time, accel_time(iLoc,:), ...
         'LineWidth', 1.5, ...
         'Color', colors(iLoc,:));
end
xlabel('Time [s]');
ylabel('Nodal acceleration [m/s^2]');
legend(arrayfun(@(z) sprintf('z = %.1f m', z), z_acc, ...
               'UniformOutput', false), ...
       'Location', 'Best');
title('Nodal Acceleration at Selected Locations');

end