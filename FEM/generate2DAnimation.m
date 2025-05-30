function generate2DAnimation(time, U, FHyd_history, FAero_history, elemZ, params)
% GENERATE2DANIMATION Creates a 2D animation of the monopile with forces
%
% Inputs:
%   time          - Time vector
%   U             - Displacement array
%   FHyd_history  - History of hydrodynamic forces
%   FAero_history - History of aerodynamic forces
%   elemZ         - Z-coordinates of nodes
%   params        - Structure containing simulation parameters

% Extract parameters
d_factor = params.d_factor;
Le = params.Le;
L_total = params.L_total;
nSteps = params.nSteps;

% Create video writer object
v = VideoWriter('morisonLoads_deflected.mp4', 'MPEG-4');
v.FrameRate = 10;
open(v);

figure('Color', 'w');
for i = 1:10:nSteps+1
    clf; hold on; grid on;
    
    % Plot deflected tower (scaled)
    U_def = d_factor * U(1:2:end, i);  % U at each node, scaled
    plot(U_def, elemZ, 'k-', 'LineWidth', 2);
    
    % Plot morison loads as arrows at the deflected tower
    scale = Le*2 / max(abs([FHyd_history(:); FAero_history(:)]));
    for j = 1:length(elemZ)
        Fji = FHyd_history(2*(j-1)+1, i) + FAero_history(2*(j-1)+1, i);
        quiver(U_def(j), elemZ(j), Fji*scale, 0, 0, 'r', ...
               'MaxHeadSize', 3, 'LineWidth', 1);
    end
    
    title(sprintf('Loads & Deflection (×%d) at t = %.2f s', d_factor, time(i)));
    xlabel('Horizontal [m]'); ylabel('Elevation [m]');
    axis tight; ylim([0 L_total+20]); xlim([-2 2]);
    
    % Capture frame and write to video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
disp('2D Animation saved to morisonLoads_deflected.mp4');

end