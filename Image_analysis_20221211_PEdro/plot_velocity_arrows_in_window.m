function plot_velocity_arrows_in_window(im, particle_positions, fps, pixel_size, arrow_scale, window_center)
% PLOT_VELOCITY_ARROWS_IN_WINDOW - plots velocity magnitudes with arrows in a specified window
%
% Syntax:
%   plot_velocity_arrows_in_window(im, particle_positions, fps, pixel_size, arrow_scale, window_center)
%
% Inputs:
%   im - grayscale image
%   particle_positions - a cell array of particle positions over time
%   fps - frames per second
%   pixel_size - pixel size in micrometers
%   arrow_scale - a scalar defining the arrow scale factor
%   window_center - a 2-element vector specifying the center of the window in pixels
%
% Outputs:
%   None (plots the grayscale image and velocity arrows in the specified window)

% Convert pixel size to millimeters
pixel_size_mm = pixel_size / 1000; % mm

% Define window size and extract particles within window
window_size = 10; % pixels
window_half_size = floor(window_size / 2);
window_top_left = window_center - window_half_size;
window_bottom_right = window_center + window_half_size;
particle_positions_window = {};
for i = 1:numel(particle_positions)
    particle_positions_i = particle_positions{i}(:,2:3);
    in_window = particle_positions_i(:,1) >= window_top_left(1) & ...
                particle_positions_i(:,1) <= window_bottom_right(1) & ...
                particle_positions_i(:,2) >= window_top_left(2) & ...
                particle_positions_i(:,2) <= window_bottom_right(2);
    if any(in_window)
        particle_positions_window{end+1} = particle_positions{i}(in_window,:);
    end
end

% Compute average velocity in window
vel_window = [];
for i = 1:numel(particle_positions_window)
    pos = particle_positions_window{i}(:,2:3);
    dt = 1/fps;
    t = (0:size(pos,1)-2)*dt;
    vel = diff(pos)./dt;
    vel_mag = sqrt(sum((vel .* repmat([pixel_size_mm pixel_size_mm], size(vel,1),1)).^2,2));
    vel_mag(isnan(vel_mag)) = 0;
    vel_window = [vel_window; vel];
end
if isempty(vel_window)
    warning('No particles found in window');
    return;
end
avg_vel = mean(vel_window, 1);

% Plot the grayscale image
imshow(im);
hold on;

% Plot the velocity arrows in window
if ~isempty(avg_vel)
    arrow_pos = window_center;
    arrow_mag = mean(vel_mag) * arrow_scale;
    arrow_dir = avg_vel / norm(avg_vel);
    quiver(arrow_pos(1), arrow_pos(2), arrow_mag*arrow_dir(1), arrow_mag*arrow_dir(2), 0, 'color', 'r', 'linewidth', 1, 'MaxHeadSize', 0.8);
end

% Flip the y-axis back to its normal orientation
set(gca, 'YDir', 'reverse');

% Add title and axes labels
title('Velocity magnitudes with arrows');
xlabel('X (pixels)');
ylabel('Y (pixels)');

end
