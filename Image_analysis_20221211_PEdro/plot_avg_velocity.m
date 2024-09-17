function plot_avg_velocity_2(im, particle_positions, fps, pixel_size, window_size)
% PLOT_AVG_VELOCITY_2 - plots the average velocity magnitude in a sliding window for all trajectories
%
% Syntax:
%   plot_avg_velocity_2(im, particle_positions, fps, pixel_size, window_size)
%
% Inputs:
%   im - grayscale image
%   particle_positions - a cell array of particle positions over time
%   fps - frames per second
%   pixel_size - pixel size in micrometers
%   window_size - size of the sliding window in pixels
%
% Outputs:
%   None (plots the grayscale image with velocity magnitude color-coded)

% Convert pixel size to millimeters
pixel_size_mm = pixel_size / 1000; % mm

% Compute the sliding window half-size in pixels
half_window_size = floor(window_size/2);

% Compute the sliding window half-size in millimeters
half_window_size_mm = half_window_size * pixel_size_mm;

% Compute time interval between consecutive frames
dt = 1/fps; % seconds

% Compute velocity and velocity magnitude for all particles
velocities = cellfun(@(pos) diff(pos(:,2:3)) ./ dt, particle_positions, 'UniformOutput', false);
vel_magnitudes = cellfun(@(vel) sqrt(sum((vel .* repmat([pixel_size_mm pixel_size_mm], size(vel,1),1)).^2,2)), velocities, 'UniformOutput', false);

% Replace NaN values with 0
for i = 1:numel(vel_magnitudes)
    vel_magnitudes{i}(isnan(vel_magnitudes{i})) = 0;
end

% Compute the average velocity magnitude in a sliding window for all particles
avg_vel_mag = zeros(size(im));
for i = 1:numel(particle_positions)
    if size(particle_positions{i},1) > 30
        pos = particle_positions{i}(:,2:3);
        vel_mag = vel_magnitudes{i};
        for j = 1:size(pos,1)
            x_min = max(1, pos(j,1) - half_window_size);
            x_max = min(size(im,2), pos(j,1) + half_window_size);
            y_min = max(1, pos(j,2) - half_window_size);
            y_max = min(size(im,1), pos(j,2) + half_window_size);
            window_vel_mag = vel_mag(pos(:,1) >= x_min & pos(:,1) <= x_max & pos(:,2) >= y_min & pos(:,2) <= y_max);
            avg_vel_mag(y_min:y_max, x_min:x_max) = avg_vel_mag(y_min:y_max, x_min:x_max) + mean(window_vel_mag);
        end
    end
end

% Normalize the average velocity magnitude by the number of contributing particles
particle_counts = cellfun(@(pos) size(pos,1), particle_positions);
particle_counts(particle_counts == 0) = 1; % prevent divide by zero
avg_vel_mag = avg_vel_mag ./ particle_counts;

% Plot the average velocity magnitude color-coded on the grayscale image
imshow(im);
hold on;
imagesc(avg_vel_mag);
colormap(jet);
caxis([0 3]);
alpha(0.5);

% Reverse the y-axis
set(gca, 'YDir', 'reverse');

% Add colorbar for velocity magnitude
c = colorbar;
c.Label.String = 'Average velocity magnitude (mm/s)';

% Add colorbar for velocity magnitude
c = colorbar;
c.Label.String = 'Average velocity magnitude (mm/s)';

end