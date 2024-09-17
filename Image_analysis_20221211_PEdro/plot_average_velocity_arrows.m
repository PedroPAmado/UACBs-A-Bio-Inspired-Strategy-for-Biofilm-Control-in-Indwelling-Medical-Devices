function plot_average_velocity_arrows(im, particle_positions, fps, pixel_size, arrow_scale, window_size, overlap)
% PLOT_AVERAGE_VELOCITY_ARROWS - plots average velocity magnitudes with arrows in small windows
%
% Syntax:
%   plot_average_velocity_arrows(im, particle_positions, fps, pixel_size, arrow_scale, window_size, overlap)
%
% Inputs:
%   im - grayscale image
%   particle_positions - a cell array of particle positions over time
%   fps - frames per second
%   pixel_size - pixel size in micrometers
%   arrow_scale - a scalar defining the arrow scale factor
%   window_size - size of the windows for computing average velocity (in pixels)
%   overlap - overlap between adjacent windows (in pixels)
%
% Outputs:
%   None (plots the grayscale image and velocity arrows)

% Convert pixel size to millimeters
pixel_size_mm = pixel_size / 1000; % mm

% Get the image dimensions
im_height = size(im, 1);
im_width = size(im, 2);

% Initialize the array for storing the average velocities
avg_velocities = zeros(im_height/window_size, im_width/window_size, 2);

% Loop over all the windows in the image
for y = 1:overlap:im_height-window_size+1
    for x = 1:overlap:im_width-window_size+1
        
        % Get the particle positions that fall within the current window
        particle_positions_in_window = {};
        for i = 1:numel(particle_positions)
            pos = particle_positions{i}(:,2:3);
            indices = pos(:,1) >= x & pos(:,1) <= x+window_size-1 & pos(:,2) >= y & pos(:,2) <= y+window_size-1;
            if any(indices)
                particle_positions_in_window{end+1} = pos(indices,:);
            end
        end
        
        % Compute the velocities for all the particles within the window
        all_velocities = []; % initialize the array for storing the velocities
        for i = 1:numel(particle_positions_in_window)
            pos = particle_positions_in_window{i};
            dt = 1/fps;
            t = (0:size(pos,1)-2)*dt;
            vel = diff(pos)./dt;
            vel_mag = sqrt(sum((vel .* repmat([pixel_size_mm pixel_size_mm], size(vel,1),1)).^2,2));
            vel_mag(isnan(vel_mag)) = 0;
            if ~isempty(vel)
                all_velocities = [all_velocities; vel];
            end
        end
        
        % Compute the average velocity for all the particles within the window
        if ~isempty(all_velocities)
            avg_velocities(y+window_size/2-1, x+window_size/2-1,:) = mean(all_velocities, 1);
        end
        
    end
end


% Plot the average velocities as vectors at the center of each window
imshow(im);
hold on;
for y = window_size/2:overlap:im_height-window_size/2+1
    for x = window_size/2:overlap:im_width-window_size/2+1
        if any(avg_velocities(y,x,:))
            quiver(x, y, avg_velocities(y,x,1)*arrow_scale, avg_velocities(y,x,2)*arrow_scale, 0, 'color', 'k', 'linewidth', 1, 'MaxHeadSize', 0.8);
        end
    end
end

% Flip the y-axis back to its normal orientation
set(gca, 'YDir', 'reverse');

% Add title and axes labels
title('Average velocity magnitudes with arrows in small windows');
xlabel('X (pixels)');
ylabel('Y (pixels)');

end
