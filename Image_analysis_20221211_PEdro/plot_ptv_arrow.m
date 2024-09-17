function plot_ptv_arrow(folder_postions,image_stack_position,image_mask_position,fps,pixel_size,arrow_scale,window_size, scale_factor)

%FROM_PTV_TO_ARROWS  Convert PTV particle positions to velocity arrows and plot average velocities in each window
%
%   Inputs:
%   - folder_postions: string with the path to the folder containing the PTV particle positions
%   - image_stack_position: string with the path to the folder containing the image stack
%   - image_mask_position: string with the path to the folder containing the image mask
%   - fps: scalar with the frames per second of the video
%   - pixel_size: scalar with the pixel size in micrometers
%   - arrow_scale: scalar with the scaling factor for the arrow magnitude
%   - window_size: scalar with the window size in pixels
%   - scale_factor: scalar with the scaling factor for the quiver plot arrow length
%
%   Outputs:
%   - None (outputs are displayed in figure windows)
%
%   Example:
%   from_ptv_to_arrows('C:\PTV_positions', 'C:\image_stack', 'C:\image_mask', 500, 0.32, 1, 10, 10);



%% load positions
fileName = 'positions_final.mat';
load(fullfile(folder_postions, fileName));


%% save the velocity arrows into a matrix

directory = image_stack_position;
file_list = dir([directory, '*.tif']);
im = imread([directory, file_list(1).name]);

% Convert pixel size to millimeters
pixel_size_mm  = pixel_size / 1000; % mm

% Initialize matrix for storing positions and velocities
pos_vel_matrix = zeros(numel(particlepositionfinale_changes), 5);

% Iterate over all particles
for i = 1:numel(particlepositionfinale_changes)
    if size(particlepositionfinale_changes{1,i},1) > 3
        % Extract x-y coordinates for current particle
        pos = particlepositionfinale_changes{i}(:,2:3);

        % Compute time interval between consecutive frames
        dt = 1/fps; % seconds
        t = (0:size(pos,1)-2)*dt;

        % Compute velocity
        vel = diff(pos)./dt;
        vel_mag = sqrt(sum((vel .* repmat([pixel_size_mm pixel_size_mm], size(vel,1),1)).^2,2)); % mm/s

        % Compute arrow positions and magnitudes
        arrow_pos = pos(1:end-1,:);
        arrow_mag = vel_mag.*arrow_scale;

        % Compute arrow directions
        arrow_dir = vel./repmat(sqrt(sum(vel.^2,2)),1,2);
        
        nan_mask = isnan(arrow_dir);
        arrow_dir(nan_mask) = zeros(sum(nan_mask(:)),1);
        
        % Check for NaNs in the arrow_dir matrix
        nan_idx = find(any(isnan(arrow_dir), 2));
        if ~isempty(nan_idx)
            disp(['Found NaNs in arrow_dir matrix for particle ', num2str(i), ' at indices:']);
            disp(nan_idx);
        end

        % Add positions and velocities to pos_vel_matrix
        n_arrows = size(arrow_pos, 1);
        pos_vel_matrix((i-1)*n_arrows+1:i*n_arrows,1:2) = arrow_pos;
        pos_vel_matrix((i-1)*n_arrows+1:i*n_arrows,3:4) = arrow_mag.*arrow_dir;
        pos_vel_matrix((i-1)*n_arrows+1:i*n_arrows,5) = vel_mag;
    end
end

disp('%%%%%%%%%%%% saving the velocity arrows into a matrix %%%%%%%%%%%%')

% imshow(im);
% hold on;
% quiver(pos_vel_matrix(:,1), pos_vel_matrix(:,2), pos_vel_matrix(:,3), pos_vel_matrix(:,4), 0, 'color', 'k', 'linewidth', 1, 'MaxHeadSize', 0.8);


%% calculate the arrows

% Load the image
directory = image_stack_position;
file_list = dir([directory, '*.tif']);
im = imread([directory, file_list(1).name]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% coment if no mask %%%%%%%%%%%%%%%%%%%%%%%%%%%
% directory = image_mask_position;
% file_list = dir([directory, '*.tiff']);
% binary_mask = imread([directory, file_list(1).name]);
% 
% % Convert the binary mask to a logical array
% binary_mask = logical(binary_mask);
% 
% % Set the grayscale value to 255 (white) where the binary mask is true
% im(binary_mask) = 120;



% Compute the number of rows and columns in the grid
[im_height, im_width, ~] = size(im);
n_rows = ceil(im_height / window_size);
n_cols = ceil(im_width / window_size);

% Preallocate the average velocities matrix
avg_velocities = nan(n_rows, n_cols, 2);

% Preallocate the standard deviation matrix
std_velocities = nan(n_rows, n_cols, 2);

% Preallocate the magnitude standard deviation matrix
std_magnitude = nan(n_rows, n_cols);

% Compute the window indices for each position
window_indices = ceil(pos_vel_matrix(:, 1:2) / window_size);

% Calculate the total number of windows
total_windows = n_rows * n_cols;

% Calculate the average velocities for all windows at once
for k = 1:total_windows
    [r, c] = ind2sub([n_rows, n_cols], k);
    in_window = window_indices(:, 1) == c & window_indices(:, 2) == r;
    
    if any(in_window)
        avg_velocities(r, c, :) = mean(pos_vel_matrix(in_window, 3:4));
        std_velocities(r, c, :) = std(pos_vel_matrix(in_window, 3:4));
        std_magnitude(r, c) = std(pos_vel_matrix(in_window, 5));
    end
    if mod(k, 10) == 0 % check if i is a multiple of 1000
        progress = k / total_windows * 100; % calculate percentage of completion
        fprintf('%.2f%% complete\n', progress); % print percentage with 2 decimal places
    end
end

% Create a meshgrid for the centers of the windows
[X, Y] = meshgrid((1:n_cols) * window_size - window_size/2, (1:n_rows) * window_size - window_size/2);

disp('%%%%%%%%%%%% Plot the window arrow velocity  %%%%%%%%%%%%')

%% plot the images
% Plot the image and the arrows representing the average velocity in each window

% %%%%%%%%%%%%%%%%%%%%%%%%%%% coment if no mask %%%%%%%%%%%%%%%%%%%%%%%%%%%
% directory = image_mask_position;
% file_list = dir([directory, '*.tiff']);
% binary_mask = imread([directory, file_list(1).name]);
% 
% % Convert the binary mask to a logical array
% binary_mask = logical(binary_mask);
% 
% % Set the grayscale value to 255 (white) where the binary mask is true
% im(binary_mask) = 120;
% %%%%%%%%%%%%%%%%%%%%%%%%%%% coment if no mask %%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plot the grayscale image and allow the user to select a rectangular box
im_rgb = cat(3, im, im, im); % convert grayscale image to RGB format

figure(1);
imshow(im_rgb);
hold on
colormap(gray);


%% plot the trajectories
h_rect = imrect; % create an interactive rectangle on the image
setResizable(h_rect, true); % allow the user to resize the rectangle

% Wait for the user to finish selecting the rectangular box
wait(h_rect);
rect_pos = getPosition(h_rect); % get the position of the rectangle

% Define the rectangular box to crop the trajectories
x_min = rect_pos(1); % left edge of the box
y_min = rect_pos(2); % top edge of the box
width = rect_pos(3); % width of the box
height = rect_pos(4); % height of the box

% Define the desired colormap based on the velocity values you want to use
velocity_range = [0 15]; % define the range of velocity values to use
cmap = colormap(jet(256)); % define a colormap (using the jet colormap as an example)
cmap = interp1(linspace(velocity_range(1),velocity_range(2),size(cmap,1)), cmap, linspace(velocity_range(1),velocity_range(2),256));


% Iterate over all particles
hold on;
for i = 1:numel(particlepositionfinale_changes)
    if size(particlepositionfinale_changes{1,i},1)>10
        % Extract x-y coordinates for current particle
        pos = particlepositionfinale_changes{i}(:,2:3);
        
        % Check if any part of the trajectory is within the rectangular box
        if any(pos(:,1) > x_min & pos(:,1) < x_min + width & pos(:,2) > y_min & pos(:,2) < y_min + height)
            
            % Crop the trajectories to the rectangular box
            idx = pos(:,1) > x_min & pos(:,1) < x_min + width & pos(:,2) > y_min & pos(:,2) < y_min + height;
            pos_crop = pos(idx,:);
            
            % Check if there are at least two points within the rectangular box
            if size(pos_crop, 1) >= 2
                % Compute time interval between consecutive frames
                dt = 1/fps; % seconds
                t = (0:size(pos_crop,1)-2)*dt;

                % Compute velocity
                vel = diff(pos_crop)./dt;
                vel_mag = sqrt(sum((vel .* repmat([pixel_size_mm pixel_size_mm], size(vel,1),1)).^2,2)); % mm/s

                % Remove any NaN values from vel_mag and clip to the specified velocity range
                vel_mag(isnan(vel_mag)) = 0;
                vel_mag(vel_mag < velocity_range(1)) = velocity_range(1);
                vel_mag(vel_mag > velocity_range(2)) = velocity_range(2);

                % Define the color values for the trajectories based on velocity magnitude
                vel_color = interp1(linspace(velocity_range(1),velocity_range(2),size(cmap,1)), cmap, vel_mag);

                % Plot the trajectories with colors based on velocity magnitude
                for j = 1:size(pos_crop,1)-1
                    x = [pos_crop(j,1) pos_crop(j+1,1)];
                    y = [pos_crop(j,2) pos_crop(j+1,2)];
                    line(x, y, 'color', vel_color(j,:), 'linewidth', 3);
                end
            end
        end
    end
     if mod(i, 10) == 0 % check if i is a multiple of 1000
        progress = i / numel(particlepositionfinale_changes) * 100; % calculate percentage of completion
        fprintf('%.2f%% complete\n', progress); % print percentage with 2 decimal places
    end
end

% Set the colormap and color limits for the colorbar
colormap(cmap);
caxis(velocity_range);

% Flip the y-axis
set(gca, 'YDir', 'reverse');

% Add colorbar for velocity magnitude
c = colorbar;
c.Label.String = 'Velocity magnitude (mm/s)';

% Add title and axes labels
title('Velocity magnitude - 400um, 20x,500fps, 17KHz, 2.0V, MidPlane');
xlabel('X (pixels)');
ylabel('Y (pixels)');


%% plot the arrows
% Find the top 10 maximum average velocities
flat_avg_velocities = sqrt(avg_velocities(:, :, 1).^2 + avg_velocities(:, :, 2).^2);
[sorted_velocities, sorted_indices] = sort(flat_avg_velocities(:), 'descend');
top_10_indices = sorted_indices(1:10);
top_10_values = sorted_velocities(1:10);

% Create a logical matrix indicating the top 10 maximum average velocity arrows
top_10_logical = false(size(flat_avg_velocities));
top_10_logical(top_10_indices) = true;

% Display the top 10 maximum average velocities in the command line
fprintf('Top 10 Maximum Average Velocities:\n');
for i = 1:10
    fprintf('%d. %f\n', i, top_10_values(i));
end

%%%%%%%%%%%%%%%%%

% Plot the green arrows representing the average velocity in each window
green_avg_velocities = avg_velocities;
green_avg_velocities(repmat(top_10_logical, [1 1 2])) = 0;
quiver(X, Y, scale_factor * green_avg_velocities(:, :, 1), scale_factor * green_avg_velocities(:, :, 2), 'color', 'w', 'linewidth', 2, 'MaxHeadSize', 5, 'AutoScale', 'off');
hold on;

% Plot the top 10 maximum average velocity arrows in red
red_avg_velocities = zeros(size(avg_velocities));
red_avg_velocities(repmat(top_10_logical, [1 1 2])) = avg_velocities(repmat(top_10_logical, [1 1 2]));
quiver(X, Y, scale_factor * red_avg_velocities(:, :, 1), scale_factor * red_avg_velocities(:, :, 2), 'color', 'w', 'linewidth', 2, 'MaxHeadSize', 5, 'AutoScale', 'off');

% Calculate a representative arrow length based on the representative velocity value
arrow_length = scale_factor * 1;

% Add a representative green arrow to the same figure as the quiver plot
hold on;
quiver(gca, 0.9*im_width, 0.5*im_height, arrow_length, 0, 'Color', 'w', 'LineWidth', 2, 'MaxHeadSize', 5, 'AutoScale', 'off');
axis(gca, 'equal');

% Add a label with the respective velocity
text(gca, 0.8*im_width + arrow_length + 10, 0.5*im_height, sprintf('Velocity: %.2f mm/s', 1), 'FontSize', 12);


disp('%%%%%%%%%%%% ploting the images  %%%%%%%%%%%%')

hold off
end