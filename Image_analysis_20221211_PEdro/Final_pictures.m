%%
clc 
clear all 
close all

%%


Numberparticle=max(particleInfo(:,6));

for ii=1:Numberparticle
    b=find(ii==particleInfo(:,6));
    
    for iii=1:length(b)
    
       particleposition_changes(iii,:)=particleInfo(b(iii),:);
 
    end 
     particlepositionfinale_changes{ii}=particleposition_changes;
 
    clear b particleposition_changes
    
end

%%

directory = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\test_script\Image_visualization\Image_stack\stack.tif';
image = imread(directory);
figure
% Plot the image
imshow(image);
hold on;

for i=1:size( particlepositionfinale_changes,2)
    if size(particlepositionfinale_changes{1,i},1)>3
        a= particlepositionfinale_changes{i};
        plot (a(:,2),a(:,3));
        hold on
        clear a
    end
end
set(gca,'YDir','reverse')


%% with no smothing and set velocity for the range we want
% Assume frame rate is 10 frames per second and corrected pixel size is 0.2 micrometers
fps = 500;
pixel_size = 0.29; % micrometers

% Convert pixel size to millimeters
pixel_size_mm  = pixel_size / 1000; % mm

% Define the desired colormap based on the velocity values you want to use
velocity_range = [0 2]; % define the range of velocity values to use
cmap = colormap(jet(256)); % define a colormap (using the jet colormap as an example)
cmap = interp1(linspace(velocity_range(1),velocity_range(2),size(cmap,1)), cmap, linspace(velocity_range(1),velocity_range(2),256));

% Load the image
directory = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\test_script\Image_visualization\Image_stack/';
file_list = dir([directory, '*.tif']);
im = imread([directory, file_list(1).name]);
im_rgb = cat(3, im, im, im); % convert grayscale image to RGB format

% Plot the grayscale image
figure;
imshow(im_rgb);
colormap(gray);

% Iterate over all particles
hold on;
for i = 1:numel(particlepositionfinale_changes)
    if size(particlepositionfinale_changes{1,i},1)>10
        % Extract x-y coordinates for current particle
        pos = particlepositionfinale_changes{i}(:,2:3);
        
        % Compute time interval between consecutive frames
        dt = 1/fps; % seconds
        t = (0:size(pos,1)-2)*dt;
        
        % Compute velocity
        vel = diff(pos)./dt;
        vel_mag = sqrt(sum((vel .* repmat([pixel_size_mm pixel_size_mm], size(vel,1),1)).^2,2)); % mm/s
        
        % Remove any NaN values from vel_mag and clip to the specified velocity range
        vel_mag(isnan(vel_mag)) = 0;
        vel_mag(vel_mag < velocity_range(1)) = velocity_range(1);
        vel_mag(vel_mag > velocity_range(2)) = velocity_range(2);
        
        % Define the color values for the trajectories based on velocity magnitude
        vel_color = interp1(linspace(velocity_range(1),velocity_range(2),size(cmap,1)), cmap, vel_mag);
        
        % Plot the trajectories with colors based on velocity magnitude
        for j = 1:size(pos,1)-1
            x = [pos(j,1) pos(j+1,1)];
            y = [pos(j,2) pos(j+1,2)];
            line(x, y, 'color', vel_color(j,:), 'linewidth', 3);
        end
        
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


%% using croped zones for ploting the trajectroies
% Assume frame rate is 10 frames per second and corrected pixel size is 0.2 micrometers
fps = 500;
pixel_size = 0.33; % micrometers

% Convert pixel size to millimeters
pixel_size_mm  = pixel_size / 1000; % mm

% Define the desired colormap based on the velocity values you want to use
velocity_range = [0 3]; % define the range of velocity values to use
cmap = colormap(jet(256)); % define a colormap (using the jet colormap as an example)
cmap = interp1(linspace(velocity_range(1),velocity_range(2),size(cmap,1)), cmap, linspace(velocity_range(1),velocity_range(2),256));

% Load the image
directory = 'H:\Cornel_ First Experiments\24.03.2023\vid_2023-03-24_15-37-16_53/';
file_list = dir([directory, '*.tiff']);
im = imread([directory, file_list(1).name]);
im_rgb = cat(3, im, im, im); % convert grayscale image to RGB format

% Plot the grayscale image and allow the user to select a rectangular box
figure;
imshow(im_rgb);
colormap(gray);
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


% Iterate over all particles
hold on;
for i = 1:numel(particlepositionfinale_changes)
    if size(particlepositionfinale_changes{1,i},1)>6
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
                    line(x, y, 'color', vel_color(j,:), 'linewidth', 2);
                end
            end
        end
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
title('Trajectories with Colors Based on Velocity Magnitude');
xlabel('X (pixels)');
ylabel('Y (pixels)');

            
%% change the image i want to use

% Load the grayscale image and binary mask
directory = 'H:\Cornel_ First Experiments\24.03.2023\vid_2023-03-24_15-09-44_44/';
file_list = dir([directory, '*.tiff']);
grayscale_im = imread([directory, file_list(1).name]);


directory = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_53\Mask/';
file_list = dir([directory, '*.tiff']);
binary_mask = imread([directory, file_list(1).name]);

% Convert the binary mask to a logical array
binary_mask = logical(binary_mask);

% Create a copy of the grayscale image to modify
modified_im = grayscale_im;

% Set the grayscale value to 255 (white) where the binary mask is true
modified_im(binary_mask) = 255;

imshow(modified_im)


%% using croped zones for ploting the trajectroies
fps = 500;
pixel_size = 0.32; % micrometers

% Convert pixel size to millimeters
pixel_size_mm  = pixel_size / 1000; % mm

% Define the desired colormap based on the velocity values you want to use
velocity_range = [0 3]; % define the range of velocity values to use
cmap = colormap(jet(256)); % define a colormap (using the jet colormap as an example)
cmap = interp1(linspace(velocity_range(1),velocity_range(2),size(cmap,1)), cmap, linspace(velocity_range(1),velocity_range(2),256));

% Load the image
directory = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_53\Image_visualization\Image_stack/';
file_list = dir([directory, '*.tif']);
im = imread([directory, file_list(1).name]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% coment if no mask %%%%%%%%%%%%%%%%%%%%%%%%%%%
directory = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_53\Mask/';
file_list = dir([directory, '*.tiff']);
binary_mask = imread([directory, file_list(1).name]);

% Convert the binary mask to a logical array
binary_mask = logical(binary_mask);

% Set the grayscale value to 255 (white) where the binary mask is true
im(binary_mask) = 120;

%%%%%%%%%%%%%%%%%%%%%%%%%%% coment if no mask %%%%%%%%%%%%%%%%%%%%%%%%%%%


im_rgb = cat(3, im, im, im); % convert grayscale image to RGB format

% Plot the grayscale image and allow the user to select a rectangular box
figure(2);
imshow(im_rgb);
colormap(gray);
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


% Iterate over all particles
hold on;
for i = 1:numel(particlepositionfinale_changes)
    if size(particlepositionfinale_changes{1,i},1)>5
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
                    line(x, y, 'color', vel_color(j,:), 'linewidth', 2);
                end
            end
        end
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

%% test of funcions

figure(5)
arrow_scale = 1;
plot_velocity_arrows(im, particlepositionfinale_changes, fps, pixel_size, arrow_scale);


%% save the velocity arrows into a matrix

directory = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_53\Image_visualization\Image_stack/';
file_list = dir([directory, '*.tif']);
im = imread([directory, file_list(1).name]);

fps = 500;
pixel_size = 0.32; % micrometers

% Convert pixel size to millimeters
pixel_size_mm  = pixel_size / 1000; % mm

arrow_scale = 1;

% Initialize matrix for storing positions and velocities
pos_vel_matrix = zeros(numel(particlepositionfinale_changes), 5);

% Iterate over all particles
for i = 1:numel(particlepositionfinale_changes)
    if size(particlepositionfinale_changes{1,i},1) > 6
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
        
        nan_indices_mag = find(isnan(arrow_mag));
        nan_indices_dir = find(any(isnan(arrow_dir), 2));
        if ~isempty(nan_indices_dir)
            disp('The arrow_mag or arrow_dir matrix contains NaNs.');
        end

        % Add positions and velocities to pos_vel_matrix
        n_arrows = size(arrow_pos, 1);
        pos_vel_matrix((i-1)*n_arrows+1:i*n_arrows,1:2) = arrow_pos;
        pos_vel_matrix((i-1)*n_arrows+1:i*n_arrows,3:4) = arrow_mag.*arrow_dir;
        pos_vel_matrix((i-1)*n_arrows+1:i*n_arrows,5) = vel_mag;
    end
end


imshow(im);
hold on;
quiver(pos_vel_matrix(:,1), pos_vel_matrix(:,2), pos_vel_matrix(:,3), pos_vel_matrix(:,4), 0, 'color', 'k', 'linewidth', 1, 'MaxHeadSize', 0.8);


%% try 1

directory = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_53\Image_visualization\Image_stack/';
file_list = dir([directory, '*.tif']);
im = imread([directory, file_list(1).name]);

fps = 500;
pixel_size = 0.32; % micrometers

% Convert pixel size to millimeters
pixel_size_mm  = pixel_size / 1000; % mm

arrow_scale = 1;

% Initialize matrix for storing positions and velocities
pos_vel_matrix = zeros(numel(particlepositionfinale_changes), 5);

% Iterate over all particles
for i = 1:numel(particlepositionfinale_changes)
    if size(particlepositionfinale_changes{1,i},1) > 6
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

imshow(im);
hold on;
quiver(pos_vel_matrix(:,1), pos_vel_matrix(:,2), pos_vel_matrix(:,3), pos_vel_matrix(:,4), 0, 'color', 'k', 'linewidth', 1, 'MaxHeadSize', 0.8);



%% Plot the window arrow velocity 

% Load the image
directory = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_53\Image_visualization\Image_stack/';
file_list = dir([directory, '*.tif']);
im = imread([directory, file_list(1).name]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% coment if no mask %%%%%%%%%%%%%%%%%%%%%%%%%%%
directory = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_53\Mask/';
file_list = dir([directory, '*.tiff']);
binary_mask = imread([directory, file_list(1).name]);

% Convert the binary mask to a logical array
binary_mask = logical(binary_mask);

% Set the grayscale value to 255 (white) where the binary mask is true
im(binary_mask) = 120;


% Manually define the size of each window in pixels 
window_size = 10; % change this value to set the desired window size

% Compute the number of rows and columns in the grid
[im_height, im_width, ~] = size(im);
n_rows = ceil(im_height / window_size);
n_cols = ceil(im_width / window_size);

% Preallocate the average velocities matrix
avg_velocities = zeros(n_rows, n_cols, 2);

% Compute the window indices for each position
window_indices = ceil(pos_vel_matrix(:, 1:2) / window_size);

% Calculate the total number of windows
total_windows = n_rows * n_cols;

% Calculate the average velocities for all windows at once
for k = 1:total_windows
    k
    [r, c] = ind2sub([n_rows, n_cols], k);
    in_window = window_indices(:, 1) == c & window_indices(:, 2) == r;
    
    if any(in_window)
        avg_velocities(r, c, :) = mean(pos_vel_matrix(in_window, 3:4));
    end
end

% Create a meshgrid for the centers of the windows
[X, Y] = meshgrid((1:n_cols) * window_size - window_size/2, (1:n_rows) * window_size - window_size/2);


%% plot the images
% Plot the image and the arrows representing the average velocity in each window
figure;
imshow(im);
hold on;

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

% Add a scale factor for the arrows
scale_factor = 10; % Change this value to make the arrows bigger or smaller

% Plot the green arrows representing the average velocity in each window
green_avg_velocities = avg_velocities;
green_avg_velocities(repmat(top_10_logical, [1 1 2])) = 0;
quiver(X, Y, scale_factor * green_avg_velocities(:, :, 1), scale_factor * green_avg_velocities(:, :, 2), 'color', 'g', 'linewidth', 2, 'MaxHeadSize', 5, 'AutoScale', 'off');
hold on;

% Plot the top 10 maximum average velocity arrows in red
red_avg_velocities = zeros(size(avg_velocities));
red_avg_velocities(repmat(top_10_logical, [1 1 2])) = avg_velocities(repmat(top_10_logical, [1 1 2]));
quiver(X, Y, scale_factor * red_avg_velocities(:, :, 1), scale_factor * red_avg_velocities(:, :, 2), 'color', 'r', 'linewidth', 2, 'MaxHeadSize', 5, 'AutoScale', 'off');

% Calculate a representative arrow length based on the representative velocity value
arrow_length = scale_factor * 1;

% Add a representative green arrow to the same figure as the quiver plot
hold on;
quiver(gca, 0.8*im_width, 0.8*im_height, arrow_length, 0, 'Color', 'g', 'LineWidth', 2, 'MaxHeadSize', 5, 'AutoScale', 'off');
axis(gca, 'equal');

% Add a label with the respective velocity
text(gca, 0.8*im_width + arrow_length + 10, 0.8*im_height, sprintf('Velocity: %.2f mm/s', 1), 'FontSize', 12);



%% smoothing with gaussion
% Assume frame rate is 10 frames per second and corrected pixel size is 0.2 micrometers
fps = 1069;
pixel_size = 0.65; % micrometers

% Convert pixel size to millimeters
pixel_size_mm  = pixel_size / 1000; % mm

% Define the desired colormap based on the velocity values you want to use
velocity_range = [0 8]; % define the range of velocity values to use
cmap = colormap(jet(256)); % define a colormap (using the jet colormap as an example)
cmap = interp1(linspace(velocity_range(1),velocity_range(2),size(cmap,1)), cmap, linspace(velocity_range(1),velocity_range(2),256));

% Load the image
directory = 'H:\Cornel_ First Experiments\24.03.2023\vid_2023-03-24_15-44-29_56/';
file_list = dir([directory, '*.tiff']);
im = imread([directory, file_list(1).name]);
im_rgb = cat(3, im, im, im); % convert grayscale image to RGB format

% Plot the grayscale image
figure;
imshow(im_rgb);
colormap(gray);

% Define smoothing window size in frames
smooth_window_frames = 10;

% Iterate over all particles
hold on;
% Plot the trajectories with colors based on velocity magnitude
for i = 1:numel(particlepositionfinale_changes)
    if size(particlepositionfinale_changes{1,i},1)>20
        % Extract x-y coordinates for current particle
        pos = particlepositionfinale_changes{i}(:,2:3);
        
        % Compute time interval between consecutive frames
        dt = 1/fps; % seconds
        t = (0:size(pos,1)-2)*dt;
        
        % Compute velocity and velocity magnitude
        vel = diff(pos)./dt;
        vel_mag = sqrt(sum((vel .* repmat([pixel_size_mm pixel_size_mm], size(vel,1),1)).^2,2)); % mm/s
        
        % Remove any rows containing NaN values
        pos(any(isnan(pos), 2), :) = [];
        
        % Smooth the x and y positions separately using a Gaussian filter
        sigma = 3; % standard deviation of Gaussian filter
        smooth_window_frames = 25; % window size for Gaussian filter (in number of frames)
        xi_smooth = imgaussfilt(pos(:,1), sigma, 'Padding', 'replicate', 'FilterSize', smooth_window_frames);
        yi_smooth = imgaussfilt(pos(:,2), sigma, 'Padding', 'replicate', 'FilterSize', smooth_window_frames);
        
        % Combine the smoothed x and y positions
        pos_smooth = [xi_smooth(:) yi_smooth(:)];
        
        % Compute velocity and velocity magnitude for smoothed points
        vel_smooth = diff(pos_smooth)./dt;
        vel_mag_smooth = sqrt(sum((vel_smooth .* repmat([pixel_size_mm pixel_size_mm], size(vel_smooth,1),1)).^2,2)); % mm/s
        
        % Remove any NaN values from vel_mag_smooth and clip to the specified velocity range
        vel_mag_smooth(isnan(vel_mag_smooth)) = 0;
        vel_mag_smooth(vel_mag_smooth < velocity_range(1)) = velocity_range(1);
        vel_mag_smooth(vel_mag_smooth > velocity_range(2)) = velocity_range(2);
        
        % Define the color values for the trajectories based on velocity magnitude
        vel_color = interp1(linspace(velocity_range(1), velocity_range(2), size(cmap,1)), cmap, vel_mag_smooth);
        
        % Plot the trajectories with colors based on velocity magnitude
        for j = 1:size(pos_smooth,1)-1
            x = [pos_smooth(j,1) pos_smooth(j+1,1)];
            y = [pos_smooth(j,2) pos_smooth(j+1,2)];
            line(x, y, 'color', vel_color(j,:), 'linewidth', 2);
        end
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

%% velocity not smooth and with interpolated points.
% with no smothing and set velocity for the range we want
% Assume frame rate is 10 frames per second and corrected pixel size is 0.2 micrometers
fps = 1069;
pixel_size = 0.65; % micrometers

% Convert pixel size to millimeters
pixel_size_mm  = pixel_size / 1000; % mm

% Define the desired colormap based on the velocity values you want to use
velocity_range = [0 8]; % define the range of velocity values to use
cmap = colormap(jet(256)); % define a colormap (using the jet colormap as an example)
cmap = interp1(linspace(velocity_range(1),velocity_range(2),size(cmap,1)), cmap, linspace(velocity_range(1),velocity_range(2),256));

% Load the image
directory = 'H:\Cornel_ First Experiments\data_analysis\01.03.2023\number_11\image_contrast/';
file_list = dir([directory, '*.tif']);
im = imread([directory, file_list(1).name]);
im_rgb = cat(3, im, im, im); % convert grayscale image to RGB format

% Plot the grayscale image
figure;
imshow(im_rgb);
colormap(gray);

% Set the number of interpolated points between each trajectory point
num_interp_points = 10;

% Iterate over all particles
hold on;
for i = 2763%1:numel(particlepositionfinale_changes)
    if size(particlepositionfinale_changes{1,i},1)>20
        % Extract x-y coordinates for current particle
        pos = particlepositionfinale_changes{i}(:,2:3);

        % Interpolate the positions with 2 points between each trajectory point
        t_original = 1:size(pos,1);
        t_interp = linspace(1, size(pos,1), (size(pos,1)-1)*num_interp_points + 1);
        pos_interp = interp1(t_original, pos, t_interp, 'spline');

        % Compute time interval between consecutive frames (including the interpolated points)
        dt_interp = 1/fps/num_interp_points; % seconds

        % Compute velocity
        vel = diff(pos_interp)./dt_interp;
        vel_mag = sqrt(sum((vel .* repmat([pixel_size_mm pixel_size_mm], size(vel,1),1)).^2,2)); % mm/s

        % Remove any NaN values from vel_mag and clip to the specified velocity range
        vel_mag(isnan(vel_mag)) = 0;
        vel_mag(vel_mag < velocity_range(1)) = velocity_range(1);
        vel_mag(vel_mag > velocity_range(2)) = velocity_range(2);

        % Define the color values for the trajectories based on velocity magnitude
        vel_color = interp1(linspace(velocity_range(1),velocity_range(2),size(cmap,1)), cmap, vel_mag);

        % Plot the trajectories with colors based on velocity magnitude
        for j = 1:size(pos_interp,1)-1
            x = [pos_interp(j,1) pos_interp(j+1,1)];
            y = [pos_interp(j,2) pos_interp(j+1,2)];
            line(x, y, 'color', vel_color(j,:), 'linewidth', 2);
        end

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



%%
image_path_1 = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_63\PTV\search_radiuos_20_ALLframes\';
image_path_2 = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_63\Images_visualization\Image_stack\';
image_path_3 = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_63\Mask\';

fps = 2800;
pixel_size = 0.64 %microns
arrow_dilation = 1;
window_size = 20;
arrow_scale = 1;

plot_ptv_arrow(image_path_1,image_path_2,image_path_3 , fps, pixel_size, arrow_dilation, window_size, arrow_scale);

