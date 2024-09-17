% Specify the folder path containing the images
folder_path = 'H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\21.07.2023\300_old_cavity_cleaning and desin\Crystal segmentation\Segmented';

% Load the ROI image
roi_image = imread('H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\21.07.2023\300_old_cavity_cleaning and desin\Crystal segmentation\Area\roi_seconds.tiff');  % Replace with the actual path

mask_image = imread('H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\21.07.2023\300_old_cavity_cleaning and desin\Crystal segmentation\Mask\mask_seconds.tiff');  % Replace with the actual path


% Get the binary mask of the ROI
roi_binary = imbinarize(roi_image);

% Get a list of all image files in the folder
image_files = dir(fullfile(folder_path, '*.tiff'));  % Modify the file extension if needed

starting_frame = 1;

% Initialize variables for storing results
num_images = numel(image_files);
num_blobs = zeros(1, num_images-starting_frame);
blob_areas = cell(1, num_images-starting_frame);

% Process each image in the folder
for i = starting_frame:num_images
    
    save_idx = i-starting_frame+1;
    % Load the image
    image_path = fullfile(folder_path, image_files(i).name);
    image = imread(image_path);
    
    % Apply the mask to the image
    masked_image = bsxfun(@times, image, cast(~mask_image, 'like', image));
    
    % Convert the masked image to binary
    binary_image = imbinarize(masked_image);
    
    % Apply the ROI as a mask to the binary image
    binary_image_roi = binary_image & roi_binary;
    
    % Perform connected component analysis on the ROI-masked binary image
    cc_roi = bwconncomp(binary_image_roi);
    
    % Get the number of connected components
    num_blobs(save_idx) = cc_roi.NumObjects;
    
    % Get the area of each connected component
    blob_areas{save_idx} = cellfun(@numel, cc_roi.PixelIdxList);
    
    % Filter blobs based on the threshold
    threshold = 5;  % Adjust the threshold as desired
    valid_blobs = blob_areas{save_idx} > threshold;
    num_blobs(save_idx) = sum(valid_blobs);
    blob_areas{save_idx} = blob_areas{save_idx}(valid_blobs);
    
    % Display the image with filled blob regions for the first image
    if i == 1400
        figure;
        imshow(image);
        hold on;
        
        % Overlay the ROI image
        roi_overlay = imfuse(image, roi_image, 'blend', 'Scaling', 'joint');
        h_roi_overlay = imshow(roi_overlay);
        set(h_roi_overlay, 'AlphaData', 0.5);  % Adjust the transparency as needed
        
        % Display blob boundaries within the ROI only
        boundaries = bwboundaries(binary_image_roi);
        for j = 1:length(boundaries)
            boundary = boundaries{j};
            plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2);
        end
        
        hold off;
        title('Blob Contours and ROI Overlay in Image');
    end
end


%% plot

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultfigurecolor',[1 1 1])

fps = 200;  % Frames per second
starting_frame = 1;
num_images = 3119;
time = (starting_frame:num_images) / fps;
mean_areas_um2 = cellfun(@(x) max(x) * (0.7^2), blob_areas);  % Compute mean area in square micrometers
time = time - time(303);

% Define a custom color palette (you can choose your preferred colors)
lineColor = [0, 0.4470, 0.7410];   % Aesthetic blue color for lines
bgColor = [1, 1, 1];  % White background color

% Apply moving average filter
window_size = 100;  % Adjust the window size as desired
num_blobs_smoothed = movmean(num_blobs, window_size);
mean_areas_um2_smoothed = movmean(mean_areas_um2, window_size);
mean_areas_um2_smoothed = movmean(mean_areas_um2_smoothed, window_size);

% Create a figure for the smoothed plot
fig_smoothed = figure('Position', [100, 100, 800, 800], 'Color', 'white');
ax_smoothed = gca;

% Plot the smoothed data with enhanced aesthetics
plot(time, mean_areas_um2_smoothed, 'Color', lineColor, 'LineWidth', 8);

% Adjust x-axis limits
xlim([-2, 13]);

% Calculate the maximum value for y-axis and add 10% buffer
yMaxSmoothed = max(mean_areas_um2_smoothed) * 1.1; 

% Adjust the y-axis limits
ylim([0, yMaxSmoothed]);

% Set custom y-axis ticks with larger spacing
ytick_interval = 10000; % Adjust this value to increase the distance between ticks
yticks(0:ytick_interval:ceil(yMaxSmoothed));

% Set custom x-axis ticks (independently from y-axis ticks)
xtick_interval = 1;  % Adjust this value to set x-axis tick spacing
xticks(-2:xtick_interval:13);

% Remove grid lines
grid off;

% Add a background color to the plot area
ax_smoothed.Color = bgColor; % White background

% Customize the axis labels and title for the smoothed plot with increased font size
ylabel(strcat('\bf Maximum Area',' [$\bf \mu m^2$]'),'Interpreter','latex', 'FontSize', 25);
xlabel(strcat('\bf Time (s)'),'Interpreter','latex', 'FontSize', 25);
title(strcat('\bf Mean Blob Area over Time'),'Interpreter','latex', 'FontSize', 25);

% Customize the legend for the smoothed plot with increased font size
legend(strcat('\bf Biggest Blob Area'),'Interpreter','latex', 'FontSize', 16);

% Re-apply the y-axis limits in case they were adjusted after setting yticks
ylim([0, yMaxSmoothed]);
yticks(0:ytick_interval:ceil(yMaxSmoothed));

% Set the aspect ratio of the plot
pbaspect([1 1 1]);  % This sets the aspect ratio to be square (1:1:1)

% Set LaTeX interpreter and font size for axes
ax_smoothed.TickLabelInterpreter = 'latex';
ax_smoothed.XAxis.FontSize = 30;
ax_smoothed.YAxis.FontSize = 30;
