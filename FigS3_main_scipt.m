%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           %                
%    Main code for cilia analysis           %
%                                           %
%          PEDRO PEREIRA AMADO              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Pre processing is done in FIJI for brigthness and contrast adjustment
% just do the auto contrast enhancement.

%%
clc
clear all
close all

%% align the images if needed DONT RUN

align_images_gui();

%% contrast_to_stack
% create a folder called Image_Contrast where you have the contrast images
% from fiji.

% select image path for pre_processing
image_path = 'H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\Data_to_PTV\24.3.2023\processed_data\Number_68\';
n_frames = 40; % number of frames for bacground removal
threshold = 25; % trheshold for binarization
contrast_to_stack(image_path,n_frames,threshold);


%% use fiji to get the TIFF file and place it on the correct folder for PTV
% import the prepare images for PTV to FIJI (Images_Contrast_LUT), save them as tiff and place
% them in the folder where the PTV videos are located to be processed.
input_path = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_63\Images_visualization\Images_Contrast_LUT\';
output_path = 'H:\Cornel_ First Experiments\data_analysis\scripts_24.05.23\PTV\Image_analysis_20221211_PEdro\videos\';
output_file = 'PTV_video_test';
saveImagesAsTIFF(input_path,output_path,output_file);


%% PTV - run the algorithm.

% in main:
% select the right fps

% in image processing: 
% use new general ROI selection
% use binarization
% use light intensity filtering
% use bacgroud removal
% select the running average frames - 50
% manual binarization - select a nice value

% in image analysis
% select PTV 
% select Visual PTV control
% calibration um/px - correct value of pixel
% uncheck set diameter automattically
% select the reprs. diameter [px]
% select the iterrogation window

% go to main and start
addpath('H:\Cornel_ First Experiments\data_analysis\scripts_24.05.23\PTV\Image_analysis_20221211_PEdro')

main();

%% save PTV and data
% move the data to the correct folder
% load the PTV

Numberparticle=max(particleInfo(:,6));

for ii=1:Numberparticle
    b=find(ii==particleInfo(:,6));
    
    for iii=1:length(b)
    
       particleposition_changes(iii,:)=particleInfo(b(iii),:);
 
    end 
     particlepositionfinale_changes{ii}=particleposition_changes;
 
    clear b particleposition_changes
    
end

% save the particleposition_changes with the name: positions_final

%% PTV visualization

fps = 50;
pixel_size = 0.30; % micrometers
velocity_range = [0, 1]; % mm/s
image_path_2 = 'H:\Projects\Cilia\Cornel_ First Experiments\data_analysis\13.04.23\200um\200_swep5\Image_stack\';

%need to set the path in the function

visualizeParticleTrajectories(fps, pixel_size, velocity_range, particlepositionfinale_changes,image_path_2);

%% from_ptv_to_arrows

image_path_1 = 'H:\Projects\Cilia\Cornel_ First Experiments\data_analysis\24.03.2023\power_swep\Number_64\\PTV\';
image_path_2 = 'H:\Projects\Cilia\Cornel_ First Experiments\data_analysis\24.03.2023\power_swep\Number_64\Image_stack\';
image_path_3 = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_44\Mask\';

fps = 2800;
pixel_size = 0.65 %microns
arrow_dilation = 1;
window_size = 10;
arrow_scale = 1;

[pos_vel_matrix, avg_velocities, std_velocities,std_magnitude, X, Y] = from_ptv_to_arrows(image_path_1,image_path_2,image_path_3 , fps, pixel_size, arrow_dilation, window_size, arrow_scale);




%% plot the data
image_path_1 = 'H:\Projects\Cilia\Cornel_ First Experiments\data_analysis\24.03.2023\power_swep\Number_64\\PTV\';
image_path_2 = 'H:\Projects\Cilia\Cornel_ First Experiments\data_analysis\24.03.2023\power_swep\Number_64\Image_stack\';
image_path_3 = 'H:\Cornel_ First Experiments\data_analysis\24.03.2023\Number_44\Mask\';

fps =2800;
pixel_size = 0.65 %microns
arrow_dilation = 1;
window_size = 40;
arrow_scale = 10;
velocity_range = [0 15];

plot_ptv_arrow(image_path_1,image_path_2,image_path_3 , fps, pixel_size, arrow_dilation, window_size, arrow_scale,velocity_range);





%% CAVITY
% Load the image for drawing the polygon
polygon_image = imread('H:\Projects\Cilia\Cornel_ First Experiments\data_analysis\24.03.2023\power_swep\Number_64\Image_stack\stack.tif');  % Replace 'your_image.png' with your image's file name

% Display the polygon image
figure;
imshow(polygon_image);
hold on;

% Prompt the user to draw a polygon
h = impoly;

% Wait for the user to finish drawing the polygon
wait(h);

% Get the coordinates of the polygon vertices
polygon_coords = getPosition(h);

% Close the figure
close(gcf);

% Initialize arrays for positions and velocity vectors
positions = [];
velocity_vectors = [];

% Iterate over avg_velocities
for k = 1:size(avg_velocities, 1)
    for p = 1:size(avg_velocities, 2)
        % Check if the position falls inside the polygon
        if inpolygon(X(k,p), Y(k,p), polygon_coords(:, 1), polygon_coords(:, 2))
            % Store the position and velocity vector
            positions = [positions; X(k, p), Y(k, p)];
            velocity_vector = reshape(avg_velocities(k, p, :), 1, []);
            velocity_vectors = [velocity_vectors; velocity_vector];
        end
    end
end
% Interpolate velocities
interp_velocity_vectors = velocity_vectors;
for dim = 1:2
    % Extract valid (non-NaN) velocities and their positions
    valid_velocities = velocity_vectors(~isnan(velocity_vectors(:,dim)),dim);
    valid_positions = positions(~isnan(velocity_vectors(:,dim)), :);

    % Extract invalid (NaN) velocities and their positions
    nan_positions = positions(isnan(velocity_vectors(:,dim)), :);

    if ~isempty(nan_positions)  % perform interpolation only when there are NaN velocities
        % Interpolate the velocities using 'cubic' method
        interp_velocities = griddata(valid_positions(:, 1), valid_positions(:, 2), valid_velocities, nan_positions(:, 1), nan_positions(:, 2), 'cubic');

        % Replace the NaN velocities with the interpolated ones
        interp_velocity_vectors(isnan(velocity_vectors(:,dim)), dim) = interp_velocities;
    end

    % Find any remaining NaN values
    remaining_nan = isnan(interp_velocity_vectors(:,dim));
    if any(remaining_nan)
        % Update nan_positions after 'cubic' interpolation
        nan_positions = positions(remaining_nan, :);

        % Create a scattered interpolant for the 'nearest' method
        F = scatteredInterpolant(valid_positions(:, 1), valid_positions(:, 2), valid_velocities, 'nearest');
        
        % Interpolate remaining NaN values with the 'nearest' method
        interp_velocity_vectors(remaining_nan, dim) = F(nan_positions(:, 1), nan_positions(:, 2));
    end
end
% Compute the number of positions inside the polygon
position_count = size(interp_velocity_vectors, 1);

% Compute the mean velocity magnitude
mean_velocity_magnitude = mean(sqrt(sum(interp_velocity_vectors.^2, 2)));

% Display the mean velocity magnitude
fprintf('Mean Velocity Magnitude inside the Polygon: %.2f mm/s\n', mean_velocity_magnitude);

% Plot the quiver of interpolated velocities
figure;
imshow(polygon_image);
hold on;
quiver(positions(:, 1), positions(:, 2), interp_velocity_vectors(:, 1), interp_velocity_vectors(:, 2), 'r', 'AutoScale', 'on', 'AutoScaleFactor', 1.5);
hold off;


%% TIP TO TIP
% Data extraction
X_1 = pos_vel_matrix(:, 1);
Y_1 = pos_vel_matrix(:, 2);
U_1 = pos_vel_matrix(:, 3);
V_1 = pos_vel_matrix(:, 4);

% Load your image
I =imread('H:\Cornel_ First Experiments\data_analysis\03.08.23\400um_normal_1stvideo\Images_Contrast_Aligned\vid_2023-08-03_15-24-130000.tif');  % Replace 'your_image.png' with your image's file name

% Display the image
imshow(I);
hold on

% Get the coordinates of the first point
[x1, y1] = ginput(1);
xlabel('Click the second point.');

% Get the coordinates of the second point
[x2, y2] = ginput(1);

% Calculate the ROI width and height
roi_width = 60;
roi_height = y2 - y1;

% Calculate the top-left coordinates of the ROI rectangle
roi_x = min(x1, x2);
roi_y = min(y1, y2);

% Define the rectangle coordinates
position = [roi_x, roi_y, roi_width, roi_height];

% Check if the points are inside the rectangle
inside = (X_1 >= roi_x) & (X_1 <= (roi_x + roi_width)) & ...
         (Y_1 >= roi_y) & (Y_1 <= (roi_y + roi_height));

% Triangulation and interpolation
F_u = scatteredInterpolant(X_1(inside), Y_1(inside), U_1(inside), 'natural');
F_v = scatteredInterpolant(X_1(inside), Y_1(inside), V_1(inside), 'natural');

% Define the grid points for the ROI
density = 20;
density2 = 80;
x_grid = linspace(roi_x, roi_x + roi_width, density);
y_grid = linspace(roi_y, roi_y + roi_height, density2);
[Xq, Yq] = meshgrid(x_grid, y_grid);

% Interpolate the velocities at the grid points
Uq = F_u(Xq, Yq);
Vq = F_v(Xq, Yq);

% Keep only the velocities inside the rectangle
inside_grid = (Xq >= roi_x) & (Xq <= (roi_x + roi_width)) & ...
              (Yq >= roi_y) & (Yq <= (roi_y + roi_height));
Uq(~inside_grid) = NaN;
Vq(~inside_grid) = NaN;

% Compute the mean U, V, and velocity magnitude inside the ROI
mean_u_tip_2_tip = nanmean(Uq(:));
mean_v_tip_2_tip = nanmean(Vq(:));
velocity_magnitude = sqrt(Uq.^2 + Vq.^2);
velocity_magnitude_tip_2_tip = nanmean(velocity_magnitude(:));


% Display the velocities and the rectangle
figure
imshow(I)
hold on
quiver(Xq, Yq, Uq, Vq, 'r', 'AutoScale', 'on', 'AutoScaleFactor', 1.5)
rectangle('Position', position, 'EdgeColor', 'b', 'LineWidth', 2)
hold off
axis equal
title('Velocity Field')

% Display the mean velocity magnitude
fprintf('Mean Velocity Magnitude inside the Polygon: %.2f mm/s\n', velocity_magnitude_tip_2_tip);

%% plot the points.

% Define the data points
voltage = [3.0, 3.0, 3.0];
cilia_length = [25, 75 , 125];
%velocity_tip_to_tip = [0.32, 1.042062237, 0.49];
%velocity_cavity = [0.35, 0.3, 0.11];

velocity_tip_to_tip = [0.60, 2.30, 0.42];
velocity_cavity = [0.15, 1.03, 0.57];


% Create the plot
figure;
plot(cilia_length, velocity_tip_to_tip, 'o-', 'LineWidth', 2);
hold on;
plot(cilia_length, velocity_cavity, 's-', 'LineWidth', 2);

% Set plot properties
xlabel('Cilia Distance (\mum)');
ylabel('Velocity (mm/s)');
title('Velocity vs. Cilia Distance');
legend('Velocity tip to tip', 'Velocity cavity');
grid on;

% Set axis limits
xlim([min(cilia_length) - 50, max(cilia_length) + 50]);
ylim([0, max([velocity_tip_to_tip, velocity_cavity]) + 0.1]);

% Add data labels
for i = 1:numel(voltage)
    text(cilia_length(i), velocity_tip_to_tip(i), sprintf('%.2f', velocity_tip_to_tip(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    text(cilia_length(i), velocity_cavity(i), sprintf('%.2f', velocity_cavity(i)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
end


%% fit for velocity sweep

% Define your data points
VPP = [1.2, 2, 3, 3.5];
velocity = [0.59, 3.02, 8.03, 10.23];


% Plot the data
figure;
scatter(VPP, velocity, 'o', 'filled');
xlabel('VPP');
ylabel('Velocity');
title('Scatter Plot of Data Points');
grid on;

% Define the function to fit (power function: f(x) = ax^b)
fitType = fittype('a*x^b', 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b'});

% Set initial guess for parameters (you can adjust these)
initialGuess = [1, 2];

% Fit the data
fitResult = fit(VPP', velocity', fitType, 'StartPoint', initialGuess);

% Get the fitted coefficients
a = fitResult.a;
b = fitResult.b;

% Create the fitted equation string
fitEquation = ['Fitted Equation: f(x) = ', num2str(a), 'x^', num2str(b)];

% Evaluate the fitted function over a range
fittedVPP = linspace(min(VPP), max(VPP), 100);
fittedVelocity = a * fittedVPP.^b;

% Plot the best-fit curve
hold on;
plot(fittedVPP, fittedVelocity, 'r-', 'LineWidth', 2);
legend('Data Points', 'Best Fit');
hold off;

% Display the fitted equation
disp(fitEquation);