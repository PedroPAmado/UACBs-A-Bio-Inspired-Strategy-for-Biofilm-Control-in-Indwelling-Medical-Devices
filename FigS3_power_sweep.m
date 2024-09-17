% Initialize global statistics
all_average_velocities = [];
all_overall_averages = [];
all_overall_stds = [];

% Load your first image to select ROI
%some_frame = imread('H:\Projects\Cilia\Cornel_ First Experiments\data_analysis\24.03.2023\power_swep\Number_61\Image_stack\stack.tif');
some_frame= imread('H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\13.04.2023\300um\300_swep5\Image_stack\stack.tif');

% Draw and get ROI position
figure; imshow(some_frame);
title('Draw a rectangle to define the ROI');
h = imrect;
roi_position = wait(h);
roi_x_min = roi_position(1);
roi_y_min = roi_position(2);
roi_x_max = roi_x_min + roi_position(3);
roi_y_max = roi_y_min + roi_position(4);

% Loop over each data set for sweep 300um 5 points
path1 = 'H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\13.04.2023\300um\300_swep1\PTV\PTV_Images_Contrast_LUT_bacground_removal_threshold.mat'
path2 = 'H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\13.04.2023\300um\300_swep2\PTV\PTV_Images_Contrast_LUT_bacground_removal_threshold.mat'
path3 = 'H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\13.04.2023\300um\300_swep3\PTV\PTV_Images_Contrast_LUT_bacground_removal_threshold.mat'
path4 = 'H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\13.04.2023\300um\300_swep4\PTV\PTV_Images_Contrast_LUT_bacground_removal_threshold.mat'
path5 = 'H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\13.04.2023\300um\300_swep5\PTV\PTV_Images_Contrast_LUT_bacground_removal_threshold.mat'

% path1 = 'H:\Projects\Cilia\Cornel_ First Experiments\data_analysis\24.03.2023\power_swep\Number_61\PTV\PTV_Images_Contrast_LUT_bacground_removal_threshold.mat'
% path2 = 'H:\Projects\Cilia\Cornel_ First Experiments\data_analysis\24.03.2023\power_swep\Number_62\PTV\PTV_Images_Contrast_LUT_bacground_removal_threshold.mat'
% path3 = 'H:\Projects\Cilia\Cornel_ First Experiments\data_analysis\24.03.2023\power_swep\Number_63\PTV\PTV_Images_Contrast_LUT_bacground_removal_threshold.mat'
% path4 = 'H:\Projects\Cilia\Cornel_ First Experiments\data_analysis\24.03.2023\power_swep\Number_64\PTV\PTV_Images_Contrast_LUT_bacground_removal_threshold.mat'


data_file_paths = {path1, path2, path3, path4, path5};
fps_values = [200, 200, 500, 500, 1000]; % Replace these with your actual fps values
pixel_to_mm = 0.30e-3; % Convert pixel to mm (0.65 micrometers = 0.65e-3 mm)

for i = 1:length(data_file_paths)
    % Load your data and image for each set
    data = load(data_file_paths{i});  % Replace this with actual loading method
    data = data.particleInfo;
    
    fps = fps_values(i);
    
    % Compute statistics for this data set
    [average_velocities_per_frame, overall_average_velocity, overall_velocity_std] = compute_velocity(data, fps, roi_x_min, roi_y_min, roi_x_max, roi_y_max,pixel_to_mm);
    
    % Aggregate statistics
    all_average_velocities = [all_average_velocities; average_velocities_per_frame];
    all_overall_averages = [all_overall_averages; overall_average_velocity];
    all_overall_stds = [all_overall_stds; overall_velocity_std];
end

%% plot

% Your previously computed all_overall_averages and all_overall_stds go here
% all_overall_averages = [your computed values];
% all_overall_stds = [your computed values];

% Define your VPP data points (replace these with your actual data)
%all_vpps = [1.2, 2, 3, 3.5];  % Or whatever VPPs correspond to your datasets

all_vpps = [22.5, 30, 37.5, 45, 52.5]; % Replace with your VPP values or similar

% Plot the data along with error bars (Standard Deviations)
figure;
errorbar(all_vpps, all_overall_averages, all_overall_stds, 'o', 'MarkerFaceColor', 'b');
xlabel('VPP');
ylabel('Average Velocity');
title('Scatter Plot of Data Points with Standard Deviations');
grid on;
hold on;

% Define the function to fit (power function: f(x) = ax^b)
fitType = fittype('a*x^b', 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b'});

% Set initial guess for parameters (you can adjust these)
initialGuess = [0.001, 2];

% Fit the data
fitResult = fit(all_vpps(:), all_overall_averages(:), fitType, 'StartPoint', initialGuess);

% Get the fitted coefficients
a = fitResult.a;
b = fitResult.b;

% Create the fitted equation string
fitEquation = ['Fitted Equation: f(x) = ', num2str(a), 'x^', num2str(b)];

% Evaluate the fitted function over a range
fittedVPP = linspace(min(all_vpps), max(all_vpps), 100);
fittedVelocity = a * fittedVPP.^b;

% Plot the best-fit curve
plot(fittedVPP, fittedVelocity, 'r-', 'LineWidth', 2);
legend('Velocity with Standard Deviations', 'Best Fit');
hold off;

% Display the fitted equation
disp(fitEquation);

%%
