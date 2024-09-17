% Load image sequence from a folder
folder = 'H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\flapping\flapping_C001H001S0001'; % Replace with the path to your image folder
%folder = 'H:\Pedro\Data_from_old_hardisk\Projects\Cilia\Processed data\flapping\flapping2_C001H001S0001'; % Replace with the path to your image folder

imageFiles = dir(fullfile(folder, '*.tif')); % Change extension to '*.tif' for TIFF images
numFrames = 50 %numel(imageFiles);
imageSequence = cell(1, numFrames);

for i = 1:numFrames
    % Load the image as uint16
    uint16Image = imread(fullfile(folder, imageFiles(i).name));
    
    % Scale the image values to the range [0, 1]
    scaledImage = mat2gray(uint16Image);
    
    % Convert the scaled image to uint8
    uint8Image = im2uint8(scaledImage);
    
    % Store the converted image in the cell array
    imageSequence{i} = uint8Image;
end

% Display the first frame and allow user to select a point using the mouse
figure;
imshow(imageSequence{1});
title('Select a point');

% Wait for user to select a point
[x, y] = ginput(1); % Get the x and y coordinates of the selected point

% Extract 8 horizontal pixels to the right of the selected point in each frame
x_start = round(x);
y_start = round(y);
x_end = x_start + 7; % Extract 8 pixels to the right

selectedPixels = zeros(7, numFrames);

for i = 1: numFrames
    frame = imageSequence{i};
    selectedPixels(:, i) = frame(y_start, x_start:x_end-1);
end

% Create a new image with the extracted pixels organized in columns
newImage = reshape(selectedPixels, [], numFrames);

% Display the new image
figure;
imshow(newImage, []);
title('Extracted and rearranged pixels');


%%
% Define the frequency and amplitude of the sine wave
frequency = 11000; % Adjust the frequency as needed
amplitude = 1.5; % Adjust the amplitude as needed

% Calculate the time points based on the frame number and fps
fps = 100000; % Change this to your actual frames per second
timePoints = (1:numFrames) / fps;

% Generate the sine wave
sineWave = amplitude * sin(2 * pi * frequency * timePoints + 5.4);

% Specify the vertical offset
verticalOffset =4.5; % Adjust the offset as needed

% Overlay the sine wave on the image
figure;
imshow(newImage, []);
hold on;
plot(sineWave + verticalOffset, 'r', 'LineWidth', 2); % Add the vertical offset
title('Image with overlaid sine wave');
hold off;

