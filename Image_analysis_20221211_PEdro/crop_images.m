% select the input image file
[filename, pathname] = uigetfile('*.tif', 'Select the input image file');

% read the TIFF file
tiff = imread(fullfile(pathname, filename), 'Index', 1);

% determine the number of frames in the file
num_frames = size(imfinfo(fullfile(pathname, filename)), 1);

% select the output folder
output_folder = uigetdir('', 'Select the output folder');

% crop every 10th frame
for i = 1:5:num_frames
    % read the ith frame
    frame = imread(fullfile(pathname, filename), 'Index', i);
    % save the cropped frame to the output folder
    imwrite(frame, fullfile(output_folder, sprintf('cropped_frame_%d.tif', i)));
end