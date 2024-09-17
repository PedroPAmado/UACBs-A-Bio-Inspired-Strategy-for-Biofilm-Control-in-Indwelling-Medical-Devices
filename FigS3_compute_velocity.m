function [average_velocities_per_frame, overall_average_velocity, overall_velocity_std] = compute_velocity(data, fps, roi_x_min, roi_y_min, roi_x_max, roi_y_max, pixel_to_mm)
    
    time_to_s = 1 / fps; % Convert frame to second
    
    % Initialize variables
    average_velocities_per_frame = [];
    
    % Loop through each frame
    unique_frames = unique(data(:, 1));
    for frame_number = unique_frames'
        % Filter particles within the current frame and ROI
        particles_in_frame = data(data(:, 1) == frame_number & ...
            data(:, 2) >= roi_x_min & data(:, 2) <= roi_x_max & ...
            data(:, 3) >= roi_y_min & data(:, 3) <= roi_y_max, :);
        
        % If there are particles in the frame and ROI, calculate average velocity
        if ~isempty(particles_in_frame)
            u_velocity = particles_in_frame(:, 4) * pixel_to_mm / time_to_s;
            v_velocity = particles_in_frame(:, 5) * pixel_to_mm / time_to_s;
            velocity_magnitudes = sqrt(u_velocity.^2 + v_velocity.^2);
            average_velocity = mean(velocity_magnitudes);
            average_velocities_per_frame = [average_velocities_per_frame; average_velocity];
        end
    end
    
    % Calculate overall average and standard deviation of velocity magnitudes
    overall_average_velocity = mean(average_velocities_per_frame);
    overall_velocity_std = std(average_velocities_per_frame);
end
