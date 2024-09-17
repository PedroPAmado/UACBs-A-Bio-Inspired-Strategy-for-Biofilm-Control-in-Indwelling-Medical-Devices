function h = quiver_with_scale(X,Y,U,V,scale,units)
%QUIVER_WITH_SCALE Create a quiver plot with a scale arrow
%   h = QUIVER_WITH_SCALE(X,Y,U,V,scale,units) creates a quiver plot of
%   velocity vectors with the specified X and Y coordinates and U and V
%   components. The scale argument specifies the length of the scale arrow
%   in the same units as the vectors. The units argument is a string that
%   specifies the units to display on the scale arrow. The function returns
%   the handle to the quiver plot.

% Plot the quiver plot
h = quiver(X,Y,U,V, 'color', 'g', 'linewidth', 2, 'MaxHeadSize', 5);

% Compute the length of the vectors
L = sqrt(U.^2 + V.^2);

% Compute the size of the scale arrow
scale_arrow_length = scale / max(L(:));

% Compute the position of the scale arrow
scale_arrow_position = [0.1, 0.1];

% Compute the direction of the scale arrow
scale_arrow_direction = [1, 0];

% Plot the scale arrow
annotation('arrow', scale_arrow_position, scale_arrow_position + scale_arrow_direction * scale_arrow_length, ...
           'LineWidth', 1.5, 'HeadWidth', 8, 'HeadLength', 10, 'Color', 'k', 'LineStyle', '-');

% Add the scale label
if nargin > 5
    annotation('textbox', [scale_arrow_position + scale_arrow_direction * scale_arrow_length + [0.01, -0.005], 0.1, 0.1], ...
               'String', sprintf('%g %s', scale, units), 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end