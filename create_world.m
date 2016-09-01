function [X, Y] = create_world(length, width, grid)
% Creates the mesh grid on which wave elevation
% is plotted and defines other variables relating to the world

% Assuming fishing on continental shelf, depth is 50 - 200 m
% If deep ocean, depth is 1000 m 
% according to http://www.gov.scot/Resource/Doc/158590/0043011.pdf

% X and Y axis 
x = 0:grid:length;   % longitudinal axis
y = 0:grid:width;    % lateral axis

[X, Y] = meshgrid(x, y);
end

