% a random walk in 2D with reflective boundaries and multiple targets to
%simulate a 2D random walk of the filament tip on surface 
%to detect a target the code generates a number of point within each circle and then 
%checks euclidean distance from the walker to each point. it runs <1% slower than code
%that checks the position of the walker vs each circle perimeter

% the walk is effectively diffusing point tied by an effective elastic spring to the center 
%elastic spring part not yet modeled 
 
clc;    % Clear the command window.
clearvars;
close all;  % Close all figs
workspace;  % Make sure the workspace panel is showing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters
radius = 20;    % filament movement radius variable
num_discs = 6;  % Number of discs or binding points 
target = radius / 10;   % Membrane disc variable size

% create a grid for cell
[X, Y] = meshgrid(linspace(-radius, radius, 1000), linspace(-radius, radius, 1000));
mask = (X.^2 + Y.^2 <= radius^2);

% plot the cell
pos_outer = [-radius, -radius, 2*radius, 2*radius]; 
rectangle('Position', pos_outer, 'Curvature', [1 1], 'LineWidth', 1);
axis square;
hold on;


%%%Disc plotting section%%%%%%%%%%%%%%%%%%%%%%%%% 
% warnings for constant Z data
warning('off', 'MATLAB:contour:ConstantData');

org_coord = zeros(num_discs, 2); % stores all the origins 
circle_coord = cell(num_discs, 2);  % stores all the coordinates 
circle_x = cell(num_discs, 1);  % Initialize as cell array
circle_y = cell(num_discs, 1);  % Initialize as cell array
theta = linspace(0, 2*pi, 100); % Angles from 0 to 2*pi

for i = 1:num_discs
    % calculate grid row and column
    row = floor((i - 1) / sqrt(num_discs)) + 1;
    col = mod(i - 1, sqrt(num_discs)) + 1;
    % calculate center coordinates of the current disc
    center_x = -radius + col * radius * 2 / (sqrt(num_discs) + 1);
    center_y = -radius + row * radius * 2 / (sqrt(num_discs) + 1);
    
    % store the origin coordinates
    org_coord(i, :) = [center_x, center_y];

    % calculate coordinates for inside red circles
    circ_points = 10;  % Number of points per circle
    spacing = 2 * pi / circ_points;  % Spacing between points

    % Initialize arrays to store x and y coordinates
    circle_x = zeros(1, circ_points);
    circle_y = zeros(1, circ_points);

    % Calculate coordinates for each point on the circle
    for j = 1:circ_points
        angle = (j - 1) * spacing;  % Compute angle
        
        % Compute x and y coordinates
        circle_x(j) = center_x + target * cos(angle);
        circle_y(j) = center_y + target * sin(angle);
    end
    
    % store the full coordinates 
    circle_coord{i} = [circle_x; circle_y];

    % plot the disc clipped by the larger circle mask
    disc_masked = mask & ((X - center_x).^2 + (Y - center_y).^2 <= target^2);
    contour(X, Y, disc_masked, [0.5 0.5], 'r', 'LineWidth', 1);
end

% restore warning state
warning('on', 'MATLAB:contour:ConstantData');

% Set plot properties
ax = gca;
ax.FontWeight = 'normal';
ax.FontSize = 8;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold on;

% Open a file for writing
fileID = fopen('hitcounterwithsteps.csv', 'a');   %w for write a for append

% write headers only if empty
if fileID == -1
    fileID = fopen('hitcounterwithsteps.csv', 'w'); % Open in write mode to create the file
    fprintf(fileID, 'Target,Hit Counter,Steps to Hit\n'); % Write headers
end
fclose(fileID);

%Random Walk Section%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%setup%%
tsteps=5000;  % number of time steps in random walk
tau=0.01; %time increment between steps units of seconds 
D=10; %diffusion coefficient for the walk units of um^2/sec
delta=sqrt(2*D*tau); % distance increment (in units of um)
num_points= 1000 % getting all points for walk
step_size = 1;  % step size for grid
theta = linspace(0, 2*pi, num_points); % angle range for points 
counter = 0 %initialize step counter 
%results = zeros(1, tsteps);  % array for counter results
hit_counter = zeros(1, num_discs); % Initialize hit counter for each target
 
%Repeat random walk%%%%% 

%results = zeros(1, num_runs);  % array for counter results

% Initialize arrays to hold the x and y positions
x = zeros(tsteps, 1);
y = zeros(tsteps, 1);

% open the file again for appending data
fileID = fopen('hitcounterwithsteps.csv', 'a');


x(1) = 0.0;
y(1) = 0.0;
hit_counter = 0;

circle_x = cell(num_discs, 1);  % Initialize as cell array
circle_y = cell(num_discs, 1);  % Initialize as cell array

for i = 1:num_discs
    % Calculate coordinates for each circle
    circle_x{i} = org_coord(i, 1) + target * cos(theta);
    circle_y{i} = org_coord(i, 2) + target * sin(theta);
end 
    
    %Random Walk
    for i = 2 :tsteps
	    % particle random step and position change 
        delx = 2 * (randi(2) - 1) - 1;   %random step 
	    xbound = x(i-1) + delx*delta;    % new position 
        dely = 2 * (randi(2) - 1) - 1;
        ybound= y(i-1) + dely*delta;
    
        %checks if particle is at any target
	    for k = 1:num_discs       
                if any(abs(xbound - circle_x{k}) <= step_size/4 & abs(ybound - circle_y{k}) <= step_size/4)
                    hit_counter = hit_counter + 1; % Increment hit counter for the target
                    fprintf('Random walk has hit target %d!\n', k);
                    %fprintf('Hit counter: %d\n', hit_counter(k));
                    % Calculate the number of steps to hit the target
                    x(i) = xbound; %update pos
                    y(i) = ybound;
                    continue;
                elseif  norm([xbound,ybound]) >= radius %is it beyond the boundary?
                    x(i)=x(i-2); %move back
                    y(i)=y(i-2);
                    continue;
                %continue;          
                else   %it's walking without hitting a target
                    x(i) = xbound; %update pos
                    y(i) = ybound;
                    counter = counter + 1 %update step counter 
                end
          end
    
end







plot(x, y);
xlabel('X Position');
ylabel('Y Position');
title('Random Walk');
axis equal; % set equal aspect ratio
grid on;



