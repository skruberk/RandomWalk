% a random walk in 2D with reflective boundaries and multiple targets to
%simulate a 2D random walk of the filament tip on surface 
% the walk is effectively diffusing point tied by an effective elastic spring to the center 
%elastic spring part not yet modeled 
 
clc;    % Clear the command window.
clearvars;
close all;  % Close all figs
workspace;  % Make sure the workspace panel is showing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters
radius = 20;    % filament movement radius variable
num_discs = 10;  % Number of discs or binding points 
target = radius / 50;   % Membrane disc variable size

% create a grid for cell
[X, Y] = meshgrid(linspace(-radius, radius, 1000), linspace(-radius, radius, 1000));
mask = (X.^2 + Y.^2 <= radius^2);

% plot the cell
pos_outer = [-radius, -radius, 2*radius, 2*radius]; 
rectangle('Position', pos_outer, 'Curvature', [1 1], 'LineWidth', 1);
axis square;
hold on;

%disc origin storage
org_coord = zeros(num_discs, 2);

% Generate and plot discs
for i = 1:num_discs
    % warnings for constant Z data
    warning('off', 'MATLAB:contour:ConstantData');

    % calculate grid row and column
    row = floor((i - 1) / sqrt(num_discs)) + 1;
    col = mod(i - 1, sqrt(num_discs)) + 1;
    
    % calculate center coordinates of the current disc
    center_x = -radius + col * radius * 2 / (sqrt(num_discs) + 1);
    center_y = -radius + row * radius * 2 / (sqrt(num_discs) + 1);

    % store the origin coordinates
    org_coord(i, :) = [center_x, center_y];

    % calculate coordinates for red circles
    theta = linspace(0, 2*pi, 100); % Angles from 0 to 2*pi
    circle_x = org_coord(1) + target * cos(theta);
    circle_y = org_coord(2) + target * sin(theta);
    
    % store the coordinates 
    circle_coordinates{i} = [circle_x; circle_y];
    
    % plot the disc clipped by the larger circle mask
    disc_masked = mask & ((X - center_x).^2 + (Y - center_y).^2 <= target^2);
    contour(X, Y, disc_masked, [0.5 0.5], 'r', 'LineWidth', 1);

    % plot masked area
    contour(X, Y, disc_masked, [0.5 0.5], 'r', 'LineWidth', 1);
    
    % restore warning state
    warning('on', 'MATLAB:contour:ConstantData');
end

% Set plot properties
ax = gca;
ax.FontWeight = 'normal';
ax.FontSize = 8;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold on;

% Open a file for writing
fileID = fopen('hitcounterwithsteps.csv', 'w');

% Write header to the CSV file
fprintf(fileID, 'Target,Hit Counter,Steps to Hit\n');


%Random Walk Section
% specify the number of time steps in our random walk
tsteps=1000;
% here is the time increment between steps
% we will assume that this is in units of seconds
tau=0.01;
% here is the diffusion coefficient for the walk
% let's assume units of um^2/sec
D=10;
% now, from the diffusion coefficient and the 
% time increment we can calculate the appropriate
% distance increment (in units of um)
delta=sqrt(2*D*tau);
num_points= 1000
% getting all points for walk
step_size = 1;  % step size for grid
theta = linspace(0, 2*pi, num_points);  % angle range for points 

counter = 0 %initialize step counter 
results = zeros(1, tsteps);  % array for counter results
hit_counter = zeros(1, num_discs); % Initialize hit counter for each target
 


% plot discs
circle_x = cell(1, num_discs); % Initialize array for circle x-coordinates for each origin 
circle_y = cell(1, num_discs); %

for i = 1:num_discs
    % Calculate coordinates for each circle
    circle_x{i} = org_coord(i, 1) + target * cos(theta);
    circle_y{i} = org_coord(i, 2) + target * sin(theta);
end

%initialize arrays to hold the x and y positions for walk
x=zeros(tsteps,1);
y=zeros(tsteps,1);
x(1) = 0.0;
y(1) = 0.0;

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
                hit_counter(k) = hit_counter(k) + 1; % Increment hit counter for the target
                fprintf('Random walk has hit target %d!\n', k);
                fprintf('Hit counter: %d\n', hit_counter(k));
                % Calculate the number of steps to hit the target
                steps_to_hit = i - sum(hit_counter(1:k-1));
                % Write data to the CSV file
                fprintf(fileID, '%d,%d,%d\n', k, hit_counter(k), steps_to_hit);
            %continue;
            elseif  norm([xbound,ybound]) >= radius %is it beyond the boundary?
                x(i)=x(i-2);
                y(i)=y(i-2);
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

