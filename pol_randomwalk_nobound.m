% a random walk in 2D with reflective boundaries and multiple targets to
%simulate a 2D random walk of the filament tip on surface 
%determines whether something has hit a target by checking distance vs
%circle perimeter
%once walker has hit target it relocates to a random starting position 
%this code runs much faster as there is no circular boundary- the walker can
%explore the whole grid 
%runs the code a number of times and generates a histogram 
% the walk is effectively a diffusing point tied by an effective elastic spring to the center 
%elastic spring part not yet modeled 

clc;    % Clear the command window.
clearvars;
close all;  % Close all figs
workspace;  % Make sure the workspace panel is showing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters here
radius = 20;    % filament movement radius variable
num_discs = 30;  % Number of discs or binding points 
target = radius / 20;   % Membrane disc variable size
target = round(target,1);
tsteps=1000;  % number of time steps in random walk
%Repeat random walk%%%%% 
num_runs = 10;%how many times to run code

% create a grid for cell
[X, Y] = meshgrid(linspace(-radius, radius, 1000), linspace(-radius, radius, 1000));
%mask = (X.^2 + Y.^2 <= radius^2);

% plot the grid
%pos_outer = [-radius, -radius, 2*radius, 2*radius]; 
%rectangle('Position', pos_outer, 'Curvature', [1 1], 'LineWidth', 1);
axis square;
hold on;

%%%Disc plotting section%%%%%%%%%%%%%%%%%%%%%%%%% 
% warnings for constant Z data
warning('off', 'MATLAB:contour:ConstantData');

% calculate coordinates for red circles
org_coord = zeros(num_discs, 2); % stores all the origins 
circle_coord = zeros(num_discs, 2);  % stores all the coordinates 
% calculate coordinates for red circles
circ_points = 50;  % Number of points per circle
spacing = 2 * pi / circ_points;  % Spacing between points
radius_points = target * sqrt(rand(1, circ_points)); % Random radius values
angle_points = linspace(0, 2*pi, circ_points); % Angles for equally spaced points

% place all circles 
for i = 1:num_discs   
    % Calculate grid row and column
    row = floor((i - 1) / sqrt(num_discs)) + 1;
    col = mod(i - 1, sqrt(num_discs)) + 1;
    % Calculate center coordinates of the current disc
    center_x = -radius + col * radius * 2 / (sqrt(num_discs) + 1);
    center_y = -radius + row * radius * 2 / (sqrt(num_discs) + 1);
    
    % Store the origin coordinates
    org_coord(i, :) = [center_x, center_y];
    
    % Plot the discs inside the loop
    disc_draw = ((X - center_x).^2 + (Y - center_y).^2 <= target^2);
    contour(X, Y, disc_draw, [0.5 0.5], 'r', 'LineWidth', 1);
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

%Random Walk Section%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%setup%%
tau=0.01; %time increment between steps units of seconds 
D=10; %diffusion coefficient for the walk units of um^2/sec
delta=sqrt(2*D*tau); % distance increment (in units of um)
num_points= 1000; % getting all points for walk
step_size = 1;  % step size for grid
theta = linspace(0, 2*pi, num_points); % angle range for points 
counter = 0; %initialize step counter 
%results = zeros(1, tsteps);  % array for counter results
hit_counter = zeros(1, num_discs); % Initialize hit counter for each target
 
% Initialize arrays to hold the x and y positions
x = zeros(tsteps, 1);
y = zeros(tsteps, 1);

circle_x = cell(num_discs, 1);  % initialize
circle_y = cell(num_discs, 1);  
for i = 1:num_discs
    % Calculate coordinates for each circle
    circle_x{i} = org_coord(i, 1) + target * cos(theta);
    circle_y{i} = org_coord(i, 2) + target * sin(theta);
end 

color_palette = lines(num_runs);

for run = 1:num_runs
    %reset each run 
    x(1) = 0.0;
    y(1) = 0.0;
    hit_counter(:) = 0;
    % Initialize logical array to mark targets as hit
    targets_hit = false(1, num_discs);

    %Random Walk
    for step = 2:tsteps
    % particle random step and position change 
    delx = 2 * (randi(2) - 1) - 1;   % random step 
    xbound = x(step - 1) + delx * delta;    % new position 
    dely = 2 * (randi(2) - 1) - 1;
    ybound = y(step - 1) + dely * delta;

         % Check if the particle is at any target
        for k = 1:num_discs       
            if any(abs(xbound - circle_x{k}) <= step_size/10 & abs(ybound - circle_y{k}) <= step_size/10)
                hit_counter(k) = hit_counter(k) + 1; % Increment hit counter for the target
                fprintf('Random walk has hit target %d!\n', k)
                % Calculate the number of steps to hit the target
                steps_to_hit = step - sum(hit_counter(1:k-1));
                % Write data to the CSV file
                %fprintf(fileID, '%d,%d,%d\n', k, hit_counter(k), steps_to_hit);
                % return walker 
                x(step) = min(X(:)) + rand() * (max(X(:)) - min(X(:)));
                y(step) = min(Y(:)) + rand() * (max(Y(:)) - min(Y(:)));
                break;
            elseif norm([xbound,ybound]) >= radius %is it beyond the boundary?
                x(step) = min(max(xbound, -radius), radius);
                %x(step) = x(step) + rand() * 2 * delta - delta;
                %x(step)=-x(step-1);    
                y(step) = min(max(ybound, -radius), radius);
                %y(step) = y(step) + rand() * 2 * delta - delta;
                %y(step)=-y(step-1);
            else   %it's walking without hitting a target
               
                x(step) = xbound; % update pos
                y(step) = ybound;
                counter = counter + 1 % update step counter 
            end
        end
    end
plot(x, y,'Color', color_palette(run,:), 'LineWidth', 1.5);
xlabel('X Position');
ylabel('Y Position');
title('Random Walk');
axis equal; % set equal aspect ratio
grid on;
hold on;
end




