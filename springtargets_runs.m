% a random walk in 2D with reflective boundaries and multiple targets to
%simulate a 2D random walk of the filament tip on surface
%walker explores a circular region on a grid
%each of the targets performs a restricted random walk and is tied to their
%origin using a variable spring constant 
%determines whether something has hit a target by checking distance vs
%circle perimeter and then sending the walker back to a random point after a
%hit
%runs the code a number of times and generates a histogram 

clc;    % Clear the command window.
clearvars;
close all;  % Close all figs
workspace;  % Make sure the workspace panel is showing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters here
radius = 20;    % filament movement radius variable
num_discs = 10;  % Number of discs or binding points 
target = radius / 20;   % Membrane disc variable size
target = round(target,1);
tsteps=1000;  % number of time steps in random walk
%Repeat random walk%%%%% 
num_runs = 10000;%how many times to run code

%%filehandling%%%%%%%%%%%%%%%%%
filename = sprintf('steps2hit_refltarget_%d_%d_%d_%d.csv', num_discs, target, tsteps,num_runs);
% Open a file for writing
fileID = fopen(filename, 'a');   %w for write a for append
% Check if the file is empty
fileEmpty = fseek(fileID, 0, 'eof') == 0;
% Write headers only if the file is empty
if fileEmpty
    fprintf(fileID, 'Target,Hit Counter,Steps to Hit\n'); % Write headers
end
fclose(fileID);

% create a grid for cell
[X, Y] = meshgrid(linspace(-radius, radius, 1000), linspace(-radius, radius, 1000));
mask = (X.^2 + Y.^2 <= radius^2);

% plot the grid
%pos_outer = [-radius, -radius, 2*radius, 2*radius]; 
%rectangle('Position', pos_outer, 'Curvature', [1 1], 'LineWidth', 1);
%axis square;
%hold on;

%%%Disc plotting section%%%%%%%%%%%%%%%%%%%%%%%%% 
% warnings for constant Z data
warning('off', 'MATLAB:contour:ConstantData');

org_coord = zeros(num_discs, 2); % stores all the origins 
circle_coord = zeros(num_discs, 2);  % stores all the coordinates 
% calculate coordinates for red circles
circ_points = 50;  % Number of points per circle
spacing = 2 * pi / circ_points;  % Spacing between points

% Calculate coordinates for points within the circle
radius_points = target * sqrt(rand(1, circ_points)); % Random radius values
angle_points = linspace(0, 2*pi, circ_points); % Angles for equally spaced points
    % calculate disc centers
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

    % plot the disc clipped by the larger circle mask
    %disc_masked = mask & ((X - center_x).^2 + (Y - center_y).^2 <= target^2);
    %contour(X, Y, disc_masked, [0.5 0.5], 'r', 'LineWidth', 1);

% restore warning state
warning('on', 'MATLAB:contour:ConstantData');

% Set plot properties
ax = gca;
ax.FontWeight = 'normal';
ax.FontSize = 8;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
%hold on;

%Random Walk Section%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%setup%%
tau=0.01; %time increment between steps units of seconds 
D=10; %diffusion coefficient for the walk units of um^2/sec
delta=sqrt(2*D*tau); % distance increment (in units of um)
num_points= 1000; % getting all points for walk
step_size = 1;  % step size for grid
theta = linspace(0, 2*pi, num_points); % angle range for points 
counter = 0; %initialize step counter 
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
% open the file again for appending data%%%%%%%%%%%%%%%%%%%%
fileID = fopen(filename, 'a');

spring_constant = 0.1;  
target_x = zeros(num_discs, 1);
target_y = zeros(num_discs, 1);
for run = 1:num_runs
    % Reset each run 
    x(1) = 0.0;
    y(1) = 0.0;
    hit_counter(:) = 0;
    % Initialize target positions
    for k = 1:num_discs
        target_x(k) = org_coord(k, 1);
        target_y(k) = org_coord(k, 2);
    end
% Random Walk
    for step = 2:tsteps
        % Particle random step and position change 
        delx = 2 * (randi(2) - 1) - 1;   % Random step 
        xbound = x(step - 1) + delx * delta;    % New position 
        dely = 2 * (randi(2) - 1) - 1;
        ybound = y(step - 1) + dely * delta;
        % add spring forces to target positions
        for k = 1:num_discs
            spring_force_x = spring_constant * (org_coord(k, 1) - target_x(k));
            spring_force_y = spring_constant * (org_coord(k, 2) - target_y(k));
            target_x(k) = target_x(k) + spring_force_x;
            target_y(k) = target_y(k) + spring_force_y;
            %limited random walk for each target
            target_x(k) = target_x(k) + randn() * 0.1;  % spring target constant 0.1-0.5 works 
            target_y(k) = target_y(k) + randn() * 0.1;
            %collision check, euclidean distance 
            if sqrt((xbound - target_x(k))^2 + (ybound - target_y(k))^2) <= target
                hit_counter(k) = hit_counter(k) + 1; % Increment hit counter for the target
                fprintf('Random walk has hit target %d!\n', k)
                % number of steps to hit the target
                steps_to_hit = step - sum(hit_counter(1:k-1));
                % Write data to the CSV file
                fprintf(fileID, '%d,%d,%d\n', k, hit_counter(k), steps_to_hit);
                % return walker 
                x(step) = rand() * (2 * radius) - radius;
                y(step) = rand() * (2 * radius) - radius;
                break;
            elseif norm([xbound,ybound]) >= radius %is it beyond the boundary?
                x(step) = min(max(xbound, -radius), radius); %reflect over boundary
                x(step)=x(step-2); %step back
                y(step) = min(max(ybound, -radius), radius);
                y(step)=y(step-2);
            else   %it's walking without hitting a target
                x(step) = xbound; % update pos
                y(step) = ybound;
                counter = counter + 1 % update step counter 
            end
        end
    end
end

% Close the file
fclose(fileID);

% Read the data from the CSV file
opts = detectImportOptions(filename);
opts.VariableNamesLine = 1; % Assuming the headers are in the first line
opts.VariableNamingRule = 'preserve'; % Preserve original column headers
data = readtable(filename, opts);
% Extract the steps to hit from the data and convert to array
steps_to_hit = table2array(data(:, 3));

% Plot the histogram
histogram(steps_to_hit, 'BinWidth', 10);  %'BinWidth', 1, 'Normalization', 'probability'
xlabel('Steps to Hit', 'FontSize', 14);
ylabel('count', 'FontSize', 14);
title(['Histogram of Steps to Hit (' , num2str(num_discs), 'discs and ' num2str(target), ' target radius)']);
%xlim([0, 10000])
% Retrieve the current y-axis tick values
%yTicks = get(gca, 'YTick');
%tick values to strings without scientific notation
%yTickLabels = arrayfun(@(x) sprintf('%d', x), yTicks, 'UniformOutput', false);
% Apply the new tick labels to the y-axis
%set(gca, 'YTickLabel', yTickLabels);

