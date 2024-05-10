% a random walk in 2D with reflective boundaries and multiple targets to
%simulate a 2D random walk of the filament tip on surface 
% the walk is effectively diffusing point tied by an effective elastic spring to the center 
%elastic spring part not yet modeled 
%this is an older version where:
%the calculation for whether the random walk reaches the target is performed by iterating through each 
%target and checking if the walker is closeby. a newer version checks the position of the walker and 
%then iterates through the origins of each target to see if its close by and then checks the coordinates,
%which should be much more efficient. 
 
clc;    % Clear the command window.
clearvars;
close all;  % Close all figs
workspace;  % Make sure the workspace panel is showing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters here
radius = 20;    % filament movement radius variable
num_discs = 20;  % Number of discs or binding points 
target = radius / 10;   % Membrane disc variable size
tsteps=10000;  % number of time steps in random walk
%Repeat random walk%%%%% 
num_runs = 1000;%how many times to run code

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

% calculate coordinates for inside red circles
circ_points = target*(radius*2.5);  % Number of points per circle
org_coord = zeros(num_discs, 2); % stores all the origins 
circle_coord = cell(num_discs, 2);  % stores all the coordinates 
  
% Calculate coordinates for points within the circle
radius_points = target * sqrt(rand(1, circ_points)); % Random radius values
angle_points = linspace(0, 2*pi, circ_points); % Angles for equally spaced points

for i = 1:num_discs   % place all circles 
    % calculate grid row and column
    row = floor((i - 1) / sqrt(num_discs)) + 1;
    col = mod(i - 1, sqrt(num_discs)) + 1;
    % calculate center coordinates of the current disc
    center_x = -radius + col * radius * 2 / (sqrt(num_discs) + 1);
    center_y = -radius + row * radius * 2 / (sqrt(num_discs) + 1);
    
    % store the origin coordinates
    org_coord(i, :) = [center_x, center_y];

    % Initialize arrays to store x and y coordinates
    %circle_x = zeros(1, circ_points);
    %circle_y = zeros(1, circ_points);
    circle_x_point = cell(1, num_discs);
    circle_y_point = cell(1, num_discs);
    % Calculate coordinates for each point on the circle
    for j = 1:circ_points
        % Generate random radius within the target radius
        r = rand(1) * target;
        % Generate random angle
        angle = rand(1) * 2 * pi;
        % Compute x and y coordinates using polar coordinates
        circle_x_point{j} = center_x + r * cos(angle);
        circle_y_point{j} = center_y + r * sin(angle);
    end
   
    % store the full coordinates 
    circle_coord{i} = [circle_x_point; circle_y_point];

    % plot the disc clipped by the larger circle mask
    disc_masked = mask & ((X - center_x).^2 + (Y - center_y).^2 <= target^2);
    contour(X, Y, disc_masked, [0.5 0.5], 'r', 'LineWidth', 1);

    % Plot the points within the circle to check them 
    %scatter(circle_x_point, circle_y_point, 10, 'b', 'filled');
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
tau=0.01; %time increment between steps units of seconds 
D=10; %diffusion coefficient for the walk units of um^2/sec
delta=sqrt(2*D*tau); % distance increment (in units of um)
num_points= 1000 % getting all points for walk
step_size = 1;  % step size for grid
theta = linspace(0, 2*pi, num_points); % angle range for points 
counter = 0 %initialize step counter 
%results = zeros(1, tsteps);  % array for counter results
hit_counter = zeros(1, num_discs); % Initialize hit counter for each target
 

%results = zeros(1, num_runs);  % array for counter results

% Initialize arrays to hold the x and y positions
x = zeros(tsteps, 1);
y = zeros(tsteps, 1);

% open the file again for appending data
fileID = fopen('hitcounterwithsteps.csv', 'a');

for run = 1:num_runs
    %reset each run 
    x(1) = 0.0;
    y(1) = 0.0;
    hit_counter(:) = 0;

%circle_x = cell(num_discs, 1);  % Initialize as cell array
%circle_y = cell(num_discs, 1);  % Initialize as cell array

for i = 1:num_discs
    % Calculate coordinates for each circle
    circle_x{i} = org_coord(i, 1) + target * cos(theta);
    circle_y{i} = org_coord(i, 2) + target * sin(theta);
    circle_x_point{i} = zeros(1, circ_points);
    circle_y_point{i} = zeros(1, circ_points);
end 
 
k=0 %initialize loop variable

    %Random Walk
    for i = 2 :tsteps
	    % particle random step and position change 
        delx = 2 * (randi(2) - 1) - 1;   %random step 
	    xbound = x(i-1) + delx*delta;    % new position 
        dely = 2 * (randi(2) - 1) - 1;
        ybound= y(i-1) + dely*delta;
    
        %checks if particle is at any target
	    for k = 1:num_discs       
                if any(abs(xbound - circle_x_point{k}) <= step_size/4 & abs(ybound - circle_y_point{k}) <= step_size/4)
                    hit_counter(k) = hit_counter(k) + 1; % Increment hit counter for the target
                    fprintf('Random walk has hit target %d!\n', k);
                    %fprintf('Hit counter: %d\n', hit_counter(k));
                    % Calculate the number of steps to hit the target
                    steps_to_hit = i - sum(hit_counter(1:k-1));
                    % Write data to the CSV file
                    fprintf(fileID, '%d,%d,%d\n', k, hit_counter(k), steps_to_hit);
                    x(i) = xbound; %update pos
                    y(i) = ybound;
                    continue;
                elseif  norm([xbound,ybound]) >= radius %is it beyond the boundary?
                    x(i)=x(i-2);
                    y(i)=y(i-2); 
                    continue;
                else   %it's walking without hitting a target
                    x(i) = xbound; %update pos
                    y(i) = ybound;
                    counter = counter + 1 %update step counter 
                end
          end
    end
end


% Close the file
fclose(fileID);

% Read the data from the CSV file
%data = readtable('hitcounterwithsteps.csv');
opts = detectImportOptions('hitcounterwithsteps.csv');
opts.VariableNamingRule = 'preserve'; % Preserve original column headers
data = readtable('hitcounterwithsteps.csv', opts);

% Extract the steps to hit from the data and convert to array
hit_counter = table2array(data(:, 2));

% Plot the histogram
histogram(hit_counter, 'BinWidth', 1);  %'BinWidth', 1, 'Normalization', 'probability'
xlabel('Steps to Hit', 'FontSize', 14);
ylabel('count', 'FontSize', 14);
title(['Histogram of Steps to Hit (' , num2str(num_discs), 'discs and ' num2str(target), ' target radius)']);
%xlim([0, 550])
