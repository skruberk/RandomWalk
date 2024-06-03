% a 2D random walk with reflective boundaries and multiple targets to
% simulate the interaction of an actin filament tip on a surface
%walker explores a circular region on a grid
%each of the targets also performs a restricted random walk and is tied to their
%origin using a variable spring constant 
%determines whether something has hit a target by checking distance vs
%circle perimeter and then sending the walker back to a random point after a
%hit
%runs the code a number of times and generates a movie of the walk, also
%plots each of the target positions 

clc;    % Clear the command window.
clearvars;
close all;  % Close all figs
workspace;  % Make sure the workspace panel is showing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters here
radius = 20;    % filament movement radius variable
num_discs = 1;  % Number of discs or binding points for movie only works for 1 target
target = radius / 6;   % Membrane disc variable size
target = round(target,1);
tsteps=10000;  % number of time steps in random walk
%Repeat random walk%%%%% 
num_runs = 1;%how many times to run code, doesn't super work for more than one run

% create a grid for cell
[X, Y] = meshgrid(linspace(-radius, radius, 1000), linspace(-radius, radius, 1000));
mask = (X.^2 + Y.^2 <= radius^2);

% plot the grid
pos_outer = [-radius, -radius, 2*radius, 2*radius]; 
rectangle('Position', pos_outer, 'Curvature', [1 1], 'LineWidth', 1);
axis square;
hold on;

%%%Disc plotting section%%%%%%%%%%%%%%%%%%%%%%%%% 
% warnings for constant Z data
warning('off', 'MATLAB:contour:ConstantData');

org_coord = zeros(num_discs, 2); % stores all the origins 
%circle_coord = zeros(num_discs, 2);  % stores all the coordinates 
% calculate coordinates for red circles
circ_points = 50;  % Number of points per circle
spacing = 2 * pi / circ_points;  % Spacing between points

% Calculate coordinates for points on the circle
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
    
    % plot the discs inside the loop THIS CHANGED 
    disc_draw = mask & ((X - center_x).^2 + (Y - center_y).^2 <= target^2);
    contour(X, Y, disc_draw, [0.5 0.5], 'r', 'LineWidth', 1);
    %disc_draw = ((X - center_x).^2 + (Y - center_y).^2 <= target^2);
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
tau=0.1; %time increment between steps units of seconds 
D=10; %diffusion coefficient for the walk units of um^2/sec
delta=sqrt(2*D*tau); % distance increment (in units of um)
num_points= 10000; % getting all points for walk
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

color_palette = lines(num_runs);
spring_constant = 0.01;  
target_x = zeros(num_discs, 1);
target_y = zeros(num_discs, 1);

for run = 1:num_runs
    figure; % Open a new figure for each run
    hold on;
    axis equal;
    axis([-radius radius -radius radius]);
    rectangle('Position', [-radius, -radius, 2*radius, 2*radius], 'Curvature', [1 1], 'LineWidth', 1);
    h = plot(x, y, 'Color', color_palette(run,:), 'LineWidth', 1.5); % Plot walker track
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
        %delx = 2 * (randi(2) - 1) - 1;   % Random step 
        %xbound = x(step - 1) + delx * delta;    % New position 
        %dely = 2 * (randi(2) - 1) - 1;
        %ybound = y(step - 1) + dely * delta;
        %other way to calculate the walk 
        angle = rand() * 2 * pi;  % Random angle between 0 and 2*pi
        delx = delta * cos(angle);  % Step size in the x direction
        dely = delta * sin(angle);  % Step size in the y direction
        xbound = x(step - 1) + delx;  % New position in x direction
        ybound = y(step - 1) + dely;  % New position in y direction
        % Apply spring forces to target positions
        for k = 1:num_discs
            spring_force_x = spring_constant * (org_coord(k, 1) - target_x(k));
            spring_force_y = spring_constant * (org_coord(k, 2) - target_y(k));
            target_x(k) = target_x(k) + spring_force_x;
            target_y(k) = target_y(k) + spring_force_y;
            % random walk for each target
            target_x(k) = target_x(k) + randn() * 0.5;  % target walk step size 0.1-0.5 works well
            target_y(k) = target_y(k) + randn() * 0.5;
            % make targets stay within the larger circle  
            if sqrt(target_x(k)^2 + target_y(k)^2) > radius
            % Adjust position to be within the circle boundary
                angle = atan2(target_y(k), target_x(k));
                distance = sqrt(target_x(k)^2 + target_y(k)^2);
                scaling_factor = radius / distance * rand() * 0.9; % random btw 0 - 0.9
                target_x(k) = target_x(k) * scaling_factor;
                target_y(k) = target_y(k) * scaling_factor;
            end
            %collision check 
            if sqrt((xbound - target_x(k))^2 + (ybound - target_y(k))^2) <= target
                fprintf('Random walk has hit target %d!\n', k)
                % Calculate the number of steps to hit the target
                % return walker 
                x(step) = rand() * (2 * radius) - radius;
                y(step) = rand() * (2 * radius) - radius;
                break;
            elseif norm([xbound,ybound]) >= radius %is it beyond the boundary?     
                % If the walker is outside and xbound is positive and
                % ybound is positive etc
                % wrap walker around to the other side
                    if xbound > 0 && ybound > 0
                        x(step) = (-abs(xbound))+2; % wrap to left
                        y(step) = (-abs(ybound))+2; % wrap to bottom  
                    elseif xbound < 0 && ybound > 0
                        x(step) = (abs(xbound))-2; % wrap to the right 
                        y(step) = (-abs(ybound))+2;% wrap to bottom 
                    elseif xbound > 0 && ybound < 0
                        x(step) = (-abs(xbound))+2; % wrap to left
                        y(step) = (abs(ybound))-2; % wrap to top 
                    elseif xbound < 0 && ybound < 0
                        x(step) =  (abs(xbound))-2; % wrap to the right 
                        y(step) =  (abs(ybound))-2; % wrap to top 
                    end
            else   %it's walking without hitting a target
                x(step) = xbound; % update pos
                y(step) = ybound;    
            end
            % walker's track
            set(h, 'XData', x(1:step), 'YData', y(1:step));
            % plot updated targets
            for m = 1:num_discs
                % Calculate the new disc region
                disc_draw = mask & ((X - target_x(m)).^2 + (Y - target_y(m)).^2 <= target^2);
                % Delete the previous disc plot
                delete(findobj(gca, 'type', 'contour'));
                % Plot the new disc
                contour(X, Y, disc_draw, [0.5 0.5], 'k', 'LineWidth', 1);
            end
            pause(0.0001); % Pause to make the update visible
            counter = counter + 1 % update step counter 
    end
end
        end
        
%end
% Plot the walker's final positions
    %plot(x, y,'Color', color_palette(run,:), 'LineWidth', 1.5);
    
    % Plot the updated positions of the discs
    %for k = 1:num_discs
     %   disc_draw = ((X - target_x(k)).^2 + (Y - target_y(k)).^2 <= target^2);
      %  contour(X, Y, disc_draw, [0.5 0.5], 'k', 'LineWidth', 1); % Use a different color to distinguish from initial positions
    %end
%end
xlabel('X Position');
ylabel('Y Position');
title('Random Walk');
axis equal; % set equal aspect ratio
grid on;
hold on;

        