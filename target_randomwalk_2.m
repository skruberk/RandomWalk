% a random walk in 2D with reflective boundaries and one target 
clc;    % Clear the command window.
clearvars;
close all;  % Close all figures (except those of imtool.)
workspace;  % Make sure the workspace panel is showing.
% specify the number of time steps in our random walk
tsteps=50000
% here is the time increment between steps
% we will assume that this is in units of seconds
tau=0.01
% here is the diffusion coefficient for the walk
% let's assume units of um^2/sec
D=10
% now, from the diffusion coefficient and the 
% time increment we can calculate the appropriate
% distance increment (in units of um)
delta=sqrt(2*D*tau);

% make the cell with circular disc
radius = 10;    %cell radius variable
target=radius/1.5;   %membrane disc variable size
pos_outer = [-radius, -radius, 2*radius, 2*radius]; 
rectangle('Position',pos_outer,'Curvature',[1 1], 'LineWidth', 3)
size_inner = [-target, -target, 2*target, 2*target];
%make these coordinates a variable 
viscircles([radius,0], target, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 3);
axis square;
ax = gca;
ax.FontWeight = 'normal';
ax.FontSize = 8;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold on;

% getting all points for inner circle
step_size = 1;  % step size for grid
theta = linspace(0, 2*pi, tsteps);  % angle range for points 
% Calculate the coordinates of the points within the inner circle
x_inner = radius + target * cos(theta);  %change coordinate to variable
y_inner = 0 + target * sin(theta);  %change coordinate to variable

%array for time info (not currently used but could be for time to target)
time=0:tau:(tsteps*tau);


%for i = 2:tsteps
%delx = 2 * (randi(2) - 1) - 1;
  %  dely = 2 * (randi(2) - 1) - 1;
   % y(i) = y(i - 1) + dely * delta;
    %x(i) = x(i - 1) + delx * delta;
%end


%how many times to run code
num_runs = 100;
results = zeros(1, num_runs);  % array for counter results

for run = 1:num_runs
    counter = 0 %initialize target hitting counter 
    %initialize arrays to hold the x and y positions
    x=zeros(tsteps,1);
    y=zeros(tsteps,1);
    x(1) = 0.0;
    y(1) = 0.0;


    for i = 2 :tsteps
	    % Check if particle hits boundary in x
        delx = 2 * (randi(2) - 1) - 1;   %random step 
	    xbound = x(i-1) + delx*delta;    % new position 
        dely = 2 * (randi(2) - 1) - 1;
        ybound= y(i-1) + dely*delta;
         %checks if at target
	    if any(abs(x - x_inner(i)) <= step_size/3) && any(abs(y - y_inner(i)) <= step_size/3)
            %any(sqrt((x - x_inner).^2 + (y - y).^2) <= (step_size));
            disp('random walk has hit target!');
            disp(num2str(counter))
            break;
        elseif  norm([xbound,ybound]) >= radius 
            x(i)=x(i-2);
            y(i)=y(i-2);
            continue;          
        else
            x(i) = xbound; %update pos
            y(i) = ybound;
            counter = counter + 1
        end
    end
    results(run) = counter;  % Store the final counter value for this run
end

% plot each random walk
%plot(x, y);
%xlabel('X Position');
%ylabel('Y Position');
%title('Random Walk');
%axis equal; % set equal aspect ratio
%grid on;


% Plot histogram of counter results
figure;
histogram(results);
xlim([0, 50000])  %change to tsteps
xlabel('Step Number');
ylabel('Frequency');
title(radius);
xticks(get(gca, 'XTick'));  % Ensure all tick values are visible
xticklabels(arrayfun(@(x) sprintf('%d', x), get(gca, 'XTick'), 'UniformOutput', false));










