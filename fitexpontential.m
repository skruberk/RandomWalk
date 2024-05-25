clc;    % Clear the command window.
clearvars;
close all;  % Close all figs

% import from csv
filename = 'steps2hit_springtarget_30_1_1000_10000_5.000000e-01.csv';
dataTable = readtable(filename, 'VariableNamingRule', 'preserve'); 
printFile = strrep(filename, '_', ' ');
% Read the column of CSV file
data = dataTable{:, 3};

% make data numeric
if iscell(data)
    data = cellfun(@str2double, data);
end

%Remove any NaN values
data = data(~isnan(data));

%data not empty and numeric 
if isempty(data) || ~isnumeric(data)
    error('Data is not numeric or contains only NaNs.');
end

%histogram
[binCounts, binEdges] = histcounts(data);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

%data fitting
xData = binCenters(:);
yData = binCounts(:);

%filter out zero counts to avoid fitting issues
nonZeroIdx = yData > 0;
xData = xData(nonZeroIdx);
yData = yData(nonZeroIdx);

%two-term exponential model, anonymous in matlab allows it to go into
%array
expModel = @(b, x) b(1) * exp(b(2) * x) + b(3) * exp(b(4) * x);  %not sure this is right 
%expModel = @(b, x) b(1) * exp(b(2) * x + b(3)); %works but not well

%initial guesses with variables, change initial offset with magnitude
initialAmplitude = max(yData);
initialDecayRate1 = -0.1; % adjust?
initialOffset1 = 1000;  % adjust on magnitude, changed from PDF to 0
initialDecayRate2 = -0.01;

%initial guesses into vector
initialGuess = [initialAmplitude, initialDecayRate1, initialOffset1, initialDecayRate2];

% nlinfit nonlinear model to fit the data, tolerances
options = statset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6);
[beta, R, J, COVB, mse] = nlinfit(xData, yData, expModel, initialGuess, options); 

%fit values
fittedYData = expModel(beta, xData);

%error check for Inf or NaN values in the fitted curve
if any(~isfinite(fittedYData))
    disp('contains Inf or NaN values.');
end
% define point range for slopes 
xRange2 = [100, 1000];  % data range for slope
xRange1 = [0, 100];
% fit values 
fittedValuesAtEndpoints1 = expModel(beta, xRange1);
disp('Fitted values at endpoint1:');
disp(fittedValuesAtEndpoints1);
fittedValuesAtEndpoints2 = expModel(beta, xRange2);
disp('Fitted values at endpoint2:');
disp(fittedValuesAtEndpoints2);
% slope calculation 
slope1 = (fittedValuesAtEndpoints1(2) - fittedValuesAtEndpoints1(1)) / (xRange1(2) - xRange1(1));
slope2 =(fittedValuesAtEndpoints2(2) - fittedValuesAtEndpoints2(1)) / (xRange2(2) - xRange2(1));
disp(['Slope: ', num2str(slope1)]);
disp(['Slope: ', num2str(slope2)]);
% plot histogram and fitted curve
figure;
bar(binCenters, binCounts, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none');
hold on;
plot(xData, fittedYData, 'r-', 'LineWidth', 2);
xlabel('Steps to Hit');
ylabel('Raw Counts');
title(printFile);
legend('Histogram', 'Exponential Fit');
hold on;
%plot slope, text, maybe export as csv
text(mean(xRange1), mean(fittedValuesAtEndpoints1), sprintf('Slope: %.4f', slope1), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'blue');
text(mean(xRange2), mean(fittedValuesAtEndpoints2), sprintf('Slope: %.4f', slope2), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'blue');
hold off;
% print fit parameters, print to plot? 
disp('Fitted Parameters:');
disp(beta);
