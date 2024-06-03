clc;    % Clear the command window.
clearvars;
close all;  % Close all figs

%hits_boundarywrap_30_5.000000e-01_10000_1000_1
%hits_boundarywrap_5_5.000000e-01_10000_1000_1
%hits_boundarywrap_20_5.000000e-01_10000_1000_1
%hits_boundarywrap_10_5.000000e-01_10000_1000_1
%hits_boundarywrap_10_2_10000_1000_1
%hits_boundarywrap_15_2_10000_1000_1
%hits_boundarywrap_20_2_10000_1000_1
%hits_boundarywrap_30_2_10000_1000_1
%hits_boundarywrap_5_2_10000_1000_1
%hits_boundarywrap_3_2_10000_1000_1
%hits_boundarywrap_3_1_10000_1000_1
%hits_boundarywrap_15_1_10000_1000_1
%hits_boundarywrap_5_1_10000_1000_1
%hits_boundarywrap_30_1_10000_1000_1
%hits_boundarywrap_20_1_10000_1000_1
%hits_boundarywrap_10_1_10000_1000_1
%hits_boundarywrap_15_2_10000_1000_1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hits_boundarywrap_3_3.000000e-01_10000_1000_1
%hits_boundarywrap_30_3.000000e-01_10000_1000_1
%hits_boundarywrap_10_3.000000e-01_10000_1000_1
%hits_boundarywrap_10_3.000000e-01_10000_1000_1.000000e-01
%hits_boundarywrap_3_3.000000e-01_10000_1000_1.000000e-01
%hits_boundarywrap_30_3.000000e-01_10000_1000_1.000000e-01
%hits_boundarywrap_30_2_10000_1000_1.000000e-01
%hits_boundarywrap_20_2_10000_1000_1.000000e-01
%hits_boundarywrap_15_2_10000_1000_1.000000e-01
%hits_boundarywrap_10_2_10000_1000_1.000000e-01
%hits_boundarywrap_5_2_10000_1000_1.000000e-01
%hits_boundarywrap_3_2_10000_1000_1.000000e-01
%hits_boundarywrap_3_1_10000_1000_1.000000e-01
%hits_boundarywrap_5_1_10000_1000_1.000000e-01
%hits_boundarywrap_10_1_10000_1000_1.000000e-01
%hits_boundarywrap_15_1_10000_1000_1.000000e-01
%hits_boundarywrap_20_1_10000_1000_1.000000e-01
%hits_boundarywrap_30_1_10000_1000_1.000000e-01
%hits_boundarywrap_30_5.000000e-01_10000_1000_1.000000e-01
%hits_boundarywrap_20_5.000000e-01_10000_1000_1.000000e-01
%hits_boundarywrap_15_5.000000e-01_10000_1000_1.000000e-01
%hits_boundarywrap_10_5.000000e-01_10000_1000_1.000000e-01
%hits_boundarywrap_5_5.000000e-01_10000_1000_1.000000e-01
%hits_boundarywrap_3_5.000000e-01_10000_1000_1.000000e-01
%hits_boundarywrap_15_5.000000e-01_10000_1000_1
%hits_boundarywrap_3_5.000000e-01_10000_1000_1


%hits_boundarywrap_10_1.500000e+00_10000_1000_1
%hits_boundarywrap_5_2.500000e+00_10000_1000_1
%hits_boundarywrap_10_2.500000e+00_10000_1000_1
%hits_boundarywrap_15_2.500000e+00_10000_1000_1
%hits_boundarywrap_20_2.500000e+00_10000_1000_1
%hits_boundarywrap_30_2.500000e+00_10000_1000_1
%hits_boundarywrap_30_1.500000e+00_10000_1000_1
%hits_boundarywrap_20_1.500000e+00_10000_1000_1
%hits_boundarywrap_15_1.500000e+00_10000_1000_1
%hits_boundarywrap_5_1.500000e+00_10000_1000_1
%
%
%
% import from csv
filename = 'hits_boundarywrap_5_1.500000e+00_10000_1000_1.csv';
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

% check where x is 0, where the fitted y values first fall below the threshold
threshold = 2.1;
belowThresholdIdx = find(fittedYData < threshold, 1);
% Check if we found a point below the threshold
if ~isempty(belowThresholdIdx)
    % If the exact point is needed, perform linear interpolation to find a more precise x value
    if belowThresholdIdx > 1
        x1 = xData(belowThresholdIdx - 1);
        x2 = xData(belowThresholdIdx);
        y1 = fittedYData(belowThresholdIdx - 1);
        y2 = fittedYData(belowThresholdIdx);
        % Linear interpolation formula
        xzero = x1 + (threshold - y1) * (x2 - x1) / (y2 - y1);
    else
    disp('The fitted y values do not cross zero within the range of xData.');
    end
end
% define point range for slopes 
xRange2 = [xzero,10000];  % data range for slope
xRange1 = [0,xzero];

%fit data break into two ranges
Yfit1 = xData >= xRange1(1) & xData <= xRange1(2);
Yfit2 = xData >= xRange2(1) & xData <= xRange2(2);
xDataRange1 = xData(Yfit1);
fittedYData1 = fittedYData(Yfit1);
xDataRange2 = xData(Yfit2);
fittedYData2 = fittedYData(Yfit2);

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
%derivative
expDx = @(b, x) b(1) * b(2) * exp(b(2) * x) + b(3) * b(4) * exp(b(4) * x);
%xEval= [50, 500]; %evaluate at those points
%Dx = expDx(beta, xEval);
% Dx xRange1
xValues1 = linspace(xRange1(1), xRange1(2), 1000);
dxValues1 = expDx(beta, xValues1);
% Dx xRange2
xValues2 = linspace(xRange2(1), xRange2(2), 1000);
dxValues2 = expDx(beta, xValues2);
% ave Dx xRange1
aveDx1 = mean(dxValues1);
% ave Dx xRange2
aveDx2 = mean(dxValues2);

% Find where the fitted y values first hit zero or cross zero
%zeroCrossingIdx = find(diff(sign(fittedYData)) ~= 0, 1);


% plot histogram and fitted curve%
figure;
bar(binCenters, binCounts, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none');
hold on;
%plot(xData, fittedYData, 'r-','Linewidth', 2)
plot(xDataRange1, fittedYData1, 'b-', 'LineWidth', 2);
hold on;
plot(xDataRange2, fittedYData2, 'r-', 'LineWidth', 2 )
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
%text(sprintf('Dx: %.4f', Dx), ...
  %  'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', 'red');
hold off;
% print fit parameters, print to plot? 
disp('Fitted Parameters:');
disp(beta);
disp('dDx:');
disp(aveDx1);
disp(aveDx2);

%%filehandling%%%%%%%%%%%%%%%%%

% Open a file for writing
fileID = fopen('OutputSlope3.csv', 'a');%w for write a for append   
% Check if the file is empty and write headers if it is
fileInfo = dir('OutputSlope3.csv');
if fileInfo.bytes == 0
    fprintf(fileID, 'filename,slope1,slope2,beta\n');
else
    % Move to the end of the file and add a newline if the file is not empty
    fseek(fileID, 0, 'eof');
    fprintf(fileID, '\n');
end

%fileID = fopen(filename, 'a');
% Write data to the CSV file
fprintf(fileID, '%s,%.4f,%.4f,%s\n', printFile, slope1, slope2, num2str(beta, '%.4f'));

%fprintf(fileID, '%s,%.4f,%.4f,%s\n\n', printFile, slope1, slope2, num2str(beta, '%.4f '));
fclose(fileID);
