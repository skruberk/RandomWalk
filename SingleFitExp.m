clc;    % Clear the command window.
clearvars;
close all;  % Close all figs

filelist = readtable('filename.csv','VariableNamingRule', 'preserve');

% Loop over each filename
for i = 1:size(filelist, 1)
     % Get the filenames
    filename = filelist.filename{i}
    filename = strrep(filename, '"', ''); % Remove double quotes
    % import from csv
    %filename = 'hits_boundarywrap_15_1.500000e+00_10000_1000_1.csv';
    dataTable = readtable(filename, 'VariableNamingRule', 'preserve'); 
    printFile = strrep(filename, '_', ' '); %for output
    
    % Read the column of each CSV file
    data = dataTable{:, 3};
    
    % Make data numeric
    if iscell(data)
        data = cellfun(@str2double, data);
    end
    
    % Remove any NaN values
    data = data(~isnan(data));
    
    % Check if data is not empty and numeric 
    if isempty(data) || ~isnumeric(data)
        error('Data is not numeric or contains only NaNs.');
    end
    
    % Histogram
    [binCounts, binEdges] = histcounts(data);
    binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
    
    % Data fitting
    xData = binCenters(:);
    yData = binCounts(:);
    
    % Filter out zero counts to avoid fitting issues
    nonZeroIdx = yData > 0;
    xData = xData(nonZeroIdx);
    yData = yData(nonZeroIdx);
    
    % Initial guesses into vector
    expModel = @(b, x) b(1) * exp(b(2) * x) + b(3);
    initialAmplitude = max(yData);
    initialDecayRate = -0.01;
    initialOffset = 10000;
    initialGuess = [initialAmplitude, initialDecayRate, initialOffset];
    
    % Nonlinear model to fit the data, with tolerances
    options = statset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6);
    [beta, R, J, COVB, mse] = nlinfit(xData, yData, expModel, initialGuess, options); 
    %beta(1)is coefficient multiplying the exponential term
    %beta(2) is exp of the term
    %beta(3) is constant term
    
    % Fit values
    fittedYData = expModel(beta, xData);
    
    % Error check for Inf or NaN values in the fitted curve
    if any(~isfinite(fittedYData))
        disp('contains Inf or NaN values.');
    end
    
    
    % Define point range for slopes 
      % Data range for slope
    xRange = [0, 10000];
    Yfit = true(size(xData));
    
    % Fit values at endpoints 
    fittedValuesAtEndpoints = expModel(beta, xRange);
    disp('Fitted values at endpoint1:');
    
    % Slope calculation 
    slope = (fittedValuesAtEndpoints(2) - fittedValuesAtEndpoints(1)) / (xRange(2) - xRange(1));
    disp(['Slope: ', num2str(slope)]);
    
    % Derivative for the single-term exponential model
    expDx = @(b, x) b(1) * b(2) * exp(b(2) * x);
    
    % Dx xRange
    xValues = linspace(xRange(1), xRange(2), 10000);
    dxValues = expDx(beta, xValues);
    
    % Average Dx
    aveDx = mean(dxValues);
    
    % Plot histogram and fitted curve
    figure;
    bar(binCenters, binCounts, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none');
    hold on;
    plot(xData, fittedYData, 'b-', 'LineWidth', 2);
    xlabel('Steps to Hit');
    ylabel('Raw Counts');
    title(printFile);
    legend('Histogram', 'Exponential Fit');
    text(mean(xRange), mean(fittedValuesAtEndpoints), sprintf('Slope: %.4f', slope), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'blue');
    hold off;
    
    % Print fit parameters
    disp('Fitted Parameters:');
    disp(beta);
    disp('dDx:');
    disp(aveDx);
    
    exponent = beta(2)
    
    % File handling
    % Open a file for writing
    fileID = fopen('OutputExp.csv', 'a'); % 'w' for write, 'a' for append   
    % Check if the file is empty and write headers if it is
    fileInfo = dir('OutputExp.csv');
    if fileInfo.bytes == 0
        fprintf(fileID, 'filename,slope,exponent,beta1,beta2,beta3\n');
    else
        % Move to the end of the file and add a newline if the file is not empty
        fseek(fileID, 0, 'eof');
        fprintf(fileID, '\n');
    end
    
    % Write data to the CSV file
    %fprintf(fileID, '%s,%.4f,%.4f,%.4f,%.4f\n', printFile,slope,num2str(beta(1),'%.4f'),num2str(beta(2),'%.4f')num2str(beta(3),'%.4f')num2str(beta(4),'%.4f'));
    fprintf(fileID, '%s,%.4f,%.4f,%.4f,%.4f,%.4f', printFile, slope, exponent,beta(1), beta(2), beta(3));
    fclose(fileID);
end
