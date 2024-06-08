%save all figs
% Get handles to all open figures
figHandles = findall(0, 'Type', 'figure');
% Directory to save the figures
saveDir = pwd
% Loop through each figure and save it
for i = 1:length(figHandles)
    % Get the current figure
    fig = figHandles(i);
    
    % Create a filename for the figure
    figName = sprintf('expfig_%d', fig.Number);
    
    
    % Save the figure as .png
    saveas(fig, fullfile(saveDir, [figName '.png']));
    
    
    % saveas(fig, fullfile(saveDir, [figName '.pdf']));
    % exportgraphics(fig, fullfile(saveDir, [figName '.pdf']), 'ContentType', 'vector');
end