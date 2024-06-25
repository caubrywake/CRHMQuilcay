%% Import CRHM outputs and save as matlab format for the main simulations
% edited by Caroline Aubry-Wake, 2023-12-12
% contact at caroline.aubrywake@gmail.com

% This script pimports the CRHM outputs 9either .txt or .csv) and organises
% them and converts to a matlab format, which is used in further analysis.

%% Set-up
close all
clear all
cd 'G:\11_CRHM_cuchi\CRHM\20240527\output\' ;% set to working folder
addpath 'G:\11_CRHM_cuchi\script\' ;% local directory with CRHM model results
figdir = 'G:\11_CRHM_cuchi\fig\sensitivity\';

%% Load measured sreamflow
d = readtable('G:\11_CRHM_cuchi\data\Q\cda_lev-q-t.csv'); % Read the table from the file
d = table2timetable(d);
d = retime(d,'daily', 'mean');
Q = d.Discharge_cms_;
Qt = d.Datetime;
Qt - Qt - days(2);

% cut to 2013-2020
a = find(Qt == '26-June-2014');
Qt(1:a)=[];
Q(1:a)=[];

% set anomaly period to nan
a = find(Qt == '04-Aug-2016');
b = find(Qt == '07-Dec-2016');
Q(a:b)=nan;


%% updated
% Initialize arrays to store data and filenames
allData1 = [];
allData2 = [];

% Display the filtered file list
fileList = dir('scenarioquilcay_routing_*.txt'); % Get list of files

fileNames = strings(numel(fileList), 1);

for i = 1:numel(fileList)
    fn = fileList(i).name; % Get file name
    % Extract characters between positions 20 and 30
    extractedName = erase(strrep(fn, 'scenarioquilcay_routing_', ''), '.txt');
    fileNames(i) = extractedName; % Store the modified filename
    
    
    H = importdata(fn,' ',2); % Import headers
    D = importdata(fn); % Import data
    headers = regexp(H(1, 1), '\s+', 'split'); % Split headers
    headers = string(vertcat(headers{:})); % Split headers
    idxvar = [1, find(contains(headers,'(1)'))]; % Select the number of variables in the file
    % Time is always the first one, followed by 2 or 3 variables
    numvar = numel(idxvar); % Number of variables

    % Initialize an array to store the current file's data
    currentData = [];

    for ii = 1:numvar % For each variable, select all the column with data corresponding to that name
        varname = char(headers(idxvar(ii)));
        varname = strcat(varname(1:end-3)); % Remove hru number from name
        if varname == 't' % Exception is time - it does not have a number
            varname = 'time';
        end 
        Index = strfind(headers, varname);
        Index = find(not(cellfun('isempty', Index)));
        currentData = [currentData, D.data(:, Index)]; % Append data for the current variable
    end 

 % Append the current file's data to the allData1 and allData2 arrays
    if isempty(allData1) && numvar > 1
        allData1 = currentData(:, 2); % Second column
        if numvar > 2
            allData2 = currentData(:, 3); % Third column
        end
    else
        if numvar > 1
            allData1 = [allData1, currentData(:, 2)]; % Append second column
        end
        if numvar > 2
            allData2 = [allData2, currentData(:, 3)]; % Append third column
        end
    end
end 

% Convert the time to a datetime format
time = datetime(datevec(currentData(:, 1) + 693960));
time = dateshift(time, 'start', 'hour', 'nearest');

basinflow_lag_stor = allData1;
gwflow_lag_stor = allData2;
lag_stor_columnname = fileNames

% Save the combined data and filenames
save('Cuchi_routing_sensitivity.mat', 'basinflow_lag_stor', 'gwflow_lag_stor', 'lag_stor_columnname' );

% Clear variables
clear D H fileList fn headers i idxvar ii Index numvar pattern splitcells str varname

%% retime to daily
X = timetable (time, allData1, allData2);
XX = retime(X, 'daily', 'mean');
timed = XX.time;
allData1d =XX.allData1/3600;
allData2d = XX.allData2/3600;
%% Plot the variables
% Define color order
lineColors = [
    [173, 216, 230] / 255; % Pale blue
    [0, 0, 255] / 255;     % Mid blue
    [0, 0, 0] / 255;       % Black
    [139, 0, 0] / 255;     % Dark red
    [255, 0, 0] / 255;     % Red
    [255, 182, 193] / 255  % Light red
];

% Define factors
factors = [0.1, 0.5, 1, 1.5, 2, 10];

% Line width for plots
lw = 0.8;

% Letters for subplot titles
letters = 'abcdefghijklm';

% Remove the last 2 characters from each filename to get the base name
baseNames = extractBefore(fileNames, strlength(fileNames) - 1);
% Find unique base names
uniqueBaseNames = unique(baseNames);

% Loop through each unique base name to create a figure for each group
for k = 1:length(uniqueBaseNames)
    baseName = uniqueBaseNames{k};
    
    % Find indices of the variables in fileNames that match the current base name
    indices = find(strcmp(baseNames, baseName));
    
    % Create a new figure
    figure('Units', 'inches', 'Position', [0, 0, 8, 3]); % Set figure size to 8 inches wide, 2 inches high
    hold on; % Hold on to the current plot
    
    % Collect the plot data
    plotData = struct('total', {}, 'gw', {}, 'color', {}, 'labelTotal', {}, 'labelGw', {});
    
    for idx = 1:length(indices)
        i = indices(idx);
        
        % Extract the parameter value for the current base name and index
        param_value = factors(idx);
        
        % Store plot data
        plotData(idx).total = allData1d(:, i) + allData2d(:, i);
        plotData(idx).gw = allData2d(:, i);
        plotData(idx).color = lineColors(idx, :);
        plotData(idx).labelTotal = ['Streamflow, * ', num2str(param_value)];
        plotData(idx).labelGw = ['GW, * ', num2str(param_value)];
    end
    
    % Plot the data, ensuring the line with factor value 1 is plotted last
    for idx = [1, 2, 4, 5, 6, 3] % Order: all except the third, then the third
        plot(timed, plotData(idx).total, 'DisplayName', plotData(idx).labelTotal, 'Color', plotData(idx).color, 'LineWidth', lw);
        plot(timed, plotData(idx).gw, '--', 'DisplayName', plotData(idx).labelGw, 'Color', plotData(idx).color, 'LineWidth', lw);
    end
    
    % Plot baseline data with "Measured" in the legend
    plot(Qt, Q, ':k', 'DisplayName', 'Measured', 'LineWidth', lw);
    
    % Add title using text for left justification
    xlimits = get(gca, 'XLim');
    ylimits = get(gca, 'YLim');
    text(datetime(2016, 7, 15), ylimits(2) - 0.05 * (ylimits(2) - ylimits(1)), ...
        ['(', letters(k), ') ', strrep(baseName, '_', '\_')], ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontWeight', 'normal');

    % Add labels and legend
    ylabel('Flow (m^3 s^{-1})');

    xlim([datetime(2016, 7, 1), datetime(2019, 10, 1)]);
    box on
    % Add legend outside the plot
    legend('Location', 'eastoutside'); % Set the legend to be outside on the right side of the plot
    % Save the figure as .pdf and .png in the specified folder
    saveas(gcf, fullfile(figdir, [baseName, '.pdf']));
    saveas(gcf, fullfile(figdir, [baseName, '.png']));
    
    hold off; % Release the plot
end

%%
% Calculate NSE and KGE for each line
nse_values = zeros(size(allData1, 2), 1);
kge_values = zeros(size(allData1, 2), 1);

for i = 1:size(allData1, 2)
    % Extract simulated discharge data
    Q_sim = allData1(:, i) + allData2(:, i);
    
    % Calculate NSE
    nse_values(i) = calculateNSE(Q*3600, Qt, Q_sim, time);
    
    % Calculate KGE
    kge_values(i) = calculateKGE(Q*3600, Qt, Q_sim, time);
end



% Plot the NSE values

%% plot 
% Repeat the factors to match the length of your base names
repeatedFactors = repmat(factors, 1, length(baseNames) / length(factors));

% Add the baseline entry manually
baselineName = "baseline";
baselineFactor = 1;
baselineNseValue = nse_values(3); % Example NSE value for baseline
baselineKgeValue = kge_values(3); % Example KGE value for baseline

% Filter out the entries where the factor is 1
validIndices = repeatedFactors ~= 1;

% Filter the data
filteredBaseNames = baseNames(validIndices);
filteredFactors = repeatedFactors(validIndices)';
filteredNseValues = nse_values(validIndices);
filteredKgeValues = kge_values(validIndices);

% Create x-axis labels by combining base names with factors
xLabels = strcat(filteredBaseNames, ' * ', string(filteredFactors));

% Prepend the baseline entry to the filtered data
filteredBaseNames = [baselineName; filteredBaseNames];
filteredFactors = [baselineFactor; filteredFactors];
filteredNseValues = [baselineNseValue; filteredNseValues];
filteredKgeValues = [baselineKgeValue; filteredKgeValues];
xLabels = ['Baseline'; xLabels];
% Replace underscores with spaces in xLabels
xLabels = strrep(xLabels, '_', ' ');

% Create a new figure
figure('Units', 'inches', 'Position', [0, 0, 14 8]); % Adjust size as needed

% Plot NSE values in the first subplot
subplot(2, 1, 1);
bar(filteredNseValues, 'FaceColor', [0.2, 0.2, 0.5]); % Customize color as needed
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels, 'XTickLabelRotation', 90);
title('NSE Values');
ylabel('NSE');
grid on;

% Plot KGE values in the second subplot
subplot(2, 1, 2);
bar(filteredKgeValues, 'FaceColor', [0.5, 0.2, 0.2]); % Customize color as needed
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels, 'XTickLabelRotation', 90);
title('KGE Values');
ylabel('KGE');
grid on;

% Adjust figure for better visualization
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [12, 8]);
set(gcf, 'PaperPosition', [0, 0, 14, 8]);

% Save the figure
saveas(gcf, fullfile(figdir, 'nse_kge_bar_plots_with_baseline.pdf'));
saveas(gcf, fullfile(figdir, 'nse_kge_bar_plots_with_baseline.png'));

%% plot the difference form the baseline instead
% Define factors
factors = [0.1, 0.5, 1, 1.5, 2, 10];

% Remove the last 2 characters from each filename to get the base name
baseNames = extractBefore(fileNames, strlength(fileNames) - 1);

% Repeat the factors to match the length of your base names
repeatedFactors = repmat(factors, 1, length(baseNames) / length(factors));

% Filter out the entries where the factor is 1
validIndices = repeatedFactors ~= 1;

% Filter the data
filteredBaseNames = baseNames(validIndices);
filteredFactors = repeatedFactors(validIndices);
filteredNseValues = nse_values(validIndices);
filteredKgeValues = kge_values(validIndices);

% Calculate baseline values
baselineIndices = repeatedFactors == 1;
baselineNseValues = nse_values(baselineIndices);
baselineKgeValues = kge_values(baselineIndices);

% Create x-axis labels by combining base names with factors
xLabels = strcat(filteredBaseNames, ' * ', string(filteredFactors'));

% Calculate differences from baseline
nseDiffs = filteredNseValues - baselineNseValues(1);
kgeDiffs = filteredKgeValues - baselineKgeValues(1);

% Replace underscores with spaces in xLabels
xLabels = strrep(xLabels, '_', ' ');

% Create a new figure for differences
figure('Units', 'inches', 'Position', [0, 0,14, 8]); % Adjust size as needed

% Plot NSE differences in the first subplot
subplot(2, 1, 1);
bar(nseDiffs, 'FaceColor', [0.2, 0.2, 0.5]); % Customize color as needed
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels, 'XTickLabelRotation', 90);
title('NSE Differences from Baseline');
ylabel('NSE Difference');
grid on;

% Plot KGE differences in the second subplot
subplot(2, 1, 2);
bar(kgeDiffs, 'FaceColor', [0.5, 0.2, 0.2]); % Customize color as needed
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels, 'XTickLabelRotation', 90);
title('KGE Differences from Baseline');
ylabel('KGE Difference');
grid on;

% Adjust figure for better visualization
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [12, 8]);
set(gcf, 'PaperPosition', [0, 0, 14, 8]);

% Save the figure
saveas(gcf, fullfile(figdir, 'nse_kge_diff_from_baseline.pdf'));
saveas(gcf, fullfile(figdir, 'nse_kge_diff_from_baseline.png'));

a = numel(find(abs(nseDiffs) > 0.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% for the leakage
allData1 = [];
allData2 = [];
fileList = dir('scenarioquilcay_leak*.txt'); % Get list of files
fileNames = strings(numel(fileList), 1);


% Display the filtered file list

fileNames = strings(numel(fileList), 1);

for i = 1:numel(fileList)
    fn = fileList(i).name; % Get file name
    % Extract characters between positions 20 and 30
    extractedName = erase(strrep(fn, 'scenarioquilcay_leakage_', ''), '.txt');
    fileNames(i) = extractedName; % Store the modified filename
    
    
    H = importdata(fn,' ',2); % Import headers
    D = importdata(fn); % Import data
    headers = regexp(H(1, 1), '\s+', 'split'); % Split headers
    headers = string(vertcat(headers{:})); % Split headers
    idxvar = [1, find(contains(headers,'(1)'))]; % Select the number of variables in the file
    % Time is always the first one, followed by 2 or 3 variables
    numvar = numel(idxvar); % Number of variables

    % Initialize an array to store the current file's data
    currentData = [];

    for ii = 1:numvar % For each variable, select all the column with data corresponding to that name
        varname = char(headers(idxvar(ii)));
        varname = strcat(varname(1:end-3)); % Remove hru number from name
        if varname == 't' % Exception is time - it does not have a number
            varname = 'time';
        end 
        Index = strfind(headers, varname);
        Index = find(not(cellfun('isempty', Index)));
        currentData = [currentData, D.data(:, Index)]; % Append data for the current variable
    end 

 % Append the current file's data to the allData1 and allData2 arrays
    if isempty(allData1) && numvar > 1
        allData1 = currentData(:, 2); % Second column
        if numvar > 2
            allData2 = currentData(:, 3); % Third column
        end
    else
        if numvar > 1
            allData1 = [allData1, currentData(:, 2)]; % Append second column
        end
        if numvar > 2
            allData2 = [allData2, currentData(:, 3)]; % Append third column
        end
    end
end 

% Convert the time to a datetime format
time = datetime(datevec(currentData(:, 1) + 693960));
time = dateshift(time, 'start', 'hour', 'nearest');

basinflow_leak = allData1;
gwflow_leak = allData2;
leak_columnname = fileNames

% Save the combined data and filenames
save('Cuchi_leakage_sensitivity.mat', 'basinflow_leak', 'gwflow_leak', 'leak_columnname' );

% Clear variables
clear D H fileList fn headers i idxvar ii Index numvar pattern splitcells str varname

%% retime to daily
X = timetable (time, allData1, allData2);
XX = retime(X, 'daily', 'mean');
timed = XX.time;
allData1d =XX.allData1/3600;
allData2d = XX.allData2/3600;
%% Plot the variables

uniqueNames = [
    "abl glac";
    "acc glac";
    "all glac";
    "all";
    "base";
    "no leak, with springs";
    "no leak, no springs";
    "lake";
    "upper valley";
    "rock, lake valley, glac";
    "rock, lake, valley";
];

order =[7, 6, 2, 1, 3, 8, 9, 11, 10, 4, 5];

% Define the number of unique names
numUniqueNames = numel(uniqueNames);

% Get a color map
cmap = jet(numUniqueNames);

% Create a new figure
figure('Units', 'inches', 'Position', [0, 0, 8, 3]); % Set figure size to 8 inches wide, 3 inches high
hold on; % Hold on to the current plot
lw = 0.8;
% Loop through each unique name to plot the data according to the order
for idxIdx = 1:length(order)
    idx = order(idxIdx);
    baseName = uniqueNames{idx};
    
    % Find indices of the variables in fileNames that match the current base name
    indices = find(ismember(uniqueNames, baseName));
    
    % Determine the color: use black for the last index in the loop
    if idxIdx == length(order)
        color = [0, 0, 0]; % Black
    else
        color = cmap(idx, :); % Color from colormap
    end

    for i = indices'
        % Plot the sum of allData1 and allData2 with "total - fileName" as the DisplayName
        plot(timed, allData1d(:, i) + allData2d(:, i), 'DisplayName', uniqueNames{idx}, 'Color', color, 'LineWidth', lw);
        % Plot allData2 with "gw - fileName" as the DisplayName
        plot(timed, allData2d(:, i), '--', 'DisplayName', [uniqueNames{idx} ' GW'], 'Color', color, 'LineWidth', lw);
    end
end

% Plot baseline data with "Measured" in the legend
plot(Qt, Q, 'DisplayName', 'Measured', 'LineStyle', ':', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);

ylimits = get(gca, 'YLim');

% Add title, labels, and legend
text(datetime(2016, 7, 15), ylimits(2) - 0.05 * (ylimits(2) - ylimits(1)), ...
    ['(a)', strrep(' gw_whereto ', '_', '\_')], ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontWeight', 'normal');

ylabel('Flow (m^3 s^{-1})');

xlim([datetime(2016, 7, 1), datetime(2019, 10, 1)]);
box on
legend('Location', 'northeastoutside'); % Set the legend to be outside on the right side of the plot
hold off; % Release the plot
legendHandle = legend('Location', 'northeastoutside'); % Set the legend to be outside on the right side of the plot
set(legendHandle, 'FontSize', 8); % Shrink the legend font size
legendHandle.ItemTokenSize = [8, 8]; % Shrink legend items size


% Save the figure
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8, 3]);
set(gcf, 'PaperPosition', [0, 0, 8, 3]);

% Save the figure
saveas(gcf, fullfile(figdir, ['leakage', '.pdf']));
saveas(gcf, fullfile(figdir, ['leakage', '.png']));
%% v2
uniqueNames = [
    "abl glac";
    "acc glac";
    "all glac";
    "all";
    "base";
    "no leak, with springs";
    "no leak, no springs";
    "lake";
    "upper valley";
    "rock, lake valley, glac";
    "rock, lake, valley";
];

order = [7, 6, 2, 1, 3, 8, 9, 11, 10, 4, 5];

% Split the order array into two parts
order1 = [7, 6, 2, 1, 3, 5];
order2 = [ 8, 9, 11, 10, 4, 5];

% Define the number of unique names
numUniqueNames = numel(uniqueNames);

% Get a color map
cmap = jet(numUniqueNames);

%% Plot the first half
% Create a new figure
figure('Units', 'inches', 'Position', [0, 0, 8, 3]); % Set figure size to 8 inches wide, 3 inches high
hold on; % Hold on to the current plot
lw = 0.8;

% Loop through each unique name to plot the data according to the order
for idxIdx = 1:length(order1)
    idx = order1(idxIdx);
    baseName = uniqueNames{idx};
    
    % Find indices of the variables in fileNames that match the current base name
    indices = find(ismember(uniqueNames, baseName));
    
    % Determine the color: use black for the last index in the loop
    if idxIdx == length(order1)
        color = [0, 0, 0]; % Black
    else
        color = cmap(idx, :); % Color from colormap
    end

    for i = indices'
        % Plot the sum of allData1 and allData2 with "total - fileName" as the DisplayName
        plot(timed, allData1d(:, i) + allData2d(:, i), 'DisplayName', uniqueNames{idx}, 'Color', color, 'LineWidth', lw);
        % Plot allData2 with "gw - fileName" as the DisplayName
        plot(timed, allData2d(:, i), '--', 'DisplayName', [uniqueNames{idx} ' GW'], 'Color', color, 'LineWidth', lw);
    end
end

% Plot baseline data with "Measured" in the legend
plot(Qt, Q, 'DisplayName', 'Measured', 'LineStyle', ':', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);

ylimits = get(gca, 'YLim');

% Add title, labels, and legend
text(datetime(2016, 7, 15), ylimits(2) - 0.05 * (ylimits(2) - ylimits(1)), ...
    ['(a)', strrep(' gw_whereto ', '_', '\_')], ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontWeight', 'normal');

ylabel('Flow (m^3 s^{-1})');

xlim([datetime(2016, 7, 1), datetime(2019, 10, 1)]);
box on
legend('Location', 'northeastoutside'); % Set the legend to be outside on the right side of the plot
hold off; % Release the plot
legendHandle = legend('Location', 'northeastoutside'); % Set the legend to be outside on the right side of the plot
set(legendHandle, 'FontSize', 8); % Shrink the legend font size
legendHandle.ItemTokenSize = [8, 8]; % Shrink legend items size

% Save the figure
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8, 3]);
set(gcf, 'PaperPosition', [0, 0, 8, 3]);

saveas(gcf, fullfile(figdir, 'leakage_1.pdf'));
saveas(gcf, fullfile(figdir, 'leakage_1.png'));

%% Plot the second half
% Create a new figure
figure('Units', 'inches', 'Position', [0, 0, 8, 3]); % Set figure size to 8 inches wide, 3 inches high
hold on; % Hold on to the current plot
lw = 0.8;

% Loop through each unique name to plot the data according to the order
for idxIdx = 1:length(order2)
    idx = order2(idxIdx);
    baseName = uniqueNames{idx};
    
    % Find indices of the variables in fileNames that match the current base name
    indices = find(ismember(uniqueNames, baseName));
    
    % Determine the color: use black for the last index in the loop
    if idxIdx == length(order2)
        color = [0, 0, 0]; % Black
    else
        color = cmap(idx, :); % Color from colormap
    end

    for i = indices'
        % Plot the sum of allData1 and allData2 with "total - fileName" as the DisplayName
        plot(timed, allData1d(:, i) + allData2d(:, i), 'DisplayName', uniqueNames{idx}, 'Color', color, 'LineWidth', lw);
        % Plot allData2 with "gw - fileName" as the DisplayName
        plot(timed, allData2d(:, i), '--', 'DisplayName', [uniqueNames{idx} ' GW'], 'Color', color, 'LineWidth', lw);
    end
end

% Plot baseline data with "Measured" in the legend
plot(Qt, Q, 'DisplayName', 'Measured', 'LineStyle', ':', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);

ylimits = get(gca, 'YLim');

% Add title, labels, and legend
text(datetime(2016, 7, 15), ylimits(2) - 0.05 * (ylimits(2) - ylimits(1)), ...
    ['(b)', strrep(' gw_whereto ', '_', '\_')], ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontWeight', 'normal');

ylabel('Flow (m^3 s^{-1})');

xlim([datetime(2016, 7, 1), datetime(2019, 10, 1)]);
box on
legend('Location', 'northeastoutside'); % Set the legend to be outside on the right side of the plot
hold off; % Release the plot
legendHandle = legend('Location', 'northeastoutside'); % Set the legend to be outside on the right side of the plot
set(legendHandle, 'FontSize', 8); % Shrink the legend font size
legendHandle.ItemTokenSize = [8, 8]; % Shrink legend items size

% Save the figure
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8, 3]);
set(gcf, 'PaperPosition', [0, 0, 8, 3]);

saveas(gcf, fullfile(figdir, 'leakage_2.pdf'));
saveas(gcf, fullfile(figdir, 'leakage_2.png'));


%%
% Calculate NSE and KGE for each line
nse_values = zeros(size(allData1, 2), 1);
kge_values = zeros(size(allData1, 2), 1);

for i = 1:size(allData1, 2)
    % Extract simulated discharge data
    Q_sim = allData1(:, i) + allData2(:, i);
    
    % Calculate NSE
    nse_values(i) = calculateNSE(Q*3600, Qt, Q_sim, time);
    
    % Calculate KGE
    kge_values(i) = calculateKGE(Q*3600, Qt, Q_sim, time);
end



% Plot the NSE values

%% plot 
% Define the order of unique names
order = [7, 6, 2, 1, 3, 8, 9, 11, 10, 4, 5];

% Extract the unique names in the specified order
orderedUniqueNames = uniqueNames(order);

% Extract the NSE and KGE values in the specified order
orderedNSEValues = nse_values(order);
orderedKGEValues = kge_values(order);

% Extract the baseline values
baselineIndex = find(order == find(strcmp(uniqueNames, 'base')));
baselineNSE = orderedNSEValues(baselineIndex);
baselineKGE = orderedKGEValues(baselineIndex);

% Calculate differences from baseline
nseDifferences = orderedNSEValues - baselineNSE;
kgeDifferences = orderedKGEValues - baselineKGE;


% Create a figure for NSE and KGE values
figure('Units', 'inches', 'Position', [0, 0,8, 8]); % Adjust size as needed

subplot(2, 1, 1);
bar(orderedNSEValues, 'FaceColor', [0.2, 0.4, 0.6]); % Adjust color as needed
set(gca, 'XTick', 1:length(orderedUniqueNames), 'XTickLabel', orderedUniqueNames, 'FontSize', 14);
xtickangle(30);

ylabel('NSE Values', 'FontSize', 14);
title('(a)', 'FontSize', 14);
grid on;

subplot(2, 1, 2);
bar(orderedKGEValues, 'FaceColor',  [0.5, 0.2, 0.2]); % Adjust color as needed
set(gca, 'XTick', 1:length(orderedUniqueNames), 'XTickLabel', orderedUniqueNames, 'FontSize', 14);
xtickangle(30);
xlabel('Model Configuration', 'FontSize', 14);
ylabel('KGE Values', 'FontSize', 14);
title('(b)', 'FontSize', 14);
grid on;

% Save the figure
saveas(gcf, fullfile(figdir, 'nse_kge_values_leak.pdf'));
saveas(gcf, fullfile(figdir, 'nse_kge_values_leak.png'));

% Create a figure for differences from baseline
figure('Units', 'inches', 'Position', [0, 0,8, 8]); % Adjust size as needed
subplot(2, 1, 1);
bar(nseDifferences, 'FaceColor', [0.2, 0.4, 0.6]); % Adjust color as needed
set(gca, 'XTick', 1:length(orderedUniqueNames), 'XTickLabel', orderedUniqueNames, 'FontSize', 14);
xtickangle(30);

ylabel('Difference from baseline (NSE)', 'FontSize', 14);
title('(a)', 'FontSize', 14);
ylim ([-0.11 0.025])
grid on;

subplot(2, 1, 2);
bar(kgeDifferences, 'FaceColor',  [0.5, 0.2, 0.2]); % Adjust color as needed
set(gca, 'XTick', 1:length(orderedUniqueNames), 'XTickLabel', orderedUniqueNames, 'FontSize', 14);
xtickangle(30);
xlabel('Model Configuration', 'FontSize', 14);
ylabel('Difference from Baseline (KGE)', 'FontSize', 14);
title('(b)', 'FontSize', 14);
ylim ([-0.07 0.02])

grid on;

% Save the figure
saveas(gcf, fullfile(figdir, 'nse_kge_differences_leak.pdf'));
saveas(gcf, fullfile(figdir, 'nse_kge_differences_leak.png'));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% for the soil moisture

allData1 = [];
allData2 = [];
fileList = dir('scenarioquilcay_soilparams_*.txt'); % Get list of files
fileNames = strings(numel(fileList), 1);


% Display the filtered file list

fileNames = strings(numel(fileList), 1);

for i = 1:numel(fileList)
    fn = fileList(i).name; % Get file name
    % Extract characters between positions 20 and 30
    extractedName = erase(strrep(fn, 'scenarioquilcay_', ''), '.txt');
    fileNames(i) = extractedName; % Store the modified filename
    
    
    H = importdata(fn,' ',2); % Import headers
    D = importdata(fn); % Import data
    headers = regexp(H(1, 1), '\s+', 'split'); % Split headers
    headers = string(vertcat(headers{:})); % Split headers
    idxvar = [1, find(contains(headers,'(1)'))]; % Select the number of variables in the file
    % Time is always the first one, followed by 2 or 3 variables
    numvar = numel(idxvar); % Number of variables

    % Initialize an array to store the current file's data
    currentData = [];

    for ii = 1:numvar % For each variable, select all the column with data corresponding to that name
        varname = char(headers(idxvar(ii)));
        varname = strcat(varname(1:end-3)); % Remove hru number from name
        if varname == 't' % Exception is time - it does not have a number
            varname = 'time';
        end 
        Index = strfind(headers, varname);
        Index = find(not(cellfun('isempty', Index)));
        currentData = [currentData, D.data(:, Index)]; % Append data for the current variable
    end 

 % Append the current file's data to the allData1 and allData2 arrays
    if isempty(allData1) && numvar > 1
        allData1 = currentData(:, 2); % Second column
        if numvar > 2
            allData2 = currentData(:, 3); % Third column
        end
    else
        if numvar > 1
            allData1 = [allData1, currentData(:, 2)]; % Append second column
        end
        if numvar > 2
            allData2 = [allData2, currentData(:, 3)]; % Append third column
        end
    end
end 

% Convert the time to a datetime format
time = datetime(datevec(currentData(:, 1) + 693960));
time = dateshift(time, 'start', 'hour', 'nearest');

basinflow_soilmoist = allData1;
gwflow_soilmoist = allData2;
soilmoist_columnname = fileNames

% Save the combined data and filenames
save('Cuchi_soilmois_sensitivity.mat', 'basinflow_soilmoist', 'gwflow_soilmoist', 'soilmoist_columnname' );

% Clear variables
clear D H fileList fn headers i idxvar ii Index numvar pattern splitcells str varname

%% retime to daily
X = timetable (time, allData1, allData2);
XX = retime(X, 'daily', 'mean');
timed = XX.time;
allData1d =XX.allData1/3600;
allData2d = XX.allData2/3600;

%% Plot the variables
% Define color order

% Define color order
lineColors = [
    [173, 216, 230] / 255; % Pale blue
    [0, 0, 255] / 255;     % Mid blue
    [0, 0, 0] / 255;       % Black
    [139, 0, 0] / 255;     % Dark red
    [255, 0, 0] / 255;     % Red
    [255, 182, 193] / 255  % Light red
];

% Define factors
factors = [0.1, 0.5, 1, 1.5, 2, 10];

% Names of the variables to be plotted
baseNames = {
    'soilparams_1'
    'soilparams_2'
    'soilparams_3'
    'soilparams_4'
    'soilparams_5'
    'soilparams_6'
};

% Create a new figure
figure('Units', 'inches', 'Position', [0, 0, 8, 3]); % Set figure size to 8 inches wide, 3 inches high
hold on; % Hold on to the current plot

lw = 0.8;

% Collect the plot data
plotData = struct('Streamflow', {}, 'GW', {}, 'color', {}, 'labelTotal', {}, 'labelGw', {});

for idx = 1:length(baseNames)
    % Find index of the variable in fileNames that matches the current base name
    baseName = baseNames{idx};
    i = find(strcmp(fileNames, baseName));
    
    if ~isempty(i)
        % Extract the parameter value for the current base name and index
        param_value = factors(idx);

        % Store plot data
        plotData(idx).total = allData1d(:, i) + allData2d(:, i);
        plotData(idx).gw = allData2d(:, i);
        plotData(idx).color = lineColors(idx, :);
        plotData(idx).labelTotal = ['Streamflow, * ', num2str(param_value)];
        plotData(idx).labelGw = ['GW, * ', num2str(param_value)];
    end
end

% Plot the data, ensuring the line with factor value 1 is plotted last
plotOrder = [1, 2, 4, 5, 6, 3]; % Order: all except the third, then the third

for idx = plotOrder
    if idx <= length(plotData) && ~isempty(plotData(idx).total)
        plot(timed, plotData(idx).total, 'DisplayName', plotData(idx).labelTotal, 'Color', plotData(idx).color, 'LineWidth', lw);
        plot(timed, plotData(idx).gw, '--', 'DisplayName', plotData(idx).labelGw, 'Color', plotData(idx).color, 'LineWidth', lw);
    end
end

% Plot baseline data with "Measured" in the legend
plot(Qt, Q, ':k', 'DisplayName', 'Measured', 'LineWidth', 1);

ylimits = get(gca, 'YLim');
text(datetime(2016, 7, 15), ylimits(2) - 0.05 * (ylimits(2) - ylimits(1)), ...
        ['(a)', strrep(' soilmoist_max ', '_', '\_')], ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontWeight', 'normal');

xlim([datetime(2016, 7, 1), datetime(2019, 10, 1)]);

% Add title, labels, and legend
ylabel('Flow (m^3 s^{-1})');

box on;

% Add legend outside the plot
legend('Location', 'eastoutside');

hold off; % Release the plot

% Save the figure
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8, 3]);
set(gcf, 'PaperPosition', [0, 0, 8, 3]);

saveas(gcf, fullfile(figdir, ['soilmoist', '.pdf']));
saveas(gcf, fullfile(figdir, ['soilmoist', '.png']));
%% Calculate and plot nse and kge values

%%
% Calculate NSE and KGE for each line
nse_values = zeros(size(allData1, 2), 1);
kge_values = zeros(size(allData1, 2), 1);

for i = 1:size(allData1, 2)
    % Extract simulated discharge data
    Q_sim = allData1(:, i) + allData2(:, i);
    
    % Calculate NSE
    nse_values(i) = calculateNSE(Q*3600, Qt, Q_sim, time);
    
    % Calculate KGE
    kge_values(i) = calculateKGE(Q*3600, Qt, Q_sim, time);
end



% Plot the NSE values

%% plot 
% Repeat the factors to match the length of your base names
repeatedFactors = repmat(factors, 1, length(baseNames) / length(factors));

% Add the baseline entry manually
baselineName = "baseline";
baselineFactor = 1;
baselineNseValue = nse_values(3); % Example NSE value for baseline
baselineKgeValue = kge_values(3); % Example KGE value for baseline

% Filter out the entries where the factor is 1
validIndices = repeatedFactors ~= 1;

% Filter the data
filteredBaseNames = baseNames(validIndices);
filteredFactors = repeatedFactors(validIndices)';
filteredNseValues = nse_values(validIndices);
filteredKgeValues = kge_values(validIndices);

% Create x-axis labels by combining base names with factors
xLabels = strcat(filteredBaseNames, ' * ', string(filteredFactors));

% Prepend the baseline entry to the filtered data
filteredBaseNames = [baselineName; filteredBaseNames];
filteredFactors = [baselineFactor; filteredFactors];
filteredNseValues = [baselineNseValue; filteredNseValues];
filteredKgeValues = [baselineKgeValue; filteredKgeValues];
xLabels = ['Baseline'; xLabels];
% Replace underscores with spaces in xLabels
xLabels = strrep(xLabels, '_', ' ');

% Create a new figure
figure('Units', 'inches', 'Position', [0, 0,6,4]); % Adjust size as needed

% Plot NSE values in the first subplot
subplot(2, 1, 1);
bar(filteredNseValues, 'FaceColor', [0.2, 0.2, 0.5]); % Customize color as needed
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels, 'XTickLabelRotation', 90);
title('(a)');
ylabel('NSE');
ylim ([0 1]);
grid on;

% Plot KGE values in the second subplot
subplot(2, 1, 2);
bar(filteredKgeValues, 'FaceColor', [0.5, 0.2, 0.2]); % Customize color as needed
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels, 'XTickLabelRotation', 90);
title('(b)');
ylabel('KGE');
ylim ([0 1]);
grid on;
% Save the figure
saveas(gcf, fullfile(figdir, 'nse_kge_bar_plots_with_baseline_moisture.pdf'));
saveas(gcf, fullfile(figdir, 'nse_kge_bar_plots_with_baseline_moisture.png'));

%% plot the difference form the baseline instead
% Define factors
factors = [0.1, 0.5, 1, 1.5, 2, 10];

% Remove the last 2 characters from each filename to get the base name
baseNames = extractBefore(fileNames, strlength(fileNames) - 2);

% Repeat the factors to match the length of your base names
repeatedFactors = repmat(factors, 1, length(baseNames) / length(factors));

% Filter out the entries where the factor is 1
validIndices = repeatedFactors ~= 1;

% Filter the data
filteredBaseNames = baseNames(validIndices);
filteredFactors = repeatedFactors(validIndices);
filteredNseValues = nse_values(validIndices);
filteredKgeValues = kge_values(validIndices);

% Calculate baseline values
baselineIndices = repeatedFactors == 1;
baselineNseValues = nse_values(baselineIndices);
baselineKgeValues = kge_values(baselineIndices);

% Create x-axis labels by combining base names with factors
xLabels = strcat(filteredBaseNames, ' * ', string(filteredFactors'));

% Calculate differences from baseline
nseDiffs = filteredNseValues - baselineNseValues(1);
kgeDiffs = filteredKgeValues - baselineKgeValues(1);

% Replace underscores with spaces in xLabels
xLabels = strrep(xLabels, '_', ' ');

% Create a new figure for differences
figure('Units', 'inches', 'Position', [0, 0,6,4]); % Adjust size as needed

% Plot NSE differences in the first subplot
subplot(2, 1, 1);
bar(nseDiffs, 'FaceColor', [0.2, 0.2, 0.5]); % Customize color as needed
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels, 'XTickLabelRotation', 30);
title('(a)');
ylabel('NSE difference');
ylim([-0.5 0.05])
grid on;

% Plot KGE differences in the second subplot
subplot(2, 1, 2);
bar(kgeDiffs, 'FaceColor', [0.5, 0.2, 0.2]); % Customize color as needed
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels, 'XTickLabelRotation', 30);
title('(b)');
ylabel('KGE difference');
ylim([-0.7 0.05])
grid on;

% Save the figure
saveas(gcf, fullfile(figdir, 'nse_kge_diff_from_baseline_moisture.pdf'));
saveas(gcf, fullfile(figdir, 'nse_kge_diff_from_baseline_moisture.png'));





