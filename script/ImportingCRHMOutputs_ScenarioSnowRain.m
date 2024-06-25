%% Import CRHM outputs and save as matlab format for the main simulations
% edited by Caroline Aubry-Wake, 2023-12-12
% contact at caroline.aubrywake@gmail.com

% This script pimports the CRHM outputs 9either .txt or .csv) and organises
% them and converts to a matlab format, which is used in further analysis.

%% Set-up
close all
clear all
cd 'G:\11_CRHM_cuchi\CRHM\output\scenario\v9\SnowRain' ;% set to working folder
addpath 'G:\11_CRHM_cuchi\script\' ;% local directory with CRHM model results
figdir = 'G:\11_CRHM_cuchi\fig\scenario\';

%% Initialize arrays to store data for each simulation
hru_rain1 = [];
hru_rain2 = [];
hru_rain3 = [];
hru_rain4 = [];

hru_snow1 = [];
hru_snow2 = [];
hru_snow3 = [];
hru_snow4 = [];

SWE1 = [];
SWE2 = [];
SWE3 = [];
SWE4 = [];

time = [];

%% Display the filtered file list
fileList = dir('quilcay_noice_RainSnow_*.txt'); % Get list of files

for i = 1:numel(fileList)
    fn = fileList(i).name; % Get file name
    
    H = importdata(fn,' ',2); % Import headers
    D = importdata(fn); % Import data
    headers = regexp(H(1, 1), '\s+', 'split'); % Split headers
    headers = string(vertcat(headers{:})); % Split headers
    
    % Extract the data
    timeData = D.data(:, 1); % Time is always the first column
    hru_rain = D.data(:, contains(headers, 'hru_rain'));
    hru_snow = D.data(:, contains(headers, 'hru_snow'));
    SWE = D.data(:, contains(headers, 'SWE'));
    
    % Store time data (assuming it's the same for all simulations)
    if isempty(time)
        time = timeData;
    end
    
    % Assign data to respective variables
    switch i
        case 1
            hru_rain1 = hru_rain;
            hru_snow1 = hru_snow;
            SWE1 = SWE;
        case 2
            hru_rain2 = hru_rain;
            hru_snow2 = hru_snow;
            SWE2 = SWE;
        case 3
            hru_rain3 = hru_rain;
            hru_snow3 = hru_snow;
            SWE3 = SWE;
        case 4
            hru_rain4 = hru_rain;
            hru_snow4 = hru_snow;
            SWE4 = SWE;
    end
end
time= datetime(datevec(time+ 693960));
time = dateshift(time,'start','hour', 'nearest');

%% Now, make plot of snowfall ratio
hru_elev = [5276 5495 5020 4793 5251 4893 4290 4045 5152 5138 4762 4358 4778 4450 4976 4666 4299 4780 4298 ];
load('G:\11_CRHM_cuchi\CRHM\output\Cuchi_20230823.mat', 'SWE', 'hru_rain', 'hru_snow')

%% Calculate Total Snowfall to Precipitation Ratio for Each Simulation
% Function to calculate the ratio for each simulation
calculate_ratio = @(hru_snow, hru_rain) sum(hru_snow) ./ (sum(hru_snow) + sum(hru_rain));

% Calculate ratios for each simulation
ratio1 = calculate_ratio(hru_snow1, hru_rain1);
ratio2 = calculate_ratio(hru_snow2, hru_rain2);
ratio3 = calculate_ratio(hru_snow3, hru_rain3);
ratio4 = calculate_ratio(hru_snow4, hru_rain4);
ratiobase = calculate_ratio(hru_snow, hru_rain)

%% Plotting SNOW ration
figure;
scatter(hru_elev, ratiobase * 100, 'k', 'filled', 'DisplayName', 'Baseline');
hold on;
scatter(hru_elev, ratio3 * 100, 'filled', 'DisplayName', '+ 4\circC');
scatter(hru_elev, ratio1 * 100, 'filled', 'DisplayName', '+ 5\circC');
% Add labels and title
xlabel('Elevation (m. a.s.l.)');
ylabel('Snowfall to Precipitation Ratio (%)');
% Set x-axis ticks every 200 meters
xtickmin = floor(min(hru_elev) / 200) * 200;
xtickmax = ceil(max(hru_elev) / 200) * 200;
xticks(xtickmin:200:xtickmax);
% Add legend
legend ('Location', 'best');
grid on
box on
% Save the figure
saveas(gcf, fullfile(figdir, 'snowfall_to_precipitation_ratio_scenario.png'));

%% Plotting SWE
figure('Units', 'inches', 'Position', [0, 0, 6, 8]);

% Columns to plot
columnsToPlot = [1, 2, 3, 4];
hru_col = hru_elev(columnsToPlot);

% Labels for the subplots
labels = {'(a)', '(b)', '(c)', '(d)'};

% Loop through each column
for i = 1:length(columnsToPlot)
    subplot(4, 1, i);
    col = columnsToPlot(i);
    
    % Plot SWE data for all scenarios for the current column
    plot(time, SWE(:, col), 'DisplayName', 'Baseline');
    hold on;
    plot(time, SWE1(:, col), 'DisplayName', '+ 5\circC, +10% ppt');
    plot(time, SWE2(:, col), 'DisplayName', '+ 5\circC, +20% ppt');
    plot(time, SWE3(:, col), 'DisplayName', '+ 4\circC, +10% ppt');
    plot(time, SWE4(:, col), 'DisplayName', '+ 4\circC, +20% ppt');
    
    % Custom text in the upper left corner
    text(0.01, 0.9, sprintf('%s HRU %d: %d m.a.s.l.', labels{i}, col, hru_col(i)), 'Units', 'normalized', 'FontSize', 10);

    % Set labels
    ylabel('SWE (mm)');
    
    % Set x-axis ticks every 6 months (0.5 years)
    xticks(datetime('26-Jun-2014'):calmonths(6):datetime('01-Apr-2020'));
    xtickformat('MMM yyyy');  % Format x-axis tick labels
    
    % Add legend only on the first subplot
    if i == 1
        legend('Location', 'best');
    end
    
    % Add grid and box
    grid on;
    box on;
    
    % Set x-axis limits
    xlim([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
end

% Save the figure
saveas(gcf, fullfile(figdir, 'SWE_timeseries_HRUscenarios.png'));