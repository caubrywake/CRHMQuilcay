% Import shallow well data and plot against modelled groundwater storage

% edited by Caroline Aubry-Wake, 2023-12-12
% contact at caroline.aubrywake@gmail.com

% This script provides well water level and CRHM model GW storage by the well data collected in the Quilcay valley to GW in the valley bottom HRU. comparing
% CRHM outputs to measurements. It also analyzed streamflow components,
% both montlhy and seasonally. 
% This scripts also generates F2 in the
% manuscript in addition to figures to support the analysis.
% 
% List of tables generated and exported to fig directory (as .csv): 
%	Stat_PearsonR_GWstorage_WellLevel.csv
% 
% List of figures generated and saved (as .pdf and .png)
%    F2_WaterLevel_GWstorageComparison
%    WaterLevel_GWstorageComparison_notnorm
%    WaterLevel_GWstorageComparison_notnorm_withporosit

%% Set-Up
close all
clear all

cd 'F:\11_CRHM_cuchi\'
figdir ='f:\11_CRHM_cuchi\fig\modeleval\'
addpath 'F:\11_CRHM_cuchi\functions\' 

%% Import well data
folderPath = 'data\well\';  % Specify the folder path where the CSV files are located
filePattern = fullfile(folderPath, 'gw*.csv');  % Specify the file pattern to match CSV files

csvFiles = dir(filePattern);  % Get a list of CSV files in the folder

startTime = datetime(2015, 7, 25);      % Start date for resampling
endTime = datetime(2019, 6, 25);        % End date for resampling
timearray = (startTime:days(1):endTime)';  % Creating time array


close all
figure;
for i = 1:numel(csvFiles)
    filePath = fullfile(folderPath, csvFiles(i).name);  % Get the full file path
    data = readtable(filePath);  % Read the CSV file using readtable function
    
    % Extract the desired columns
    timestamp = data.TimeStamp;  % Assuming the column name is 'TimeStamp'
    waterlevel_bgs = table2array(data(:, contains(lower(data.Properties.VariableNames), 'bgs')));  % Select columns containing 'bgs' in their names

    % Create a timetable for the current file
    fileTable = timetable(timestamp, -waterlevel_bgs);
    
     % Resample the data to a common daily timestep
    fileTable = retime(fileTable, timearray, 'fillwithmissing');
    
    % Add the well to the syncrhonize well data
    wl(:,i) = table2array(fileTable);
       
     subplot(2,3,i)
     plot(timearray, wl(:,i));
     ylabel (' Water Level (cm)')
     title(csvFiles(i).name) 
end


%% Compare with GW - normalsed

load('CRHM\output\Cuchi_20230823.mat', 'gw', 'time')
t = find(time==datetime(2014, 10, 01));
time(1:t)=[];
gwcrhm = gw(:,8);
gwcrhm(1:t)=[];

% Normalize
gwnorm = normalize(gwcrhm, 'range');
wlnorm = normalize(wl,'range');

figure('Units', 'inches', 'Position', [1, 1, 6, 3]);
plot(time, gwnorm, 'b', 'linewidth',1);hold on
plot(timearray, mean(wlnorm,2), '-k', 'linewidth',1)
plot(timearray, wlnorm, ':k')
legend ('Mod. GW storage', 'Avg Meas. Water level', 'Meas. Water Level', 'location','northoutside', 'orientation','horizontal')
ylabel ('Normalized Water Level')

figname ='F2_WaterLevel_GWstorageComparison';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%% not Normalized
T = timetable(time,gwcrhm);
TT = retime(T,timearray);
gwcrhm = TT.gwcrhm;
wl_mean = nanmean(wl, 2)*10 % in mm

% Create the figure

figure('Units', 'inches', 'Position', [1, 1, 6, 3]);

% Set the color of the left axis to black
yyaxis left
ax1 = gca; % Get the current axes handle
ax1.YColor = 'b'; % Set the Y-axis color to black
% Plot the data for the left axis
p1 = plot(timearray, (gwcrhm - max(gwcrhm)), 'b', 'linewidth', 1) % in mm
ylabel('HRU 8 Water Storage (mm)')
ylim([-38 1])

% Set the color of the right axis to black
yyaxis right
ax2 = gca; % Get the current axes handle
ax2.YColor = 'k'; % Set the Y-axis color to black

% Plot the data for the right axis
hold on
p2 = plot(timearray, (wl_mean - max(wl_mean)), '-k', 'linewidth', 1.1) % in mm
p3 = plot(timearray, (wl * 10 - max(wl_mean)), ':k') % in mm
ylabel({'Measured Water Level';'in Piexometers (mm)'})
ylim([-2000 0])
legend ([p1(1), p2(1) p3(1)], 'HRU8', 'Average piezo WL','Piezo WL', 'orientation','horizontal', 'Location','northoutside')

% expoerting correlation coefficient
r = corrcoef(wl_mean, gwcrhm)
pearsonr_gw_wl = table(r(2))
pearsonr_gw_wl.Properties.VariableNames = {'PearsonR'}
writetable(pearsonr_gw_wl, strcat(figdir, 'Stat_PearsonR_GWstorage_WellLevel.csv'))

figname ='WaterLevel_GWstorageComparison_notnorm';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%% Plot with porosity

figure('Units', 'inches', 'Position', [1, 1, 6, 3]);
poro = [0.025:0.005:0.03]
% Set the color of the left axis to black
% Plot the data for the left axis

p1 = plot(timearray, (gwcrhm - max(gwcrhm)), 'b', 'linewidth', 1) 
hold on
p2 = plot(timearray, (wl_mean - max(wl_mean))*0.025, '--k', 'linewidth', 1.1) % in mm
p2 = plot(timearray, (wl_mean - max(wl_mean))*0.030, ':k', 'linewidth', 1.1) % in mm
legend ('CRHM','Piexzometer * 0.025', 'Piezometer * 0.03', 'Orientation' ,'Horizontal','Location', 'NorthOutside')
ylabel('HRU 8 Water Storage (mm)')
ylim([-45 1])


figname ='WaterLevel_GWstorageComparison_notnorm_withporosity';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
