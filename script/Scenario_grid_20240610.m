%% Analysing the scenario runs

%%%%%%%%% FOR GWand surface water, not wet and dry seaosn???

% 9 subplot: GW, SW, and total, for annual values, not wet vs dry

% edited by Caroline Aubry-Wake, 2023-12-12
% contact at caroline.aubrywake@gmail.com

% This script compiles changes in streamflow and groundwater in the case of
% the precipitation and temeperature sensitivity scenarios. 
% 
% It generates figures (as pdf and png)
% 	scenario
% 	monthlyFlow_NoIce
% 	Scenario_NoIce_6subplot
% 	Scenario_25Ice_6subplot
% 	Scenario_2Ice_6subplot
% 	Scenario_GridPercentChange_WetDry_6subplot
% 	F4_Scenario_GridNumberChange_WetDry_6subplot

%% Set-up
close all
clear all
cd G:\11_CRHM_cuchi\
figdir ='fig\scenario\'
folderPath = 'CRHM\output\scenario\v9\';  % Replace with the actual folder path
mainpath = ' G:\11_CRHM_cuchi\'
%  glacier per scenario
hruarea = [1.528 1.322 1.237 0.839 7.366 4.056 3.492 3.663 3.201 4.633 9.491 2.721 2.418 2.008 6.249 10.62 ...
1.903 0.14 0.45];
sumarea = sum(hruarea)
hru_glaciercover = [sum(hruarea([2,3,5,6,9])./sumarea) 0, hruarea(2)./sumarea hruarea(5)./sumarea sum(hruarea([2,5])./sumarea)]

%% import the current simulations
load('CRHM\output\v9\Cuchi_20230823.mat', 'basinflow','basingw','time')
T = timetable(time, basinflow/3600, basingw/3600);
TT = retime(T, 'daily','mean');
b = TT.Var1;
g = TT.Var2;
t = TT.time;
s = b+g;

% import measured flow
d = readtable('data\Q\cda_lev-q-t.csv'); % Read the table from the file
d = table2timetable(d);
d = retime(d,'daily', 'mean');
Q = d.Discharge_cms_;
Qt = d.Datetime;

% cut to 2013-2020
a = find(Qt == '26-June-2014');
Qt(1:a)=[];
Q(1:a)=[];

% set that weird period to nan
t1 = find(Qt == '04-Aug-2016');
t2 = find(Qt == '07-Dec-2016');
Q(t1:t2)=nan;

%% Import scenarios
% Import sceanrios the names of the files
files = dir(fullfile(folderPath, '*.txt'));  % Replace '*.txt' with the appropriate file extension
% Extract file names and remove the first 8 letters
fileName = {files.name};
fileName = extractBetween(fileName, 9, 28);


clear B G S snowmelt firnmelt icemelt gw ssr hru_actet
for i = 1:numel(files)
fn = files(i).name; % file to import
H = importdata(strcat(folderPath, fn),' ',2); %  % import headers 
D = importdata(strcat(folderPath, fn)) ; % import data
B(:,i)= D.data(:,2);
G(:,i)= D.data(:,3);
S(:,i)= B(:,i)+G(:,i);
if i == 1
time = D.data(:,1);
end 
end

%convert the time to a datetime format
time= datetime(datevec(time+ 693960));
time = dateshift(time,'start','hour', 'nearest');
clear D H fileList fn i 

% Save scenarios at daily time step
save ('ImportedScenarios_hourly_v9.mat', 'time', 'B', 'G', 'S', 'b','g','s');

% Retime everything to daily averages
T = timetable(time, B/3600, G/3600, S/3600);
TT = retime(T, 'daily','mean');
B = TT.Var1;
G = TT.Var2;
S = TT.Var3;
timem = TT.time;

save ('ImportedScenarios_daily_v9.mat', 'timem', 'B', 'G', 'S', 'b','g','s', 'fileName')

%% compile numbers
close all
% load scenarios
load('ImportedScenarios_daily_v9.mat')

% set all values before to nan
a = find(timem == datetime('26-Jun-2014'));
b(1:a)=[];
B(1:a, :)=[];
g(1:a)=[];
G(1:a, :)=[];
s(1:a)=[];
S(1:a, :)=[];
timem(1:a)=[];


months = month(timem);
% Find monthly values
for i = 1:12
mthIndices = find(months == i);
% Extract the monthly data and average
Sm(i, :) = mean(S(mthIndices,:));
sm(i, :) = mean(s(mthIndices,:));

Gm(i, :) = mean(G(mthIndices,:));
gm(i, :) = mean(g(mthIndices,:));

Bm(i, :) = mean(B(mthIndices,:));
bm(i, :) = mean(b(mthIndices,:));
end 

wetseason = [1:3];
dryseason = [6:8];
annseason = [1:12];

Sm_wet = mean(Sm(wetseason, :));
Sm_dry = mean(Sm(dryseason, :));
Sm_ann = mean(Sm(annseason, :));

sm_dry = mean(sm(dryseason, :));
sm_wet = mean(sm(wetseason, :));
sm_ann = mean(sm(annseason, :));

Bm_wet = mean(Bm(wetseason, :));
Bm_dry = mean(Bm(dryseason, :));
Bm_ann = mean(Bm(annseason, :));

bm_dry = mean(bm(dryseason, :));
bm_wet = mean(bm(wetseason, :));
bm_ann = mean(bm(annseason, :));

Gm_wet = mean(Gm(wetseason, :));
Gm_dry = mean(Gm(dryseason, :));
Gm_ann = mean(Gm(annseason, :));

gm_dry = mean(gm(dryseason, :));
gm_wet = mean(gm(wetseason, :));
gm_ann = mean(gm(annseason, :));


%% Compile numbers for streamflow change for all simulations

% make a grid
temperatureChanges = 0:1:5;   % Temperature changes from 0 to +4°C
precipitationChanges = 80:10:120;   % Precipitation changes from -20% to +10%
numTemperatureChanges = numel(temperatureChanges);
numPrecipitationChanges = numel(precipitationChanges);

% Streamflow
streamflowTable_Sm_dry_noice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Sm_wet_noice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Sm_ann_noice = zeros(numTemperatureChanges, numPrecipitationChanges);

streamflowTable_Sm_dry_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Sm_wet_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Sm_ann_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);

streamflowTable_Sm_dry_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Sm_wet_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Sm_ann_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);

% GW
streamflowTable_Gm_dry_noice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Gm_wet_noice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Gm_ann_noice = zeros(numTemperatureChanges, numPrecipitationChanges);

streamflowTable_Gm_dry_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Gm_wet_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Gm_ann_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);

streamflowTable_Gm_dry_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Gm_wet_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Gm_ann_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);

% SW
streamflowTable_Bm_dry_noice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Bm_wet_noice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Bm_ann_noice = zeros(numTemperatureChanges, numPrecipitationChanges);

streamflowTable_Bm_dry_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Bm_wet_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Bm_ann_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);

streamflowTable_Bm_dry_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Bm_wet_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Bm_ann_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);


%% Iterate over the filenames and populate the tables
% for streamflow
for i = 1:numel(fileName)
    filename = fileName{i};
    
    % Check if the filename contains '25glac'
    if contains(filename, '25glac')
        glacString = '25glac';
    % Check if the filename contains '2glac'
    elseif contains(filename, '2glac')
        glacString = '2glac';
    % Check if the filename contains 'noice'
    elseif contains(filename, 'noice')
        glacString = 'noice';
    end
       
    % Extract the temperature change and precipitation change values from the filename
    temperatureChange = str2double(regexp(filename, '(?<=_t_)\d+', 'match'));
    precipitationChange = str2double(regexp(filename, '(?<=precip_)\d+', 'match'));

    if precipitationChange < 70 
        precipitationChange = precipitationChange*10;
    end 

    % Find the indices for the temperature and precipitation changes in the tables
    temperatureIndex = temperatureChange + 1;
    precipitationIndex = find(precipitationChanges == precipitationChange);

    % Calculate the streamflow changes for the scenarios
    streamflowChange_Sm_dry = Sm_dry(i) - sm_dry;
    streamflowChange_Sm_wet = Sm_wet(i) - sm_wet;
    streamflowChange_Sm_ann = Sm_ann(i) - sm_ann;

    streamflowChange_Bm_dry = Bm_dry(i) - bm_dry;
    streamflowChange_Bm_wet = Bm_wet(i) - bm_wet;
    streamflowChange_Bm_ann = Bm_ann(i) - bm_ann;

    streamflowChange_Gm_dry = Gm_dry(i) - gm_dry;
    streamflowChange_Gm_wet = Gm_wet(i) - gm_wet;
    streamflowChange_Gm_ann = Gm_ann(i) - gm_ann;

    % Calculate the streamflow changes as percentages
    streamflowChange_Sm_dry_p = Sm_dry(i) * 100 / sm_dry;
    streamflowChange_Sm_wet_p = Sm_wet(i) * 100 / sm_wet;
    streamflowChange_Sm_ann_p = Sm_ann(i) * 100 / sm_ann;

    streamflowChange_Bm_dry_p = Bm_dry(i) * 100 / bm_dry;
    streamflowChange_Bm_wet_p = Bm_wet(i) * 100 / bm_wet;
    streamflowChange_Bm_ann_p = Bm_ann(i) * 100 / bm_ann;

    streamflowChange_Gm_dry_p = Gm_dry(i) * 100 / gm_dry;
    streamflowChange_Gm_wet_p = Gm_wet(i) * 100 / gm_wet;
    streamflowChange_Gm_ann_p = Gm_ann(i) * 100 / gm_ann;
    
    % Update the tables with the streamflow change values based on the glacString
    switch glacString
        case '25glac'
            streamflowTable_Sm_dry_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry_p;
            streamflowTable_Sm_wet_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet_p;
            streamflowTable_Sm_ann_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_ann_p;

            streamflowTable_Sm_dry_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry;
            streamflowTable_Sm_wet_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet;
            streamflowTable_Sm_ann_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_ann;

            streamflowTable_Bm_dry_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry_p;
            streamflowTable_Bm_wet_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet_p;
            streamflowTable_Bm_ann_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_ann_p;

            streamflowTable_Bm_dry_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry;
            streamflowTable_Bm_wet_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet; 
            streamflowTable_Bm_ann_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_ann; 

            streamflowTable_Gm_dry_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry_p;
            streamflowTable_Gm_wet_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet_p;
            streamflowTable_Gm_ann_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_ann_p;

            streamflowTable_Gm_dry_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry;
            streamflowTable_Gm_wet_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet;
            streamflowTable_Gm_ann_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_ann;

        case '2glac'
            streamflowTable_Sm_dry_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry_p;
            streamflowTable_Sm_wet_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet_p;
            streamflowTable_Sm_ann_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_ann_p;

            streamflowTable_Sm_dry_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry;
            streamflowTable_Sm_wet_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet;
            streamflowTable_Sm_ann_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_ann;

            streamflowTable_Bm_dry_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry_p;
            streamflowTable_Bm_wet_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet_p;
            streamflowTable_Bm_ann_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_ann_p;

            streamflowTable_Bm_dry_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry;
            streamflowTable_Bm_wet_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet; 
            streamflowTable_Bm_ann_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_ann;

            streamflowTable_Gm_dry_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry_p;
            streamflowTable_Gm_wet_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet_p;
            streamflowTable_Gm_ann_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_ann_p;

            streamflowTable_Gm_dry_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry;
            streamflowTable_Gm_wet_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet;
            streamflowTable_Gm_ann_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_ann;

        case 'noice'
            streamflowTable_Sm_dry_noice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry_p;
            streamflowTable_Sm_wet_noice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet_p;
            streamflowTable_Sm_ann_noice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_ann_p;

            streamflowTable_Sm_dry_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry;
            streamflowTable_Sm_wet_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet;
            streamflowTable_Sm_ann_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_ann;

            streamflowTable_Bm_dry_noice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry_p;
            streamflowTable_Bm_wet_noice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet_p;
            streamflowTable_Bm_ann_noice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_ann_p;

            streamflowTable_Bm_dry_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry;
            streamflowTable_Bm_wet_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet; 
            streamflowTable_Bm_ann_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_ann;

            streamflowTable_Gm_dry_noice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry_p;
            streamflowTable_Gm_wet_noice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet_p;
            streamflowTable_Gm_ann_noice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_ann_p;

            streamflowTable_Gm_dry_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry;
            streamflowTable_Gm_wet_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet;
            streamflowTable_Gm_ann_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_ann;
    end
end

%% For B, G, and S, annual

close all
figure('Units', 'inches', 'Position',[1, 1, 11, 7]);

darkblue = [10 46 146]./255;
lightblue = [80 110 193]./255;
lightred = [191/255 84/255 84/255];
darkred = [107/255, 0, 0];

% Define custom colormap points
colorPoints = [lightred; 1, 1, 1; lightblue]; % Blue-White-Red
centerPosition = 0.5; % Define the position of the center color (white)
x = linspace(0, 1, 256); % Create a vector of color indices
newPositions = [0, centerPosition, 1]; % Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x); % Interpolate colormap using interp1

% Plotting for 'noice'
% Subplot 1: Surface water (B)
sb1 = subplot(3, 3, 1);
imagesc(streamflowTable_Bm_ann_noice-100);
title('(a) No Ice, Ov + Vad','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb1, customColormap);
caxis(sb1,[50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Bm_ann_noice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 2: Groundwater (G)
sb2 = subplot(3, 3, 2);
imagesc(streamflowTable_Gm_ann_noice-100);
title('(b) No Ice, GW','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb2, customColormap);
caxis(sb2, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Gm_ann_noice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 3: Streamflow annual (S)
sb3 = subplot(3, 3, 3);
imagesc(streamflowTable_Sm_ann_noice-100);
title('(c) No Ice, Flow','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb3, customColormap);
caxis(sb3, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_noice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Plotting for '2glac'
% Subplot 4: Surface water (B)
sb4 = subplot(3, 3, 4);
imagesc(streamflowTable_Bm_ann_2ice-100);
title('(d) 2% GC, Ov + Vad','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb4, customColormap);
caxis(sb4, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Bm_ann_2ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 5: Groundwater (G)
sb5 = subplot(3, 3, 5);
imagesc(streamflowTable_Gm_ann_2ice-100);
title('(e) 2% GC, GW','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb5, customColormap);
caxis(sb5, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Gm_ann_2ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 6: Streamflow annual (S)
sb6 = subplot(3, 3, 6);
imagesc(streamflowTable_Sm_ann_2ice-100);
title('(f) 2% GC, Flow','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb6, customColormap);
caxis(sb6, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_2ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Plotting for '25glac'
% Subplot 7: Surface water (B)
sb7 = subplot(3, 3, 7);
imagesc(streamflowTable_Bm_ann_25ice-100);
title('(g) 13% GC, Ov + Vad','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb7, customColormap);
caxis(sb7, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Bm_ann_25ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 8: Groundwater (G)
sb8 = subplot(3, 3, 8);
imagesc(streamflowTable_Gm_ann_25ice-100);
title('(h) 13% GC, GW','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb8, customColormap);
caxis(sb8, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Gm_ann_25ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 9: Streamflow annual (S)
sb9 = subplot(3, 3, 9);
imagesc(streamflowTable_Sm_ann_25ice-100);
title('(i) 13% GC, Flow','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb9, customColormap);
caxis(sb9, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_25ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end


% Set common properties for all subplots
 precipitationChangesD = [-20 -10 0 10 20]
  temperatureChangesD = [0 1 2 3 4 5]
for i = 1:9
    subplot(3, 3, i);
    set(gca, 'XTick', 1:numPrecipitationChanges, 'XTickLabel', precipitationChangesD);
    set(gca, 'YTick', 1:numTemperatureChanges, 'YTickLabel', temperatureChangesD);
    xlabel('Δ Precip (%)');
    ylabel('Δ Temp (°C)');

end

%% With percentages

close all
figure('Units', 'inches', 'Position', [1, 1, 11, 7]);

darkblue = [10 46 146]./255;
lightblue = [80 110 193]./255;
lightred = [191/255 84/255 84/255];
darkred = [107/255, 0, 0];

% Define custom colormap points
colorPoints = [lightred; 1, 1, 1; lightblue]; % Blue-White-Red
centerPosition = 0.5; % Define the position of the center color (white)
x = linspace(0, 1, 256); % Create a vector of color indices
newPositions = [0, centerPosition, 1]; % Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x); % Interpolate colormap using interp1

% Plotting for 'noice'
% Subplot 1: Surface water (B)
sb1 = subplot(3, 3, 1);
imagesc(streamflowTable_Bm_ann_noice-100);
title('(a) No Ice, Ov + Vad','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb1, customColormap);
caxis(sb1, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Bm_ann_noice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 2: Groundwater (G)
sb2 = subplot(3, 3, 2);
imagesc(streamflowTable_Gm_ann_noice-100);
title('(b) No Ice, GW','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb2, customColormap);
caxis(sb2, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Gm_ann_noice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 3: Streamflow annual (S)
sb3 = subplot(3, 3, 3);
imagesc(streamflowTable_Sm_ann_noice-100);
title('(c) No Ice, Streamflow','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb3, customColormap);
caxis(sb3, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_noice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Plotting for '2glac'
% Subplot 4: Surface water (B)
sb4 = subplot(3, 3, 4);
imagesc(streamflowTable_Bm_ann_2ice-100);
title('(d) 2% GC, Ov + Vad','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb4, customColormap);
caxis(sb4, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Bm_ann_2ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 5: Groundwater (G)
sb5 = subplot(3, 3, 5);
imagesc(streamflowTable_Gm_ann_2ice-100);
title('(e) 2% GC, GW','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb5, customColormap);
caxis(sb5, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Gm_ann_2ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 6: Streamflow annual (S)
sb6 = subplot(3, 3, 6);
imagesc(streamflowTable_Sm_ann_2ice-100);
title('(f) 2% GC, Streamflow','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb6, customColormap);
caxis(sb6, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_2ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Plotting for '25glac'
% Subplot 7: Surface water (B)
sb7 = subplot(3, 3, 7);
imagesc(streamflowTable_Bm_ann_25ice-100);
title('(g) 13% GC, Ov + Vad','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb7, customColormap);
caxis(sb7, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Bm_ann_25ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 8: Groundwater (G)
sb8 = subplot(3, 3, 8);
imagesc(streamflowTable_Gm_ann_25ice-100);
title('(h) 13% GC, GW','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb8, customColormap);
caxis(sb8, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Gm_ann_25ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 9: Streamflow annual (S)
sb9 = subplot(3, 3, 9);
imagesc(streamflowTable_Sm_ann_25ice-100);
title('(i) 13% GC, Streamflow','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb9, customColormap);
caxis(sb9, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_25ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end



% Set common properties for all subplots
 precipitationChangesD = [-20 -10 0 10 20]
for i = 1:9
    subplot(3, 3, i);
    set(gca, 'XTick', 1:numPrecipitationChanges, 'XTickLabel', precipitationChangesD);
    set(gca, 'YTick', 1:numTemperatureChanges, 'YTickLabel', temperatureChanges);
    xlabel('Δ Precip (%)');
    ylabel('Δ Temp (°C)');

end
figname ='Scenario_GridNumberChangePercent_SW_GW_Ann_subplot';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))


%%

close all
figure('Units', 'inches', 'Position', [1, 1, 11, 7]);

darkblue = [10 46 146]./255;
lightblue = [80 110 193]./255;
lightred = [191/255 84/255 84/255];
darkred = [107/255, 0, 0];

% Define custom colormap points
colorPoints = [lightred; 1, 1, 1; lightblue]; % Blue-White-Red
centerPosition = 0.5; % Define the position of the center color (white)
x = linspace(0, 1, 256); % Create a vector of color indices
newPositions = [0, centerPosition, 1]; % Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x); % Interpolate colormap using interp1

% Plotting for 'noice'
% Subplot 1: Surface water (B)
sb1 = subplot(3, 3, 1);
imagesc(streamflowTable_Bm_ann_noice-100);
title('(a) No Ice, Ov + Vad','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb1, customColormap);
caxis(sb1, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Bm_ann_noice_n(i, j),2)),'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 2: Groundwater (G)
sb2 = subplot(3, 3, 2);
imagesc(streamflowTable_Gm_ann_noice-100);
title('(b) No Ice, GW','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb2, customColormap);
caxis(sb2, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Gm_ann_noice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 3: Streamflow annual (S)
sb3 = subplot(3, 3, 3);
imagesc(streamflowTable_Sm_ann_noice-100);
title('(c) No Ice, Streamflow','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb3, customColormap);
caxis(sb3, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_noice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Plotting for '2glac'
% Subplot 4: Surface water (B)
sb4 = subplot(3, 3, 4);
imagesc(streamflowTable_Bm_ann_2ice-100);
title('(d) 2% GC, Ov + Vad','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb4, customColormap);
caxis(sb4, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Bm_ann_2ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 5: Groundwater (G)
sb5 = subplot(3, 3, 5);
imagesc(streamflowTable_Gm_ann_2ice-100);
title('(e) 2% GC, GW','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb5, customColormap);
caxis(sb5, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Gm_ann_2ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 6: Streamflow annual (S)
sb6 = subplot(3, 3, 6);
imagesc(streamflowTable_Sm_ann_2ice-100);
title('(f) 2% GC, Streamflow','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb6, customColormap);
caxis(sb6, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_2ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Plotting for '25glac'
% Subplot 7: Surface water (B)
sb7 = subplot(3, 3, 7);
imagesc(streamflowTable_Bm_ann_25ice-100);
title('(g) 13% GC, Ov + Vad','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb7, customColormap);
caxis(sb7, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Bm_ann_25ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 8: Groundwater (G)
sb8 = subplot(3, 3, 8);
imagesc(streamflowTable_Gm_ann_25ice-100);
title('(h) 13% GC, GW','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb8, customColormap);
caxis(sb8, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Gm_ann_25ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 9: Streamflow annual (S)
sb9 = subplot(3, 3, 9);
imagesc(streamflowTable_Sm_ann_25ice-100);
title('(i) 13% GC, Streamflow','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb9, customColormap);
caxis(sb9, [50 150]-100);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_25ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end



% Set common properties for all subplots
 precipitationChangesD = [-20 -10 0 10 20]
for i = 1:9
    subplot(3, 3, i);
    set(gca, 'XTick', 1:numPrecipitationChanges, 'XTickLabel', precipitationChangesD);
    set(gca, 'YTick', 1:numTemperatureChanges, 'YTickLabel', temperatureChanges);
    xlabel('Δ Precip (%)');
    ylabel('Δ Temp (°C)');

end
figname ='Scenario_GridNumberChange_SW_GW_Ann_subplot';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))


%% For wet and dry and annual streamflow 
close all
figure('Units', 'inches', 'Position', [1, 1, 11, 7]);

darkblue = [10 46 146]./255;
lightblue = [80 110 193]./255;
lightred = [191/255 84/255 84/255];
darkred = [107/255, 0, 0];

% Define custom colormap points
colorPoints = [lightred; 1, 1, 1; lightblue]; % Blue-White-Red
centerPosition = 0.5; % Define the position of the center color (white)
x = linspace(0, 1, 256); % Create a vector of color indices
newPositions = [0, centerPosition, 1]; % Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x); % Interpolate colormap using interp1

% Define common properties for the color scale
caxis_values = [50 150]-100;

% Define common properties for all subplots
precipitationChangesD = [-20 -10 0 10 20];
temperatureChanges = [0 1 2 3 4 5];
numTemperatureChanges = length(temperatureChanges);
numPrecipitationChanges = length(precipitationChangesD);

% Plotting for 'noice'
% Subplot 1: Wet season
sb1 = subplot(3, 3, 1);
imagesc(streamflowTable_Sm_wet_noice-100);
title('(a) No Ice, Wet Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb1, customColormap);
caxis(sb1, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_noice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 2: Dry season
sb2 = subplot(3, 3, 2);
imagesc(streamflowTable_Sm_dry_noice-100);
title('(b) No Ice, Dry Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb2, customColormap);
caxis(sb2, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_noice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 3: Annual
sb3 = subplot(3, 3, 3);
imagesc(streamflowTable_Sm_ann_noice-100);
title('(c) No Ice, Annual','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb3, customColormap);
caxis(sb3, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_noice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Plotting for '2glac'
% Subplot 4: Wet season
sb4 = subplot(3, 3, 4);
imagesc(streamflowTable_Sm_wet_2ice-100);
title('(d) 2% GC, Wet Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb4, customColormap);
caxis(sb4, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_2ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 5: Dry season
sb5 = subplot(3, 3, 5);
imagesc(streamflowTable_Sm_dry_2ice-100);
title('(e) 2% GC, Dry Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb5, customColormap);
caxis(sb5, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_2ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 6: Annual
sb6 = subplot(3, 3, 6);
imagesc(streamflowTable_Sm_ann_2ice-100);
title('(f) 2% GC, Annual','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb6, customColormap);
caxis(sb6, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_2ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Plotting for '25glac'
% Subplot 7: Wet season
sb7 = subplot(3, 3, 7);
imagesc(streamflowTable_Sm_wet_25ice-100);
title('(g) 13% GC, Wet Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb7, customColormap);
caxis(sb7, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_25ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 8: Dry season
sb8 = subplot(3, 3, 8);
imagesc(streamflowTable_Sm_dry_25ice-100);
title('(h) 13% GC, Dry Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb8, customColormap);
caxis(sb8, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_25ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 9: Annual
sb9 = subplot(3, 3, 9);
imagesc(streamflowTable_Sm_ann_25ice-100);
title('(i) 13% GC, Annual','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb9, customColormap);
caxis(sb9, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_25ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Set common properties for all subplots
for i = 1:9
    subplot(3, 3, i);
    set(gca, 'XTick', 1:numPrecipitationChanges, 'XTickLabel', precipitationChangesD);
    set(gca, 'YTick', 1:numTemperatureChanges, 'YTickLabel', temperatureChanges);
    xlabel('Δ Precip (%)');
    ylabel('Δ Temp (°C)');
end

figname ='Scenario_GridNumberChangePercent_WetDryAnn_subplot';
saveas(gcf, strcat(figdir, figname, '.pdf'));
saveas(gcf, strcat(figdir, figname, '.png'));
savefig(gcf, strcat(figdir, figname));


%% with  number
close all
figure('Units', 'inches', 'Position', [1, 1, 11, 7]);

darkblue = [10 46 146]./255;
lightblue = [80 110 193]./255;
lightred = [191/255 84/255 84/255];
darkred = [107/255, 0, 0];

% Define custom colormap points
colorPoints = [lightred; 1, 1, 1; lightblue]; % Blue-White-Red
centerPosition = 0.5; % Define the position of the center color (white)
x = linspace(0, 1, 256); % Create a vector of color indices
newPositions = [0, centerPosition, 1]; % Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x); % Interpolate colormap using interp1

% Define common properties for the color scale
caxis_values = [50 150]-100;

% Define common properties for all subplots
precipitationChangesD = [-20 -10 0 10 20];
temperatureChanges = [0 1 2 3 4 5];
numTemperatureChanges = length(temperatureChanges);
numPrecipitationChanges = length(precipitationChangesD);

% Plotting for 'noice'
% Subplot 1: Wet season
sb1 = subplot(3, 3, 1);
imagesc(streamflowTable_Sm_wet_noice-100);
title('(a) No Ice, Wet Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb1, customColormap);
caxis(sb1, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_noice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 2: Dry season
sb2 = subplot(3, 3, 2);
imagesc(streamflowTable_Sm_dry_noice-100);
title('(b) No Ice, Dry Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb2, customColormap);
caxis(sb2, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_noice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 3: Annual
sb3 = subplot(3, 3, 3);
imagesc(streamflowTable_Sm_ann_noice-100);
title('(c) No Ice, Annual','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb3, customColormap);
caxis(sb3, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_noice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Plotting for '2glac'
% Subplot 4: Wet season
sb4 = subplot(3, 3, 4);
imagesc(streamflowTable_Sm_wet_2ice-100);
title('(d) 2% GC, Wet Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb4, customColormap);
caxis(sb4, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_2ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 5: Dry season
sb5 = subplot(3, 3, 5);
imagesc(streamflowTable_Sm_dry_2ice-100);
title('(e) 2% GC, Dry Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb5, customColormap);
caxis(sb5, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_2ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 6: Annual
sb6 = subplot(3, 3, 6);
imagesc(streamflowTable_Sm_ann_2ice-100);
title('(f) 2% GC, Annual','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb6, customColormap);
caxis(sb6, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_2ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Plotting for '25glac'
% Subplot 7: Wet season
sb7 = subplot(3, 3, 7);
imagesc(streamflowTable_Sm_wet_25ice-100);
title('(g) 13% GC, Wet Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb7, customColormap);
caxis(sb7, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_25ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 8: Dry season
sb8 = subplot(3, 3, 8);
imagesc(streamflowTable_Sm_dry_25ice-100);
title('(h) 13% GC, Dry Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb8, customColormap);
caxis(sb8, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_25ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 9: Annual
sb9 = subplot(3, 3, 9);
imagesc(streamflowTable_Sm_ann_25ice-100);
title('(i) 13% GC, Annual','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colormap(sb9, customColormap);
caxis(sb9, caxis_values);
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_ann_25ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Set common properties for all subplots
for i = 1:9
    subplot(3, 3, i);
    set(gca, 'XTick', 1:numPrecipitationChanges, 'XTickLabel', precipitationChangesD);
    set(gca, 'YTick', 1:numTemperatureChanges, 'YTickLabel', temperatureChanges);
    xlabel('Δ Precip (%)');
    ylabel('Δ Temp (°C)');
end

figname ='Scenario_GridNumberChange_WetDryAnn_subplot';
saveas(gcf, strcat(figdir, figname, '.pdf'));
saveas(gcf, strcat(figdir, figname, '.png'));
savefig(gcf, strcat(figdir, figname));