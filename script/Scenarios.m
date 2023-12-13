%% Analysing the scenario runs
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
cd F:\11_CRHM_cuchi\
figdir ='fig\scenario\'
folderPath = 'CRHM\output\scenario\';  % Replace with the actual folder path
mainpath = ' F:\11_CRHM_cuchi\'
%  glacier per scenario
hruarea = [1.528 1.322 1.237 0.839 7.366 4.056 3.492 3.663 3.201 4.633 9.491 2.721 2.418 2.008 6.249 10.62 ...
1.903 0.14 0.45];
sumarea = sum(hruarea)
hru_glaciercover = [sum(hruarea([2,3,5,6,9])./sumarea) 0, hruarea(2)./sumarea hruarea(5)./sumarea sum(hruarea([2,5])./sumarea)]

%% import the current simulations
load('CRHM\output\Cuchi_20230823.mat', 'basinflow','basingw', 'icemelt','time')
T = timetable(time, basinflow/3600, basingw/3600, sum(icemelt,2));
TT = retime(T, 'daily','mean');
b = TT.Var1;
g = TT.Var2;
t = TT.time;
s = b+g;
icecur = TT.Var3;

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
snowmelt(:,i) =sum(D.data(:,4:22),2);
firnmelt(:,i) = sum(D.data(:,23:41),2);
icemelt(:,i)  = sum(D.data(:,42:60),2);
gw(:,i)       = sum(D.data(:,61:79),2);
ssr(:,i)      = sum(D.data(:,80:98),2);
et(:, i)      = sum(D.data(:,99:end), 2);
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
save ('ImportedScenarios_hourly.mat', 'time', 'B', 'G', 'S', 'b','g','s', 'snowmelt', 'firnmelt','icemelt','gw','ssr','et');

% Retime everything to daily averages
T = timetable(time, B/3600, G/3600, S/3600);
TT = retime(T, 'daily','mean');
B = TT.Var1;
G = TT.Var2;
S = TT.Var3;
timem = TT.time;

save ('ImportedScenarios_daily.mat', 'timem', 'B', 'G', 'S', 'b','g','s', 'fileName')
%% Glacier contribution to groundwater and surface flow under current conditions
% just glacier impact
close all
clear all
load('CRHM\output\scenario\ImportedScenarios_daily.mat')
figdir ='F:\11_CRHM_cuchi\fig\scenario\'

% load('E:\11_CRHM_cuchi\CRHM\prjfile\scenarios\v5\ImportedScenarios_monthly.mat')
id_noice = find(contains(fileName, 'noi') & contains(fileName, 't_0') & contains(fileName, 'precip_10'));
id_25 = find(contains(fileName, '25glac') & contains(fileName, 't_0') & contains(fileName, 'precip_10'));
id_2 = find(contains(fileName, '2glac') & contains(fileName, 't_0') & contains(fileName, 'precip_10'));
id_allice = find(contains(fileName, 'allice') & contains(fileName, 't_0') & contains(fileName, 'precip_10'));
% under current condition, hru 2 and and no ice behave the same because hru
% 2 in in the accumulation zone and so it doesnt melt

figure('Units', 'inches', 'Position', [1, 1, 7,4]);
plot(timem, s, 'k', 'linewidth', 1); hold on
plot(timem, S(:, id_noice), 'r', 'linewidth', 1);
plot(timem, g, ':k', 'linewidth', 1); hold on
plot(timem, G(:, id_noice), ':r', 'linewidth', 1);
xlim([datetime('26-Jun-2014') datetime('06-Apr-2020')])
legend ('Streamflow, 25%', 'Streamflow, 0%', 'Sat. Zone, 25%', 'Sat. Zone, 0%',  'location','northoutside', 'orientation','horizontal')
ylabel ({'Monthly-averaged';' Flow (m^3 s^{-1})'})

figname ='monthlyFlow_NoIce';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%% compile numbers
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

Sm_wet = mean(Sm(wetseason, :));
Sm_dry = mean(Sm(dryseason, :));
sm_dry = mean(sm(dryseason, :));
sm_wet = mean(sm(wetseason, :));

Bm_wet = mean(Bm(wetseason, :));
Bm_dry = mean(Bm(dryseason, :));
bm_dry = mean(bm(dryseason, :));
bm_wet = mean(sm(wetseason, :));

Gm_wet = mean(Gm(wetseason, :));
Gm_dry = mean(Gm(dryseason, :));
gm_dry = mean(gm(dryseason, :));
gm_wet = mean(gm(wetseason, :));

% For the specific indices of glacier change
Smwet_diffice = [sm_wet  Sm_wet(id_25) Sm_wet(id_noice)] 
Smdry_diffice = [sm_dry  Sm_dry(id_25) Sm_dry(id_noice)] 
Gmwet_diffice = [gm_wet  Gm_wet(id_25) Gm_wet(id_noice)] 
Gmdry_diffice = [gm_dry  Gm_dry(id_25) Gm_dry(id_noice)] 
Bmwet_diffice = [bm_wet  Bm_wet(id_25) Bm_wet(id_noice)] 
Bmdry_diffice = [bm_dry  Bm_dry(id_25) Bm_dry(id_noice)] 

Smwet_diffice_p = [sm_wet  Sm_wet(id_25) Sm_wet(id_noice)] *100./sm_wet
Smdry_diffice_p = [sm_dry  Sm_dry(id_25) Sm_dry(id_noice)]  *100./sm_dry
Gmwet_diffice_p = [gm_wet  Gm_wet(id_25) Gm_wet(id_noice)]  *100./gm_wet
Gmdry_diffice_p = [gm_dry  Gm_dry(id_25) Gm_dry(id_noice)]  *100./gm_dry
Bmwet_diffice_p = [bm_wet  Bm_wet(id_25) Bm_wet(id_noice)]  *100./bm_wet
Bmdry_diffice_p = [bm_dry  Bm_dry(id_25) Bm_dry(id_noice)]  *100./bm_dry

% Define the variable names
variableNames = {'sm_wet', 'Sm_wet(id_25)', 'Sm_wet(id_noice)';
                 'sm_dry', 'Sm_dry(id_25)', 'Sm_dry(id_noice)';
                 'gm_wet', 'Gm_wet(id_25)', 'Gm_wet(id_noice)';
                 'gm_dry', 'Gm_dry(id_25)', 'Gm_dry(id_noice)';
                 'bm_wet', 'Bm_wet(id_25)', 'Bm_wet(id_noice)';
                 'bm_dry', 'Bm_dry(id_25)', 'Bm_dry(id_noice)'};

% Create the table
data = [sm_wet, Sm_wet(id_25), Sm_wet(id_noice);
        sm_dry, Sm_dry(id_25), Sm_dry(id_noice);
        gm_wet, Gm_wet(id_25), Gm_wet(id_noice);
        gm_dry, Gm_dry(id_25), Gm_dry(id_noice);
        bm_wet, Bm_wet(id_25), Bm_wet(id_noice);
        bm_dry, Bm_dry(id_25), Bm_dry(id_noice)];

tableData = array2table(data, 'RowNames', variableNames(:, 1), 'VariableNames', { '25%', '12%', '0%'});
disp(tableData);

% Calculate percentage values
data_p =data * 100 ./ data(:, 1);
tableData_p = array2table(100-round(data_p), 'RowNames', variableNames(:, 1), 'VariableNames', { '25%', '12%', '0%'});
disp(tableData_p);


%% Different scenarios
% id for current climate, with no ice

% Lets start with juste the no ice
indicesWithNoi = find(contains(fileName, 'noi'));
indicesWithNoi_t0 = find(contains(fileName, 'noi') & contains(fileName, 't_0'));
indicesWithNoi_t2 = find(contains(fileName, 'noi') & contains(fileName, 't_2'));
indicesWithNoi_t5 = find(contains(fileName, 'noi') & contains(fileName, 't_5'));
indicesWithNoi_p80 = find(contains(fileName, 'noi') & contains(fileName, 'precip_80'));
indicesWithNoi_p10 = find(contains(fileName, 'noi') & contains(fileName, 'precip_10'));
indicesWithNoi_p12 = find(contains(fileName, 'noi') & contains(fileName, 'precip_12'));

% Fopr 25 ice
indicesWith25 = find(contains(fileName, '25'));
indicesWith25_t0 = find(contains(fileName, '25') & contains(fileName, 't_0'));
indicesWith25_t2 = find(contains(fileName, '25') & contains(fileName, 't_2'));
indicesWith25_t5 = find(contains(fileName, '25') & contains(fileName, 't_5'));
indicesWith25_p80 = find(contains(fileName, '25') & contains(fileName, 'precip_80'));
indicesWith25_p10 = find(contains(fileName, '25') & contains(fileName, 'precip_10'));
indicesWith25_p12 = find(contains(fileName, '25') & contains(fileName, 'precip_12'));

% For 2 ice
indicesWith2 = find(contains(fileName, '2glac'));
indicesWith2_t0 = find(contains(fileName, '2glac') & contains(fileName, 't_0'));
indicesWith2_t2 = find(contains(fileName, '2glac') & contains(fileName, 't_2'));
indicesWith2_t5 = find(contains(fileName, '2glac') & contains(fileName, 't_5'));
indicesWith2_p80 = find(contains(fileName, '2glac') & contains(fileName, 'precip_80'));
indicesWith2_p10 = find(contains(fileName, '2glac') & contains(fileName, 'precip_10'));
indicesWith2_p12 = find(contains(fileName, '2glac') & contains(fileName, 'precip_12'));

% For all ice
indicesWithall = find(contains(fileName, 'all'));
indicesWithall_t0 = find(contains(fileName, 'all') & contains(fileName, 't_0'));
indicesWithall_t2 = find(contains(fileName, 'all') & contains(fileName, 't_2'));
indicesWithall_t5 = find(contains(fileName, 'all') & contains(fileName, 't_5'));
indicesWithall_p80 = find(contains(fileName, 'all') & contains(fileName, 'precip_80'));
indicesWithall_p10 = find(contains(fileName, 'all') & contains(fileName, 'precip_10'));
indicesWithall_p12 = find(contains(fileName, 'all') & contains(fileName, 'precip_12'));

%% Colormaps 
% Create the colros if needed
% dark red to bright orange
numColors = 6;  % Number of colors in the sequence
yellowRGB = [253, 208, 58];
redRGB = [196, 0, 26];
colorsRGB = zeros(numColors, 3);
for i = 1:numColors
     j = (i - 1) / (numColors - 1);  % Interpolation parameter
    colorsRGB_red(i, :) = ((1 - j) * yellowRGB + j * redRGB)./255;
end

% dark blue to cyan
navyBlueRGB = [0,102, 204];
cyanRGB = [153, 255, 255];
colorsRGB = zeros(numColors-1, 3);
% Interpolate between the RGB values
for i = 1:numColors-1
    j = (i - 1) / (numColors - 1);  % Interpolation parameter
    colorsRGB_blue(i, :) = ((1 - j) * navyBlueRGB + j * cyanRGB)./255;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure with just streamflow - No ice
close all
t = timem
figure('Units', 'inches', 'Position', [1, 1, 9, 7]);
subplot(3,2,1); hold on
id = indicesWithNoi_p80
plot(t, s, ':k', 'linewidth', 2)
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1);
set(p1, {'Color'}, num2cell(colorsRGB_red, 2));
title ('No Ice, P 80%, T +0-5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 5])
legend('bl','0','1','2','3','4','5','orientation','horizontal')

subplot(3,2,2); hold on
id = indicesWithNoi_t0
id = id([4,5,1,2,3])
plot(t, s, ':k', 'linewidth', 2)
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1);
set(p1, {'Color'}, num2cell(colorsRGB_blue, 2));
title ('No Ice, P 80-120%, T 0')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 5])
legend('baseline','80','90','100','110','120', 'orientation','horizontal')

subplot(3,2,3); hold on
id = indicesWithNoi_p10
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_red, 2));
title ('No Ice, P 100%, T +0-5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 5])

subplot(3,2,4); hold on
id = indicesWithNoi_t2
id = id([4,5,1,2,3])
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_blue, 2));
title ('No Ice, P 80-120%, T +2')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 5])

subplot(3,2,5); hold on
id = indicesWithNoi_p12
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_red, 2));
title ('No Ice, P 120%, T +0-5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 5])

subplot(3,2,6); hold on
id = indicesWithNoi_t5
id = id([4,5,1,2,3])

p1= plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_blue, 2));
title ('No Ice, P 80-120%, T +5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 5])

figname ='Scenario_NoIce_6subplot';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%% %% figure with just streamflow - 25Glacier
close all

figure('Units', 'inches', 'Position', [1, 1, 9, 7]);
subplot(3,2,1); hold on
id = indicesWith25_p80
plot(t, s, ':k', 'linewidth', 2)
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1);
set(p1, {'Color'}, num2cell(colorsRGB_red, 2));
title ('13% GC, P 80%, T +0-5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])
legend('bl','0','1','2','3','4','5', 'Cur, NoIce','orientation','horizontal')

subplot(3,2,2); hold on
id = indicesWith25_t0
id = id([4,5,1,2,3])
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1);
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_blue, 2));
title ('13% GC, P 80-120%, T 0')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])
legend('baseline','80','90','100','110','120','Cur, NoIce', 'orientation','horizontal')

subplot(3,2,3); hold on
id = indicesWith25_p10
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_red, 2));
title ('13% GC, P 100%, T +0-5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])

subplot(3,2,4); hold on
id = indicesWith25_t2
id = id([4,5,1,2,3])
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_blue, 2));
title ('13% GC, P 80-120%, T +2')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])

subplot(3,2,5); hold on
id = indicesWith25_p12
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_red, 2));
title ('13% GC, P 120%, T +0-5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])

subplot(3,2,6); hold on
id = indicesWith25_t5
id = id([4,5,1,2,3])

p1= plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(timem, S(:, id_25),'-k', 'linewidth',1)
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_blue, 2));
title ('13% GC, P 80-120%, T +5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])


figname ='Scenario_25Ice_6subplot';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%% 2 Glacier
close all

figure('Units', 'inches', 'Position', [1, 1, 9, 7]);
subplot(3,2,1); hold on
id = indicesWith2_p80
plot(t, s, ':k', 'linewidth', 2)
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1);
set(p1, {'Color'}, num2cell(colorsRGB_red, 2));
title ('2% GC, P 80%, T +0-5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])
legend('bl','0','1','2','3','4','5', 'Cur, NoIce','orientation','horizontal')

subplot(3,2,2); hold on
id = indicesWith2_t0
id = id([4,5,1,2,3])
plot(t, s, ':k', 'linewidth', 2)
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1);
set(p1, {'Color'}, num2cell(colorsRGB_blue, 2));
title ('2% GC, P 80-120%, T 0')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])
legend('baseline','80','90','100','110','120','orientation','horizontal')

subplot(3,2,3); hold on
id = indicesWith2_p10
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_red, 2));
title ('2% GC, P 100%, T +0-5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])

subplot(3,2,4); hold on
id = indicesWith2_t2
id = id([4,5,1,2,3])
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_blue, 2));
title ('2% GC, P 80-120%, T +2')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])

subplot(3,2,5); hold on
id = indicesWith2_p12
p1 = plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_red, 2));
title ('2% GC, P 120%, T +0-5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])

subplot(3,2,6); hold on
id = indicesWith2_t5
id = id([4,5,1,2,3])
p1= plot(timem, S(:, id), 'b', 'linewidth', 1); hold on
plot(t, s, ':k', 'linewidth', 2)
set(p1, {'Color'}, num2cell(colorsRGB_blue, 2));
title ('2% GC, P 80-120%, T +5')
ylabel ('Monthly Flow (m^3 s^{-1})')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim ([0 6])


figname ='Scenario_2Ice_6subplot';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%% Compile numbers for streamflow change for all simulations

% make a grid
temperatureChanges = 0:1:5;   % Temperature changes from 0 to +4°C
precipitationChanges = 80:10:120;   % Precipitation changes from -20% to +10%
numTemperatureChanges = numel(temperatureChanges);
numPrecipitationChanges = numel(precipitationChanges);

streamflowTable_Sm_dry_noice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Sm_wet_noice = zeros(numTemperatureChanges, numPrecipitationChanges);

streamflowTable_Sm_dry_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Sm_wet_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);

streamflowTable_Sm_dry_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);
streamflowTable_Sm_wet_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);
%% Iterate over the filenames and populate the tables

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
        precipitationChange = precipitationChange*10
    end 

    % Find the indices for the temperature and precipitation changes in the tables
    temperatureIndex = temperatureChange + 1;
    precipitationIndex = find(precipitationChanges == precipitationChange);

    % Calculate the streamflow changes for the scenarios
    streamflowChange_Sm_dry = Sm_dry(i) - sm_dry;
    streamflowChange_Sm_wet = Sm_wet(i) - sm_wet;
    streamflowChange_Bm_dry = Bm_dry(i) - bm_dry;
    streamflowChange_Bm_wet = Bm_wet(i) - bm_wet;
    streamflowChange_Gm_dry = Gm_dry(i) - gm_dry;
    streamflowChange_Gm_wet = Gm_wet(i) - gm_wet;
    % Calculate the streamflow changes as percentages
    streamflowChange_Sm_dry_p = Sm_dry(i) * 100 / sm_dry;
    streamflowChange_Sm_wet_p = Sm_wet(i) * 100 / sm_wet;
    streamflowChange_Bm_dry_p = Bm_dry(i) * 100 / bm_dry;
    streamflowChange_Bm_wet_p = Bm_wet(i) * 100 / bm_wet;
    streamflowChange_Gm_dry_p = Gm_dry(i) * 100 / gm_dry;
    streamflowChange_Gm_wet_p = Gm_wet(i) * 100 / gm_wet;
    % Update the tables with the streamflow change values based on the glacString
    switch glacString
        case '25glac'
            streamflowTable_Sm_dry_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry_p;
            streamflowTable_Sm_wet_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet_p;
            streamflowTable_Sm_dry_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry;
            streamflowTable_Sm_wet_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet;

            streamflowTable_Bm_dry_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry_p;
            streamflowTable_Bm_wet_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet_p;
            streamflowTable_Bm_dry_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry;
            streamflowTable_Bm_wet_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet; 

            streamflowTable_Gm_dry_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry_p;
            streamflowTable_Gm_wet_25ice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet_p;
            streamflowTable_Gm_dry_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry;
            streamflowTable_Gm_wet_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet;
        case '2glac'
            streamflowTable_Sm_dry_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry_p;
            streamflowTable_Sm_wet_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet_p;
            streamflowTable_Sm_dry_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry;
            streamflowTable_Sm_wet_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet;

            streamflowTable_Bm_dry_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry_p;
            streamflowTable_Bm_wet_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet_p;
            streamflowTable_Bm_dry_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry;
            streamflowTable_Bm_wet_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet; 

            streamflowTable_Gm_dry_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry_p;
            streamflowTable_Gm_wet_2ice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet_p;
            streamflowTable_Gm_dry_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry;
            streamflowTable_Gm_wet_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet;
        case 'noice'
            streamflowTable_Sm_dry_noice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry_p;
            streamflowTable_Sm_wet_noice(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet_p;
            streamflowTable_Sm_dry_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry;
            streamflowTable_Sm_wet_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet;
            
            streamflowTable_Bm_dry_noice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry_p;
            streamflowTable_Bm_wet_noice(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet_p;
            streamflowTable_Bm_dry_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_dry;
            streamflowTable_Bm_wet_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Bm_wet; 

            streamflowTable_Gm_dry_noice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry_p;
            streamflowTable_Gm_wet_noice(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet_p;
            streamflowTable_Gm_dry_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_dry;
            streamflowTable_Gm_wet_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Gm_wet;
    end
end

%% Just get the number for streamflow, not streamflow change
% % make a grid
% temperatureChanges = 0:1:5;   % Temperature changes from 0 to +4°C
% precipitationChanges = 80:10:120;   % Precipitation changes from -20% to +10%
% numTemperatureChanges = numel(temperatureChanges);
% numPrecipitationChanges = numel(precipitationChanges);
% 
% streamflowTable_Sm_dry_noice = zeros(numTemperatureChanges, numPrecipitationChanges);
% streamflowTable_Sm_wet_noice = zeros(numTemperatureChanges, numPrecipitationChanges);
% 
% streamflowTable_Sm_dry_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);
% streamflowTable_Sm_wet_25ice = zeros(numTemperatureChanges, numPrecipitationChanges);
% 
%  streamflowTable_Sm_dry_allice = zeros(numTemperatureChanges, numPrecipitationChanges);
%  streamflowTable_Sm_wet_allice = zeros(numTemperatureChanges, numPrecipitationChanges);
% 
% streamflowTable_Sm_dry_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);
% streamflowTable_Sm_wet_2ice = zeros(numTemperatureChanges, numPrecipitationChanges);
% %%
% 
% % Iterate over the filenames and populate the tables
% for i = 1:numel(fileName)
%     filename = fileName{i};
%     
%     % Check if the filename contains '25glac'
%     if contains(filename, '25glac')
%         glacString = '25glac';
%     % Check if the filename contains '2glac'
%     elseif contains(filename, '2glac')
%         glacString = '2glac';
%     % Check if the filename contains 'allglac'
%     elseif contains(filename, 'allice')
%         glacString = 'allice';
%     % Check if the filename contains 'noice'
%     elseif contains(filename, 'noice')
%         glacString = 'noice';
%     end
%        
%     % Extract the temperature change and precipitation change values from the filename
%     temperatureChange = str2double(regexp(filename, '(?<=_t_)\d+', 'match'));
%     precipitationChange = str2double(regexp(filename, '(?<=precip_)\d+', 'match'));
% 
%     if precipitationChange < 70 
%         precipitationChange = precipitationChange*10
%     end 
% 
%     % Find the indices for the temperature and precipitation changes in the tables
%     temperatureIndex = temperatureChange + 1;
%     precipitationIndex = find(precipitationChanges == precipitationChange);
% 
%     % Calculate the streamflow changes for the scenarios
%     streamflowChange_Sm_dry = Sm_dry(i);
%     streamflowChange_Sm_wet = Sm_wet(i);
% 
%     % Update the tables with the streamflow change values based on the glacString
%     switch glacString
%         case '25glac'
%             streamflowTable_Sm_dry_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry;
%             streamflowTable_Sm_wet_25ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet;
%         case '2glac'
%             streamflowTable_Sm_dry_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry;
%             streamflowTable_Sm_wet_2ice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet;
%         case 'allice'
%             streamflowTable_Sm_dry_allice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry;
%             streamflowTable_Sm_wet_allice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet;
%         case 'noice'
%             streamflowTable_Sm_dry_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_dry;
%             streamflowTable_Sm_wet_noice_n(temperatureIndex, precipitationIndex) = streamflowChange_Sm_wet;
%     end
% end

%% V2
%% Create subplots for each table
close all
figure('Units', 'inches', 'Position', [1, 1, 9, 7]);

darkblue = [10 46 146]./255
lightblue = [80 110 193]./255
lightred = [191/255 84/255 84/255]
darkred  =[107/255, 0, 0]
% Define the colormap
% Define custom colormap points

% 25 ice
% Subplot 7: Bm_dry - bm_dry
sb7 = subplot(3, 2, 5);
imagesc(streamflowTable_Sm_dry_25ice-100);
title('(e) 13% GC, Dry Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
% Add text labels% Define custom colormap points
colorPoints = [darkred; 1, 1, 1; darkblue]; % Blue-White-Red
centerPosition = 0.3;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb7, customColormap);% Apply the custom colormap to the plot
caxis(sb7, [20 280]-100)

for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_25ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 8: Bm_wet - bm_wet
sb8 = subplot(3, 2, 6);
imagesc(streamflowTable_Sm_wet_25ice-100);
title('(f) 13% GC, Wet season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
% Define custom colormap points
colorPoints = [lightred; 1, 1, 1; lightblue]; % Blue-White-Red
centerPosition = 0.46;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb8, customColormap);% Apply the custom colormap to the plot
caxis(sb8, [30 180]-100)% Add text labels

for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_25ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% 2%      
%Subplot 5: Gm_dry - gm_dry
sb5 = subplot(3, 2, 3);
imagesc(streamflowTable_Sm_dry_2ice-100);
title(' (c) 2% GC, Dry seaosn','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colorPoints = [darkred; 1, 1, 1; 1,1,1]; % Blue-White-Red
centerPosition = 0.99;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb5, customColormap);
caxis(sb5, [0 110]-100)
% Add text labels
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_2ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 6: Gm_wet - gm_wet
sb6 = subplot(3, 2, 4);
imagesc(streamflowTable_Sm_wet_2ice-100);
title('(d) 2% GC, Wet season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');colorPoints = [lightred; 1, 1, 1; 1,1,1]; % Blue-White-Red
centerPosition = 0.99;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb6, customColormap);
caxis(sb6, [30 100]-100)
% Add text labels
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_2ice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end


% Subplot 1: Gm_dry - gm_dry
sb1 = subplot(3, 2, 1);
imagesc(streamflowTable_Sm_dry_noice-100);
title('(a) 0% GC, Dry season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');colorPoints = [darkred; 1, 1, 1; 1,1 ,1]; % Blue-White-Red
centerPosition = 0.99;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb1, customColormap);
caxis(sb1, [0 100]-100)
% Add text labels
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_noice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end


% Subplot 2: Gm_wet - gm_wet
sb2 = subplot(3, 2, 2);
imagesc(streamflowTable_Sm_wet_noice-100);
title('(b) 0% GC, Wet season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colorPoints = [lightred; 1, 1, 1; 1, 1, 1]; % Blue-White-Red
centerPosition = 0.99;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb2, customColormap);
caxis(sb2, [30 100]-100)
% Add text labels
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_noice(i, j))-100), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end


% Set common properties for all subplots
 precipitationChangesD = [-20 -10 0 10 20]
for i = 1:6
    subplot(3, 2, i);
    set(gca, 'XTick', 1:numPrecipitationChanges, 'XTickLabel', precipitationChangesD);
    set(gca, 'YTick', 1:numTemperatureChanges, 'YTickLabel', temperatureChanges);
    xlabel('Δ Precip (%)');
    ylabel('Δ Temp (°C)');

end

%
figname ='Scenario_GridPercentChange_WetDry_6subplot';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))


%% Create subplots for each table with number at absolute change and color ar % change
close all
figure('Units', 'inches', 'Position', [1, 1, 9, 7]);

darkblue = [10 46 146]./255
lightblue = [80 110 193]./255
lightred = [191/255 84/255 84/255]
darkred  =[107/255, 0, 0]
% Define the colormap
% Define custom colormap points

% 25 ice
% Subplot 7: Bm_dry - bm_dry
sb7 = subplot(3, 2, 5);
imagesc(streamflowTable_Sm_dry_25ice-100);
title('(e) 13% GC, Dry Season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
% Add text labels% Define custom colormap points
colorPoints = [darkred; 1, 1, 1; darkblue]; % Blue-White-Red
centerPosition = 0.3;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb7, customColormap);% Apply the custom colormap to the plot
caxis(sb7, [20 280]-100)

for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_25ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 8: Bm_wet - bm_wet
sb8 = subplot(3, 2, 6);
imagesc(streamflowTable_Sm_wet_25ice-100);
title('(f) 13% GC, Wet season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
% Define custom colormap points
colorPoints = [lightred; 1, 1, 1; lightblue]; % Blue-White-Red
centerPosition = 0.46;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb8, customColormap);% Apply the custom colormap to the plot
caxis(sb8, [30 180]-100)% Add text labels

for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_25ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% 2%      
%Subplot 5: Gm_dry - gm_dry
sb5 = subplot(3, 2, 3);
imagesc(streamflowTable_Sm_dry_2ice-100);
title(' (c) 2% GC, Dry season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colorPoints = [darkred; 1, 1, 1; 1,1,1]; % Blue-White-Red
centerPosition = 0.99;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb5, customColormap);
caxis(sb5, [0 110]-100)
% Add text labels
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_2ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Subplot 6: Gm_wet - gm_wet
sb6 = subplot(3, 2, 4);
imagesc(streamflowTable_Sm_wet_2ice-100);
title('(d) 2% GC, Wet season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');colorPoints = [lightred; 1, 1, 1; 1,1,1]; % Blue-White-Red
centerPosition = 0.99;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb6, customColormap);
caxis(sb6, [30 100]-100)
% Add text labels
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_2ice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end


% Subplot 1: Gm_dry - gm_dry
sb1 = subplot(3, 2, 1);
imagesc(streamflowTable_Sm_dry_noice-100);
title('(a) 0% GC, Dry season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');colorPoints = [darkred; 1, 1, 1; 1,1 ,1]; % Blue-White-Red
centerPosition = 0.99;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb1, customColormap);
caxis(sb1, [0 100]-100)
% Add text labels
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_dry_noice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end


% Subplot 2: Gm_wet - gm_wet
sb2 = subplot(3, 2, 2);
imagesc(streamflowTable_Sm_wet_noice-100);
title('(b) 0% GC, Wet season','HorizontalAlignment', 'right', 'FontWeight', 'normal');
cbar = colorbar;
ylabel(cbar, '% Change from baseline');
colorPoints = [lightred; 1, 1, 1; 1, 1, 1]; % Blue-White-Red
centerPosition = 0.99;% Define the position of the center color (white)
x = linspace(0, 1, 256);% Create a vector of color indices
newPositions = [0, centerPosition, 1];% Calculate the new color positions
customColormap = interp1(newPositions, colorPoints, x);% Interpolate colormap using interp1
colormap(sb2, customColormap);
caxis(sb2, [30 100]-100)
% Add text labels
for i = 1:numTemperatureChanges
    for j = 1:numPrecipitationChanges
        text(j, i, num2str(round(streamflowTable_Sm_wet_noice_n(i, j),2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 8);
    end
end

% Set common properties for all subplots
 precipitationChangesD = [-20 -10 0 10 20]
for i = 1:6
    subplot(3, 2, i);
    set(gca, 'XTick', 1:numPrecipitationChanges, 'XTickLabel', precipitationChangesD);
    set(gca, 'YTick', 1:numTemperatureChanges, 'YTickLabel', temperatureChanges);
    xlabel('Δ Precip (%)');
    ylabel('Δ Temp (°C)');

end

figname ='Scenario_GridNumberChange_WetDry_6subplot';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
