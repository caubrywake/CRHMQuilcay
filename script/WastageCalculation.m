%% wastage vs melt in Quilcay
close all
clear all

cd 'G:\11_CRHM_cuchi\'
addpath 'G:\11_CRHM_cuchi\functions\' 

%% Import CRHM results
figdir = 'G:\11_CRHM_cuchi\fig\analysis\'
savedir = figdir


%% Contribution of wastage to basin yield
% basin yielf is the volumetric (m3) water generated in the bains:
% snowmelt, icemelt, firnmelt, rainfall runoff

% calculate the amount of icemelt + snowmelt that cmes form glacier
% compared to the amount of water generated in the whole baisn. and if we
% for water generated in glXIER HRU: WASTAGE AND  SEASONAL MELT
%for whole baisn: just seasonal melt (snowmelt+ rainfall runoff)


load('G:\11_CRHM_cuchi\CRHM\output\v9\Cuchi_20230823.mat','SWE', 'snowmelt_int', 'net_rain_org','SWEmelt', 'basinflow', 'basingw', 'firnmelt', 'icemelt', 'hru_actet', 'time', 'hru_snow','hru_rain')


rainfallrunoff = net_rain_org;
swemelt = snowmelt_int;% hourloy values
icemelt = (icemelt+firnmelt); % hourloy values


% retime stuff
T = timetable(time, icemelt);
TT = retime(T, 'daily', 'mean');
td = TT.time;
icemelt = TT.icemelt;


T = timetable(time, rainfallrunoff, hru_snow, hru_rain, swemelt);
TT = retime(T, 'daily', 'sum');
td = TT.time;
rainfallrunoff = TT.rainfallrunoff;
hru_snow = TT.hru_snow;
hru_rain=TT.hru_rain;
swemelt = TT.swemelt;

%% diagnosis figures
figure;
hru = [5,6,7,8]
for i = 1:4
    idx = hru(i);
subplot(2,2,i)
hold on
plot(td, cumsum(swemelt(:,idx)));
plot(td, cumsum(icemelt(:,idx)));
plot(td, cumsum(hru_snow(:,idx)));
plot(td, cumsum(hru_rain(:,idx)));
plot(td, cumsum(rainfallrunoff(:,idx)));hold on
legend ('swemelt','icemelt', 'snowfall', 'rainfall', 'rr');
title(num2str(idx))
end 


close all
figure;
hru = [16,17,18,19];
for i = 1:4
    idx = hru(i);
subplot(2,2,i)

%plot(td, swemelt(:,idx)); hold on
%plot(td, hru_snow(:,idx));
%plot(td, hru_rain(:,idx)); hold on
%plot(td, rainfallrunoff(:,idx))
%plot(time, SWE(:, idx))
plot(td, cumsum(hru_snow(:, idx))); hold on
plot(td, cumsum(swemelt(:, idx)))
legend ( 'snowfall','snowmelt');

%legend ('swemelt','snowfall', 'rainfall', 'rr');

%legend ( 'rainfall', 'rr', 'swe', 'snowfall','snowmelt');
title(num2str(idx))
end 

close all
figure;
hru = [7,8,11,14];
for i = 1:4
    idx = hru(i);
subplot(2,2,i)
plot(td, rainfallrunoff(:,idx)); hold on
plot(td, swemelt(:,idx));
plot(td, hru_snow(:,idx));
plot(td, hru_rain(:,idx));

legend ('rr','swemelt','snowfall', 'rainfall');
title(num2str(idx))
end 


clear net_rain_org SWEmelt firnmelt
%% everything in daily values
hru_area = [1.528 1.322 1.237 0.839 7.366 4.056 3.492 3.663 3.201 4.633 9.491 2.721 2.418 2.008 6.249 10.62 1.903 0.14 0.45];
hru_glacier = [2,3,5,6,9];

icemelt_glacier = calculate_melt_hru_m3hr(icemelt, hru_area, hru_glacier);
swemelt_glacier = calculate_melt_hru_m3hr(swemelt, hru_area, hru_glacier);
rainfallrunoff_glacier = calculate_melt_hru_m3hr(rainfallrunoff, hru_area, hru_glacier);
snowfall_glacier = calculate_melt_hru_m3hr(hru_snow, hru_area, hru_glacier);
wastage_glacier = (icemelt_glacier + swemelt_glacier) - snowfall_glacier
melt_glacier = snowfall_glacier-swemelt_glacier

total_abl = icemelt_glacier + swemelt_glacier;
total_acc = snowfall_glacier;
imbalanced = total_abl - total_acc;
balanced = total_acc;


swemelt_nonglacier= calculate_melt_hru_m3hr(swemelt, hru_area, [1,4, 7:8,10:19]);
rainfallrunoff_nonglacier = calculate_melt_hru_m3hr(rainfallrunoff, hru_area, [1,4, 7:8,10:19]);
snowfall_nonglacier = calculate_melt_hru_m3hr(hru_snow, hru_area,[1,4, 7:8,10:19]);
melt_nonglacier = snowfall_nonglacier-swemelt_nonglacier;
basinyield = wastage_glacier + melt_glacier + rainfallrunoff_glacier + rainfallrunoff_nonglacier + melt_nonglacier;

%% more diagnosis figures
figure
plot(td, cumsum(wastage_glacier)); hold on
plot(td, cumsum(melt_glacier)); hold on
plot(td, cumsum(snowfall_glacier))
plot(td, cumsum(swemelt_glacier))
plot(td, cumsum(icemelt_glacier))
legend ('wast', 'melt', 'sf','swemelt', 'icemelt')

figure
plot(td, cumsum(imbalanced)); hold on
plot(td, cumsum(balanced)); hold on
plot(td, cumsum(total_acc))
plot(td, cumsum(total_abl))
legend ('imbalanced','balanced',  'acc','abl')



figure
plot(td, cumsum(wastage_glacier)); hold on
plot(td, cumsum(melt_glacier)); hold on
plot(td, cumsum(snowfall_glacier))
plot(td, cumsum(swemelt_glacier))
plot(td, cumsum(icemelt_glacier))
legend ('non-balanced melt', 'melt', 'sf','swemelt', 'icemelt')

figure
plot(td, cumsum(snowfall_nonglacier)); hold on
plot(td, cumsum(swemelt_nonglacier))
plot(td, cumsum(rainfallrunoff_nonglacier)); hold on
plot(td, cumsum(melt_nonglacier)); hold on
legend('snowfall ng', 'swemelt ng', 'rr ng', 'melt nonglac')

figure;
plot(td, cumsum(imbalanced)); hold on
plot(td, cumsum(balanced)); hold on
plot(td, cumsum(swemelt_nonglacier))
plot(td, cumsum(rainfallrunoff_nonglacier)); hold on
legend ('imbal','bal','swemelt ng','rr ng')

%%
% lets dpo that annual
yr = (2015:2019)
for i = 1:length (yr)
 t1 = strcat ('01-Oct-', num2str(yr(i)-1)); 
 t2 = strcat ('30-Sep-', num2str(yr(i)));
   
a = find(td==datetime(t1));
b = find(td==datetime(t2));

x = imbalanced;
imbal_yr(i)  = sum(x(a:b));
x = balanced;
bal_yr(i) = sum(x(a:b));
x = swemelt_nonglacier;
swemeltng_yr(i)  = sum(x(a:b));
x = rainfallrunoff_nonglacier;
rrng_yr = sum(x(a:b));

x = swemelt_glacier;
swemeltg_yr(i)  = sum(x(a:b));
x = icemelt_glacier;
icemeltg_yr(i) = sum(x(a:b));
x = snowfall_glacier;
snowfallg_yr(i)  = sum(x(a:b));
x = rainfallrunoff_glacier;
rrg_yr(i) = sum(x(a:b))
end 

wastage_ratio = imbal_yr *100 ./ (rrng_yr + swemeltng_yr + bal_yr + imbal_yr + rrg_yr)
imbal_yr_mwe = imbal_yr/(sum(hru_area(hru_glacier)*10^6)) 
bal_yr_mwe =bal_yr/(sum(hru_area(hru_glacier)*10^6)) 
icemeltg_yr_mwe = icemeltg_yr./(sum(hru_area(hru_glacier)*10^6)) 
swemeltg_yr_mwe =swemeltg_yr./(sum(hru_area(hru_glacier)*10^6)) 
snowfallg_yr_mwe =snowfallg_yr./(sum(hru_area(hru_glacier)*10^6)) 


wastage_ratio_mb = imbal_yr .*100./ (bal_yr + imbal_yr);
cavy = [102 170 255]/255%pale blue
csnow= [0 102 204]/255%  mid blue
crain = [0 51 102]/255% dark blue
csnowm = [200 200 200]/255% light grey
cicem =  [140 140 140]/255% mid grey 
cfirnm = [90 90 90]/255% dark grey
%%
close all
fig = figure('units', 'inches', 'outerposition', [0 0 8 7]);

% Adjust marker face color and line width
lineWidth = 1; % Set line width
markerSize = 6; % Set marker size
subplot(3, 1, 1);
plot(yr, icemeltg_yr_mwe, '-x', 'Color', csnow, 'MarkerFaceColor', csnow, 'LineWidth', lineWidth, 'MarkerSize', markerSize); hold on
plot(yr, swemeltg_yr_mwe, '-o', 'Color', cicem, 'MarkerFaceColor', cicem, 'LineWidth', lineWidth, 'MarkerSize', markerSize);
plot(yr, snowfallg_yr_mwe, '-d', 'Color', crain, 'MarkerFaceColor', crain, 'LineWidth', lineWidth, 'MarkerSize', markerSize);
ylabel({'Annual Specific'; 'Glacier Mass Fluxes'; '(m w.e.)'})
ylim([0 2.5])
xlim([2014.8 2019.2])
text(2014.9, 2.3, '(a)')
grid on
box on
h = legend('Icemelt', 'Snowmelt', 'Snowfall', 'Orientation', 'Horizontal', 'Location', 'Northeast');


subplot(3, 1, 2);
plot(yr, imbal_yr_mwe, '-x', 'Color', csnow, 'MarkerFaceColor', csnow, 'LineWidth', lineWidth, 'MarkerSize', markerSize); hold on
plot(yr, bal_yr_mwe, '-o', 'Color', cicem, 'MarkerFaceColor', cicem, 'LineWidth', lineWidth, 'MarkerSize', markerSize);
plot(yr, imbal_yr_mwe+bal_yr_mwe, '-d', 'Color', crain, 'MarkerFaceColor', crain, 'LineWidth', lineWidth, 'MarkerSize', markerSize);
ylabel({'Annual Specific'; 'Mass Balance'; 'components (m.w.e.)'})
ylim([0 3.6])
xlim([2014.8 2019.2])
text(2014.9, 2.9, '(b)')
grid on
box on
h = legend('Imbalance melt (Abl. - Acc.)', 'Balanced melt (Acc.)', 'Total Ablation (Icemelt + Snowmelt)', 'Orientation', 'Horizontal', 'Location', 'Northeast');


% Define the bar width (e.g., 0.5 for half the default width)
barWidth = 0.6;

% Wastage ratio
subplot(3, 1, 3);
bar(yr, wastage_ratio_mb, 'FaceColor', cicem, 'EdgeColor', cicem, 'BarWidth', barWidth);
ylabel({'Wastage as percentage'; 'of glacier melt'});
ylim([0 90]);
xlim([2014.5 2019.4]);
text(2014.6, 85, '(c)');
grid on
box on

%
filename = 'WastageComp_2015_2019'
saveas (gcf,strcat(figdir, filename), 'png')
saveas (gcf,strcat(figdir, filename), 'pdf')
savefig (gcf,strcat(figdir, filename))

mean(wastage_ratio)


lw =1.1;
fig = figure('units','inches','outerposition',[0 0 8 4]);
bar(yr, wastage_ratio, 'FaceColor', cicem, 'EdgeColor', cicem, 'BarWidth', barWidth);
ylabel ({'Wastage ratio of basin yield (%)'});
ylim([0 42])
xlim([2014.5 2019.4]);
text(2014.7, 40, '(b)');
grid on
box on
filename = 'WastageRatioYield_2015_2019'
saveas (gcf,strcat(figdir, filename), 'png')
saveas (gcf,strcat(figdir, filename), 'pdf')
savefig (gcf,strcat(figdir, filename))
