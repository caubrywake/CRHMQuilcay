%% Assesing streamflow change and GW change for the likely scenarios
% edited by Caroline Aubry-Wake, 2023-12-12
% contact at caroline.aubrywake@gmail.com

% This script analyses the more likely scenario of th range of
% precipitation and temperature senstivity. 

% this script generates figure F5_Scenario_TimeSeries_Flowcomponent
%% Set-up
close all
clear all
cd 'F:\11_CRHM_cuchi\' % set to working folder
addpath 'F:\11_CRHM_cuchi\script\' % local directory with CRHM model results
figdir ='F:\11_CRHM_cuchi\fig\scenario\'

hru_area = [1.528 1.322 1.237 0.839 7.366 4.056 3.492 3.663 3.201 4.633 9.491 2.721 2.418 2.008 6.249 10.62 1.903 0.14 0.45];
ratio_area = hru_area./sum(hru_area);

% % glacier per scenario
hruarea = [1.528 1.322 1.237 0.839 7.366 4.056 3.492 3.663 3.201 4.633 9.491 2.721 2.418 2.008 6.249 10.62 ...
1.903 0.14 0.45];
sumarea = sum(hruarea)
hru_glaciercover = [sum(hruarea([2,3,5,6,9])./sumarea) 0, hruarea(2)./sumarea hruarea(5)./sumarea sum(hruarea([2,5])./sumarea)]

%% import the current ones
load('CRHM\output\Cuchi_20230823.mat','SWEmelt', 'icemelt','firnmelt','hru_actet', 'infil', 'runoff','basinflow','basingw', 'time')
snowmelt_cur = sum(SWEmelt.*ratio_area,2)./24;
icemelt_cur  = sum(icemelt.*ratio_area,2)./24 + sum(firnmelt.*ratio_area,2)./24;
et_cur       = sum(hru_actet.*ratio_area, 2);
infil    = sum(infil.*ratio_area, 2);
runoff   = sum(runoff.*ratio_area, 2);
rainfallrunoff_cur = (runoff + infil) -icemelt_cur;
rainfallrunoff_cur(rainfallrunoff_cur<0)=0;

T = timetable(time, basinflow/3600, basingw/3600);
TT = retime(T, 'monthly','mean');
b = TT.Var1;
g = TT.Var2;
t = TT.time;
s = b+g;
T = timetable(time, snowmelt_cur,icemelt_cur ,et_cur,rainfallrunoff_cur);
TT = retime(T, 'monthly','sum');
snowmelt_cur=TT.snowmelt_cur;
icemelt_cur=TT.icemelt_cur;
et_cur=TT.et_cur;
rainfallrunoff_cur = TT.rainfallrunoff_cur;
time_cur = TT.time;

clear basinflow basingw firnmelt hru_actet icemelt SWEmelt T TT time t rainfallrunoff infil runoff
% Import sceanrios the names of the files

%% Future scenarios
folderPath = 'CRHM\output\scenario\';  % Replace with the actual folder path
files = dir(fullfile(folderPath, '*.txt'));  % Replace '*.txt' with the appropriate file extension
% Extract file names and remove the first 8 letters
fileName = {files.name};
fileName = extractBetween(fileName, 9, 28);

%% Future 1
id_fut = find(contains(fileName, '2glac') & contains(fileName, 't_5') & contains(fileName, 'precip_11'));

clear B G S snowmelt firnmelt icemelt gw ssr hru_actet rainfallrunoff infil runoff
fn = files(id_fut).name; % file to import
H = importdata(strcat(folderPath, fn),' ',2); %  % import headers 
D = importdata(strcat(folderPath, fn)) ; % import data
B= D.data(:,2)./3600;
G= D.data(:,3)./3600;
S= B+G;
% swemelt 4:22, firnmelt 23:41, icemelt 42:60, hru_actet61:, infil, runoff, gw
snowmelt = sum(D.data(:,4:22).*ratio_area,2)./24;
icemelt  = sum(D.data(:,42:60).*ratio_area,2)./24 + sum(D.data(:,23:41).*ratio_area,2)./24;
et       = sum(D.data(:,61:79).*ratio_area, 2);
infil    = sum(D.data(:,80:98).*ratio_area, 2);
runoff   = sum(D.data(:,99: 117).*ratio_area, 2);
rainfallrunoff = (runoff + infil) -icemelt ;
rainfallrunoff(rainfallrunoff<0)=0;
time = D.data(:,1);
time= datetime(datevec(time+ 693960));
time = dateshift(time,'start','hour', 'nearest');

T = timetable(time, B, G, S);
TT = retime(T, 'monthly','mean');
B1 = TT.B;
G1 = TT.G;
S1 = TT.S;
T = timetable(time, snowmelt,icemelt ,et, rainfallrunoff);
TT = retime(T, 'monthly','sum');
snowmelt_fut=TT.snowmelt;
icemelt_fut=TT.icemelt;
et_fut=TT.et;
rainfallrunoff_fut = TT.rainfallrunoff;

time_fut = TT.time
clear basinflow basingw firnmelt hru_actet icemelt SWEmelt T TT time t  snowmelt firnmelt icemelt gw ssr hru_actet i D et

%% Future 2  
id_fut = find(contains(fileName, '2glac') & contains(fileName, 't_5') & contains(fileName, 'precip_12'));

fn = files(id_fut).name; % file to import
H = importdata(strcat(folderPath, fn),' ',2); %  % import headers 
D = importdata(strcat(folderPath, fn)) ; % import data
B= D.data(:,2)./3600;
G= D.data(:,3)./3600;
S= B+G;
% swemelt 4:22, firnmelt 23:41, icemelt 42:60, hru_actet61:, infil, runoff, gw
snowmelt = sum(D.data(:,4:22).*ratio_area,2)./24;
icemelt  = sum(D.data(:,42:60).*ratio_area,2)./24 + sum(D.data(:,23:41).*ratio_area,2)./24;
et       = sum(D.data(:,61:79).*ratio_area, 2);
infil    = sum(D.data(:,80:98).*ratio_area, 2);
runoff   = sum(D.data(:,99: 117).*ratio_area, 2);
rainfallrunoff = (runoff + infil) -icemelt ;
rainfallrunoff(rainfallrunoff<0)=0;
time = D.data(:,1);
time= datetime(datevec(time+ 693960));
time = dateshift(time,'start','hour', 'nearest');

T = timetable(time, B, G, S);
TT = retime(T, 'monthly','mean');
B2= TT.B;
G2 = TT.G;
S2 = TT.S;

T = timetable(time, snowmelt,icemelt ,et, rainfallrunoff);
TT = retime(T, 'monthly','sum');
snowmelt_fut2=TT.snowmelt;
icemelt_fut2=TT.icemelt;
et_fut2=TT.et;
rainfallrunoff_fut2 = TT.rainfallrunoff;
clear basinflow basingw firnmelt hru_actet icemelt SWEmelt T TT time t  snowmelt firnmelt icemelt gw ssr hru_actet i D et

% 
%% Future 3 
id_fut = find(contains(fileName, '2glac') & contains(fileName, 't_4') & contains(fileName, 'precip_11'));

fn = files(id_fut).name; % file to import
H = importdata(strcat(folderPath, fn),' ',2); %  % import headers 
D = importdata(strcat(folderPath, fn)) ; % import data
B= D.data(:,2)./3600;
G= D.data(:,3)./3600;
S= B+G;
% swemelt 4:22, firnmelt 23:41, icemelt 42:60, hru_actet61:, infil, runoff, gw
snowmelt = sum(D.data(:,4:22).*ratio_area,2)./24;
icemelt  = sum(D.data(:,42:60).*ratio_area,2)./24 + sum(D.data(:,23:41).*ratio_area,2)./24;
et       = sum(D.data(:,61:79).*ratio_area, 2);
infil    = sum(D.data(:,80:98).*ratio_area, 2);
runoff   = sum(D.data(:,99: 117).*ratio_area, 2);
rainfallrunoff = (runoff + infil) -icemelt ;
rainfallrunoff(rainfallrunoff<0)=0;
time = D.data(:,1);
time= datetime(datevec(time+ 693960));
time = dateshift(time,'start','hour', 'nearest');

T = timetable(time, B, G, S);
TT = retime(T, 'monthly','mean');
B3= TT.B;
G3 = TT.G;
S3 = TT.S;

T = timetable(time, snowmelt,icemelt ,et, rainfallrunoff);
TT = retime(T, 'monthly','sum');
snowmelt_fut3=TT.snowmelt;
icemelt_fut3=TT.icemelt;
et_fut3=TT.et;
rainfallrunoff_fut3 = TT.rainfallrunoff;
clear basinflow basingw firnmelt hru_actet icemelt SWEmelt T TT time t  snowmelt firnmelt icemelt gw ssr hru_actet i D et

%% Future 4 
id_fut = find(contains(fileName, '2glac') & contains(fileName, 't_4') & contains(fileName, 'precip_12'));

fn = files(id_fut).name; % file to import
H = importdata(strcat(folderPath, fn),' ',2); %  % import headers 
D = importdata(strcat(folderPath, fn)) ; % import data
B= D.data(:,2)./3600;
G= D.data(:,3)./3600;
S= B+G;
% swemelt 4:22, firnmelt 23:41, icemelt 42:60, hru_actet61:, infil, runoff, gw
snowmelt = sum(D.data(:,4:22).*ratio_area,2)./24;
icemelt  = sum(D.data(:,42:60).*ratio_area,2)./24 + sum(D.data(:,23:41).*ratio_area,2)./24;
et       = sum(D.data(:,61:79).*ratio_area, 2);
infil    = sum(D.data(:,80:98).*ratio_area, 2);
runoff   = sum(D.data(:,99: 117).*ratio_area, 2);
rainfallrunoff = (runoff + infil) -icemelt ;
rainfallrunoff(rainfallrunoff<0)=0;
time = D.data(:,1);
time= datetime(datevec(time+ 693960));
time = dateshift(time,'start','hour', 'nearest');

T = timetable(time, B, G, S);
TT = retime(T, 'monthly','mean');
B4= TT.B;
G4 = TT.G;
S4 = TT.S;

T = timetable(time, snowmelt,icemelt ,et, rainfallrunoff);
TT = retime(T, 'monthly','sum');
snowmelt_fut4=TT.snowmelt;
icemelt_fut4=TT.icemelt;
et_fut4=TT.et;
rainfallrunoff_fut4 = TT.rainfallrunoff;
clear basinflow basingw firnmelt hru_actet icemelt SWEmelt T TT time t  snowmelt firnmelt icemelt gw ssr hru_actet i D et 

%% Comapring cur and fut for most likely scenario
close all 
fig = figure('units','inches','outerposition',[0 0 8 7]);
c1 = [217 95 2]./255
c2 = [153 219 77]./255
c3 = [27 158 119]./255
c4 = [240 147 204]./255
ld =0.7;

subplot(3,2,1)
plot (time_cur, snowmelt_cur, 'k', 'LineWidth',ld); hold on
plot (time_cur, snowmelt_fut, 'Color', c1, 'LineWidth',ld);
plot (time_cur, snowmelt_fut2, 'Color', c1, 'LineStyle','--', 'LineWidth',ld);
plot (time_cur, snowmelt_fut3, 'Color', c3, 'LineWidth',ld);
plot (time_cur, snowmelt_fut4, 'Color', c3, 'LineStyle','--', 'LineWidth',ld);
lg = legend ('Present','+5^{\circ}C, +10% P','+5^{\circ}C, +20% P', '+4^{\circ}C, +10% P',  '+4^{\circ}C, +20% P', 'Orientation','Horizontal') 
lgLoc = lg.Position;
lg.Position = [0.1550    0.9387    0.7132    0.0297]
ylabel ('Snowmelt (mm/mth)')
xlim([datetime('26-Jun-2014') datetime('01-Mar-2020')])
text (datetime('30-Jul-2014'), 100, '(a)');

subplot(3,2,2)
plot (time_cur, icemelt_cur, 'k', 'LineWidth',ld); hold on
plot (time_cur, icemelt_fut, 'Color', c1, 'LineWidth',ld);
plot (time_cur, icemelt_fut2, 'Color', c1, 'LineStyle','--', 'LineWidth',ld);
plot (time_cur, icemelt_fut3, 'Color',  c3, 'LineWidth',ld);
plot (time_cur, icemelt_fut4, 'Color',  c3, 'LineStyle','--', 'LineWidth',ld);
ylabel ('Ice melt (mm/mth)')
xlim([datetime('26-Jun-2014') datetime('01-Mar-2020')])
text (datetime('30-Jul-2014'), 140, '(b)');

subplot(3,2,3)
plot (time_cur, rainfallrunoff_cur, 'k', 'LineWidth',ld); hold on
plot (time_cur, rainfallrunoff_fut, 'Color', c1, 'LineWidth',ld);
plot (time_cur, rainfallrunoff_fut2, 'Color', c1, 'LineStyle','--', 'LineWidth',ld);
plot (time_cur, rainfallrunoff_fut3, 'Color', c3, 'LineWidth',ld);
plot (time_cur, rainfallrunoff_fut4, 'Color', c3, 'LineStyle','--', 'LineWidth',ld);
ylabel ('Rainfall Runoff (mm/mth)')
xlim([datetime('26-Jun-2014') datetime('01-Mar-2020')])
text (datetime('30-Jul-2014'), 280, '(c)');

subplot(3,2,4)
plot (time_cur, et_cur, 'k', 'LineWidth',ld); hold on
plot (time_cur, et_fut, 'Color', c1, 'LineWidth',ld);
plot (time_cur, et_fut2, 'Color', c1, 'LineStyle','--', 'LineWidth',ld);
plot (time_cur, et_fut3, 'Color', c3, 'LineWidth',ld);
plot (time_cur, et_fut4, 'Color', c3, 'LineStyle','--', 'LineWidth',ld);
ylabel ('ET (mm/mth)')
xlim([datetime('26-Jun-2014') datetime('01-Mar-2020')])
text (datetime('30-Jul-2014'), 38, '(d)');

subplot(3,2,5)
plot (time_cur, b, 'k', 'LineWidth',ld); hold on
plot (time_cur, B1, 'Color', c1);
plot (time_cur, B2, 'Color', c1, 'LineStyle','--', 'LineWidth',ld);
plot (time_cur, B3, 'Color', c3);
plot (time_cur, B4, 'Color', c3, 'LineStyle','--', 'LineWidth',ld);
ylabel ('Streamflow (m^3 s^{-1})')
xlim([datetime('26-Jun-2014') datetime('01-Mar-2020')])
text (datetime('30-Jul-2014'), 3.7, '(e)');

subplot(3,2,6)
plot (time_cur, g, 'k', 'LineWidth',ld); hold on
plot (time_cur, G1, 'Color', c1, 'LineWidth',ld);
plot (time_cur, G2, 'Color', c1, 'LineStyle','--', 'LineWidth',ld);
plot (time_cur, G3, 'Color', c3, 'LineWidth',ld);
plot (time_cur, G4, 'Color', c3, 'LineStyle','--', 'LineWidth',ld);
ylabel ('Saturated zone flow (m^3 s^{-1})')
xlim([datetime('26-Jun-2014') datetime('01-Mar-2020')])
text (datetime('30-Jul-2014'), 0.55, '(f)');

%
figname ='Scenario_TimeSeries_Flowcomponent';
saveas (gcf, strcat( figdir, figname, '.pdf'));
saveas (gcf, strcat(figdir, figname, '.png'));
savefig(gcf, strcat(figdir, figname));


