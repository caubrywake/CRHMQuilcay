% CRHM model evaluation and streamflow component analysis
% edited by Caroline Aubry-Wake, 2023-12-12
% contact at caroline.aubrywake@gmail.com

% This script provides streamflow and icemelt model evaluation by comparing
% CRHM outputs to measurements. It also analyzed streamflow components,
% both montlhy and seasonally. e% This scripts also generates F1 in the
% manuscript in addition to figures to support the analysis.
% 
% List of tables generated and exported to fig directory (as .csv): 
%  	Stats_ModelStreamflow.csv
% 	Stats_ModelStreamflow_weekly.csv
% 	MonthlyFluxes_IceNoIce.csv
% 	WetDryFluxes_IceNoIce.csv
% 	Leakage_IceContribution.csv
% 	WetDryFluxes.csv
% 	WetDryFluxes_Ratios.csv
% 	WetDryFluxes_Ratios_Yr.csv
% 	ET_MonthlyFluxes_HRU8_14_mmday.csv
% 	ET_WetDryFluxes_HRU8_14.csv
% 	AnnualIceMelt_perHRU.csv
% 	IceMelt_perHRU_23Jun2014_10Jul2014.csv
% 
% List of figures generated and saved (as .pdf and .png)
%   Basinflow_scatterplot
% 	Basinflow_SufraceSubsurface_CDA
% 	RunoffComponent_m3s_comparedtostreamflow
% 	F3_Basinflow_SurfaceSubsurface_CDA_withflowcomponent
% 	DailyET
% 	SWEexample
%% Set-up
close all
clear all

cd 'F:\11_CRHM_cuchi\'
figdir = 'F:\11_CRHM_cuchi\fig\modeleval\'; % where figures are saved
addpath 'F:\11_CRHM_cuchi\functions\' % local directory with CRHM model results

%% Import measured streamflow
d = readtable('data\Q\cda_lev-q-t.csv'); % Read the table from the file
d = table2timetable(d);
d = retime(d,'daily', 'mean');
Q = d.Discharge_cms_;
Qt = d.Datetime;

% cut to 2013-2020
a = find(Qt == '26-June-2014');
Qt(1:a)=[];
Q(1:a)=[];

% set anomaly period to nan
a = find(Qt == '04-Aug-2016');
b = find(Qt == '07-Dec-2016');
Q(a:b)=nan;

%% Load modelled streamflow
load('CRHM\output\Cuchi_20230823.mat', 'basinflow','basingw', 'time')
a = find(time =='26-June-2014');

% cut the first year (spin up)
basingw(1:a)=nan;
basinflow(1:a)=nan;

%% Set measured and modleled on same time step
t = timetable(time-days(2), basinflow/3600, basingw/3600); % streamflow, m3/s (2 days as CRHM round forward and matlabd round backward in time)
tt = retime(t, Qt, 'linear'); % retime modelled streamfloe to emasured flow (daily)
bf_d = tt.Var1;
gw_d = tt.Var2;

t = timetable(Qt, Q, bf_d, gw_d); % retime modelled and measured flow to weekly timestep
tt = retime(t, 'weekly', 'mean');
bf_w = tt.bf_d;
gw_w = tt.gw_d;
Qw = tt.Q;
Qtw = tt.Qt;

a = find(isnan(Q));
bf_d(a)=nan;
gw_d(a)=nan;

figure
plot(Qt, smooth(bf_d+gw_d, 10)); hold on
plot(Qt, Q)
plot(Qt, gw_d)
xlabel ('Time')
ylabel ('Streamflow (m3/s)')
legend ('streamflow - modelled','streamflow - measured', 'gw - modelled')
ylim([0 8])
%% Calculate Performance Metrics
bf_mod = bf_d+gw_d; % groundwater and surface water
bf_meas = Q; % measuref flow
t = Qt;
% remove nan
a = find(isnan(bf_mod));
bf_mod(a)=[];
bf_meas(a)=[];
t(a)=[];

n = numel(bf_mod);
% for the 3 years
for i = 1:1
RMSE(i) =  sqrt(mean((bf_mod(:, i) - bf_meas(:, 1)).^2));
r = corrcoef(bf_mod(:, i),bf_meas(:, 1)); R2(i) = r(2)^2;
MAE(i) = mean(abs(bf_mod(:, i) - bf_meas(:, 1)));
obs = [datenum(t) bf_meas(:, 1)];
mod = [datenum(t) bf_mod(:, i)];
[NSE(i), metricid] = nashsutcliffe(obs, mod);
MB(i) = (sum(mod(:, 2)-obs(:, 2))/sum(obs(:, 2)));
[kge(i),r, relvar,bias] = klinggupta(mod(:,2), obs(:,2));
end

varname = {'All years'};
NSEall = [NSE]; 
RMSEall = [RMSE];
MAEall = [MAE];
R2all = [R2];
KGEall = [kge];
stat_crhm = table(varname, NSEall,KGEall, RMSEall, MAEall, R2all, MB)
writetable(stat_crhm, strcat(figdir, 'Stats_ModelStreamflow.csv'))

% Retime to weekly
figure
plot(Qtw, Qw); hold on
plot(Qtw, bf_w+gw_w)
plot(Qtw, gw_w)
legend ('meas','mod')
bf_mod = bf_w+gw_w;
bf_meas = Qw;
t = Qtw;
a = find(isnan(bf_meas));
bf_mod(a)=[];
bf_meas(a)=[];
t(a)=[];
n = numel(bf_mod);
% for the 3 years
for i = 1:1
RMSE(i) =  sqrt(mean((bf_mod(:, i) - bf_meas(:, 1)).^2));
r = corrcoef(bf_mod(:, i),bf_meas(:, 1)); R2(i) = r(2)^2;
MAE(i) = mean(abs(bf_mod(:, i) - bf_meas(:, 1)));
obs = [datenum(t) bf_meas(:, 1)];
mod = [datenum(t) bf_mod(:, i)];
[NSE(i), metricid] = nashsutcliffe(obs, mod);
MB(i) = (sum(mod(:, 2)-obs(:, 2))/sum(obs(:, 2)));
[ kge(i),r, relvar,bias] = klinggupta(mod(:,2), obs(:,2));
end
varname = {'All years'};
NSEall = [NSE]; 
RMSEall = [RMSE];
MAEall = [MAE];
R2all = [R2];
KGEall = [kge];
stat_crhm_week = table(varname, NSEall, KGEall, RMSEall, MAEall, R2all, MB)
writetable(stat_crhm_week, strcat(figdir, 'Stats_ModelStreamflow_weekly.csv'))

%% Figures

%% Scatterplot Mod-Meas
figure; 
sc1 = scatter(bf_meas, bf_mod, '.k'); hold on
 
ls = lsline;
rf = refline(1,0);
ls(1).Color = 'r';
rf(1).Color = 'k';
corrcoef(bf_meas, bf_mod, 'rows', 'pairwise');
xlabel ('Daily Mean Streamflow, Measured (m^3/s)');
ylabel('Daily Mean Streamflow, Modelled (m^3/s)');
legend('obs vs mod', 'least square line','1:1 line', 'location', 'Southeast') ;
grid on
text (.2,5, {strcat('R2=',num2str(R2));strcat('RMSE=',num2str(RMSE));strcat('NSE=',num2str(NSE))});

figname ='Basinflow_scatterplot';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
%% SW-GW contributions
lw =1
fig = figure('units','inches','outerposition',[0 0 8 5]);

a = area(Qtw, [gw_w,bf_w]); hold on
a(1).FaceColor = [0.7 0.7 0.7];
a(2).FaceColor = [102 178 255]./255;
a(1).EdgeColor = [0.7 0.7 0.7];
a(2).EdgeColor = [102 178 255]./255;
plot(Qtw , Qw, ':k', 'Linewidth', lw); hold on
legend ('Modelled SW', 'Modelled GW', 'Measured Streamflow');
ylabel ({'Weekly Mean' ;'streamflow, m^3/s'});
ylim([0 6])
%xlim ([datetime('01-May-2013') datetime('15-Oct-2013')])

figname ='Basinflow_SufraceSubsurface_CDA';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
%%%%%%%

%% Streamflow components
load('CRHM\output\Cuchi_20230823.mat', 'SWEmelt', 'basinflow', 'basingw', 'firnmelt', 'icemelt', 'infil', 'runoff','infil','hru_actet', 'time', 'hru_snow','hru_rain')
rainfallrunoff = (runoff + infil) -icemelt/24 ;
rainfallrunoff(rainfallrunoff<0)=0;
swemelt = SWEmelt/24;
icemelt = icemelt+firnmelt;
hru_area = [1.528 1.322 1.237 0.839 7.366 4.056 3.492 3.663 3.201 4.633 9.491 2.721 2.418 2.008 6.249 10.62 1.903 0.14 0.45];
ratio_area = hru_area./sum(hru_area);

% converting each from mm/hr to m3/sec 
a = swemelt* 0.001 /3600 ; % converting mm to m
b = a.*(hru_area*1000000);% converting to m3 - this
c = sum(b,2); % this gives m3/hr
swemelt_m3 = c;


a = icemelt/24 * 0.001 /3600; % converting mm to m
b = a.*(hru_area*1000000);% converting to m3 - this
c = sum(b,2); % this gives m3/hr
icemelt_m3 = c;

a = rainfallrunoff* 0.001 /3600; % converting mm to m
b = a.*(hru_area*1000000);% converting to m3 - this
c = sum(b,2); % this gives m3/hr
rr_m3 = c;

T = timetable(time, rr_m3, swemelt_m3, icemelt_m3);
TT = retime(T, 'daily', 'mean');
td = TT.time;
rr_m3 = TT.rr_m3;
swemelt_m3=TT.swemelt_m3;
icemelt_m3=TT.icemelt_m3;

figure
area(td, [ smooth(icemelt_m3, 5),  smooth(swemelt_m3,5), smooth(rr_m3,5)], 'Linestyle', 'none'); hold on
plot(Qt, Q, 'k', 'linewidth', 0.8)
legend ('rainfall runoff', 'snowmelt runoff', 'icemelt runoff', 'Measured streamflow')
ylabel ('Daily-Averaged Runoff Components (m^3 s^{-1})')
ylim([0 9])
xlim([datetime('26-Jun-2014') datetime('06-Apr-2020')])
hold on
figname ='RunoffComponent_m3s_comparedtostreamflow';
saveas (gcf, strcat( figdir, figname, '.pdf'));
saveas (gcf, strcat(figdir, figname, '.png'));
savefig(gcf, strcat(figdir, figname));

%% Import the no ice simualtions and preprocess 
load('output\Cuchi_NoIce_20230823.mat', 'basinflow_noice','basingw_noice', 'time')
a = find(time =='26-June-2014');
% cut the first year 9spin upP
basingw(1:a)=nan;
basinflow(1:a)=nan;

t = timetable(time, basinflow_noice/3600, basingw_noice/3600);
tt = retime(t, 'weekly', 'mean');
bf_w_noice = tt.Var1;
gw_w_noice = tt.Var2;
figure
plot(td, bf_w_noice); hold on
plot(td, bf_d)

%% Import the no groundwater simulation
load('CRHM\output\Cuchi_NoGW_20230823.mat')
time = tt_nogw;
a = find(time =='26-June-2014');
% cut the first year 9spin upP
bf_nogw(1:a)=nan;
gw_nogw(1:a)=nan;


T = timetable(time, bf_nogw/3600, gw_nogw/3600);
TT= retime(T, 'weekly', 'mean');
bf_w_nogw = TT.Var1;
gw_w_nogw = TT.Var2;
tt = TT.time;

figure
plot(tt, bf_w_nogw); hold on
plot(tt, gw_w_nogw)
%% lot streamfoloe components
lw =1.1;
fig = figure('units','inches','outerposition',[0 0 8 7]);
subplot(3,1,1)
a = area(td, [gw_d,bf_d]); hold on
a(1).FaceColor = [0.7 0.7 0.7];
a(2).FaceColor = [102 178 255]./255;
a(1).EdgeColor = [0.7 0.7 0.7];
a(2).EdgeColor = [102 178 255]./255;
plot(Qtw , Qw, ':k', 'Linewidth', lw); hold on
legend ( 'GW','Ov. + Vad.', 'Measured Streamflow', 'orientation','horizontal');
ylabel ({'Weekly Mean' ;'Streamflow, m^3/s'});
ylim([0 6]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
text (datetime('10-Jul-2014'), 5.6, '(a)');
% 
subplot(3,1,2)
a = area(td, [gw_w_noice,bf_w_noice]); hold on
a(1).FaceColor = [0.45 0.45 0.45];
a(2).FaceColor = [0 76 153]./255;
a(1).EdgeColor = [0.45 0.45 0.45];
a(2).EdgeColor = [0 76 153]./255;
plot(td, gw_w_nogw, ':k','Linewidth', lw,'color', [0.7 0.7 0.7])
plot(td, bf_w_nogw+gw_w_nogw, ':k', 'Linewidth', lw, 'color', [102 178 255]./255); hold on
plot(td, gw_d,'-k', 'Linewidth', lw,'color', [0.7 0.7 0.7]);
plot(td, bf_d+gw_d,'-k', 'Linewidth', lw, 'color', [102 178 255]./255);
legend('GW - No Ice', 'Ov. + Vad. - No Ice', 'GW - low GW', 'Ov. + Vad. - low GW','GW','Ov. + Vad.', 'orientation','horizontal');
ylabel ({'Weekly Mean' ;'Streamflow, m^3/s'});
ylim([0 6]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
text (datetime('10-Jul-2014'), 5.6, '(b)');


subplot(3,1,3)
a = area(t_week, [icebasin, swebasin, rrbasin],  'Linestyle', 'none'); hold on
b = area(t_week, -etbasin,  'Linestyle', 'none'); 
legend ('icemelt runoff','snowmelt runoff','rainfall runoff',  'ET' ,'orientation','horizontal');
ylabel ({'Weekly Basin-Averaged';'Hydrological Components (mm)'});
a(3).FaceColor = [0 76 153]./255;
a(2).FaceColor = [255 128 0]./255;
a(1).FaceColor = [153 153 255]./255;
b(1).FaceColor = [41 135 0]./255;
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim([-10 80]);
text (datetime('10-Jul-2014'), 73, '(c)');

figname ='F3_Basinflow_SurfaceSubsurface_CDA_withflowcomponent';
saveas (gcf, strcat( figdir, figname, '.pdf'));
saveas (gcf, strcat(figdir, figname, '.png'));
savefig(gcf, strcat(figdir, figname));

%% Calculate with Ice vs no Ice ratios
T = timetable(time, basinflow/3600, basingw/3600,  basinflow_noice/3600, basingw_noice/3600); % streamflow, m3/s
TT = retime(T, 'monthly', 'sum');
 % all:

months = month(TT.time);
monthNames = datestr(TT.time, 'mmm');
% Find the indices of January months
clear averageMth
for i = 1:12
mthIndices = find(months == i);
% Extract the monthly data
Data = TT(mthIndices, :);
% Calculate the average of data
averageMth (i,:) = varfun(@mean, Data(:, 1:4));
end 

% monthly values
averageMth = table2array(averageMth);
averageMth = [averageMth; sum(averageMth)]
mth = {'January'; 'February'; 'March'; 'April'; 'May'; 'June'; 'July'; 'August'; 'September'; 'October'; 'November'; 'December';' Sum'};
T = table(mth, averageMth(:,1),averageMth(:,2), averageMth(:,3), averageMth(:,4));
T.Properties.VariableNames = {'Month','SurfaceWater - Ice ','GW - Ice','SurfaceWater - No Ice ', 'GW - No Ice'}
writetable(T, strcat(figdir, 'MonthlyFluxes_IceNoIce.csv'))

% for wet and dry season insetad
wetseason = [10,11,12,1,2,3,4];
dryseason= [5:9];
peakdry = 8;
peakwet = 12;

averagewet = mean(averageMth(wetseason, :));
averagedry = mean(averageMth(dryseason, :));
peakwetseason = averageMth(peakwet, :);
peakdryseason = averageMth(peakdry, :);
x = [averagewet; averagedry; peakwetseason;peakdryseason];
label = {'wet season';'dry season';'peak wet';'peak dry'};
T = table(label, x(:,1), x(:,2), x(:,3), x(:,4));
T.Properties.VariableNames = {'Month','SurfaceWater - Ice ','GW - Ice','SurfaceWater - No Ice ', 'GW - No Ice'}
writetable(T, strcat(figdir, 'WetDryFluxes_IceNoIce.csv'))

% wet dry ratios
a1 = round(100 - (x(:,3)*100./(x(:,1))));
a2 = round(100 - (x(:,4)*100./(x(:,2))));
a3 = round(100 - ((x(:,4)+x(:,3))*100./((x(:,2)+x(:,1)))));

T = table(label, a1, a2, a3);
T.Properties.VariableNames = {'Month','SurfaceWater','GW','Streamflow'}
writetable(T, strcat(figdir, 'WetDryFluxes_IceNoIce_DecreaseRatios.csv'))

%% look at how much leakage is occuring
load('CRHM\output\Cuchi_20230823.mat', 'basingw','basinflow')
load('CRHM\output\Cuchi_NoLeak_20230823.mat', 'basingw_noleak', 'basinflow_noleak','time')
load('CRHM\output\Cuchi_NoIce_20230823.mat', 'basinflow_noice','basingw_noice', 'time')


streamflow_noleak = basinflow_noleak+ basingw_noleak;
streamflow_withleak = basinflow + basingw;
streamflow_noice =  basinflow_noice+ basingw_noice;
leakage = sum(streamflow_noleak)*100./sum(streamflow_withleak)-100
icecontrib = sum(streamflow_noice)*100./sum(streamflow_withleak)-100

x = [leakage, icecontrib];
T = table(x(:,1), x(:,2));
T.Properties.VariableNames = {'Leakage','Ice Contribution'}
writetable(T, strcat(figdir, 'Leakage_IceContribution.csv'))

%% Calculate percentages

% table of all the montlhy values: surface water, gw, icemelt, snowmelt,
% rainfall runoff, et
T = timetable(time, basinflow/3600, basingw/3600,  rr_basin, swemelt_basin, icemelt_basin, et_basin); % streamflow, m3/s
TT = retime(T, 'monthly', 'sum');

% Extract the month from the time variable
months = month(TT.time);
monthNames = datestr(TT.time, 'mmm');
% Find the indices of January months
clear averageMth
for i = 1:12
mthIndices = find(months == i);
% Extract the monthly data
Data = TT(mthIndices, :);
% Calculate the average of data
averageMth (i,:) = varfun(@mean, Data(:, 1:6));
end 

% montlhy values
averageMth = table2array(averageMth);
averageMth = [averageMth; sum(averageMth)]
mth = {'January'; 'February'; 'March'; 'April'; 'May'; 'June'; 'July'; 'August'; 'September'; 'October'; 'November'; 'December';' Sum'};
T = table(mth, averageMth(:,1),averageMth(:,2), averageMth(:,3), averageMth(:,4),averageMth(:,5), averageMth(:,6));
T.Properties.VariableNames = {'Month','SurfaceWater','GW','RainfallRunoff', 'Snowmelt','Icemelt','ET'}
writetable(T, strcat(figdir, 'MonthlyFluxes.csv'))

% for wet and dry season insetad
wetseason = [10,11,12,1,2,3,4];
dryseason= [5:9];
peakdry = 8;
peakwet = 12;

averagewet = mean(averageMth(wetseason, :));
averagedry = mean(averageMth(dryseason, :));
peakwetseason = averageMth(peakwet, :);
peakdryseason = averageMth(peakdry, :);
x = round([averagewet; averagedry; peakwetseason;peakdryseason]);
label = {'wet season';'dry season';'peak wet';' peak dry'};
T = table(label, x(:,1), x(:,2), x(:,3), x(:,4), x(:,5), x(:,6));
T.Properties.VariableNames = {'Month','SurfaceWater','GW','RainfallRunoff', 'Snowmelt','Icemelt','ET'}
writetable(T, strcat(figdir, 'WetDryFluxes.csv'))

% wet dry ratios
a = sum(x(:,1:2),2);
b1 = round(x(:,1)./a *100);
b2 = round(x(:,2)./a *100);
a =sum(x(:,3:6),2);
b3 = round(x(:,3)./a*100);
b4 = round(x(:,4)./a*100);
b5= round(x(:,5)./a*100);
b6= round(x(:,6)./a*100);

T = table(label, b1, b2, b3, b4, b5, b6);
T.Properties.VariableNames = {'Month','SurfaceWater','GW','RainfallRunoff', 'Snowmelt','Icemelt','ET'}
writetable(T, strcat(figdir, 'WetDryFluxes_Ratios.csv'))

% Annual
T = timetable(time, basinflow/3600, basingw/3600,  rr_basin, swemelt_basin, icemelt_basin, et_basin); % streamflow, m3/s
TT = retime(T, 'yearly', 'sum');
averageYr = mean(table2array(TT));
x = averageYr
a = sum(x(:,1:2),2);
b1 = round(x(:,1)./a *100);
b2 = round(x(:,2)./a *100);
a =sum(x(:,3:6),2);
b3 = round(x(:,3)./a*100);
b4 = round(x(:,4)./a*100);
b5= round(x(:,5)./a*100);
b6= round(x(:,6)./a*100);
x = round(x/12)
label = {'Ratio';'Annual (mm/month)'};

T = table(label, [b1; x(1)], [b2; x(2)], [b3;x(3)], [b4;x(5)], [b5;x(5)], [b6; x(6)]);
T.Properties.VariableNames = {'Month','SurfaceWater','GW','RainfallRunoff', 'Snowmelt','Icemelt','ET'}
writetable(T, strcat(figdir, 'WetDryFluxes_Ratios_Yr.csv'))

%% Evapotraspiration values in mm/day
% mm per day for hrub 78 and 
T = timetable(time, hru_actet(:, [5,6,8,14]));
TT = retime(T, 'daily', 'sum');
averageET = mean(table2array(TT));


ET= table2array(TT);
timed = TT.time;
figure
plot(timed, smooth(ET(:,1), 2));hold on
plot(timed, smooth(ET(:,2), 2));
plot(timed, smooth(ET(:,3), 2));
plot(timed, smooth(ET(:,4), 2));

legend ('Glacier Acc. ','Glac. Abl.','Valley','Hillside')
xlim ([datetime('01-Nov-2017') datetime('01-Nov-2019')])
ylabel('Evapotranspiration (mm day^{-1})')
figname ='DailyET';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))



%% et values montlhy
% table of all the montlhy values: surface water, gw, icemelt, snowmelt,
% rainfall runoff, et
T = timetable(time, hru_actet(:, 8), hru_actet(:, 12)); % streamflow, m3/s
TT = retime(T, 'daily', 'sum');
TTT = retime(TT, 'monthly', 'mean');
TTT = TTT(7:end, :);

plot(TTT.time, TTT.Var1); hold on
plot(TT.time, TT.Var1)
% Extract the month from the time variable
months = month(TTT.time);
monthNames = datestr(TTT.time, 'mmm');
% Find the indices of January months
clear averageMth
for i = 1:12
mthIndices = find(months == i);
% Extract the monthly data
Data = TTT(mthIndices, :);
% Calculate the average of data
averageMth (i,:) = varfun(@mean, Data(:, 1:2));
end 

% montlhy values
averageMth = table2array(averageMth);
averageMth = [averageMth; sum(averageMth)]
mth = {'January'; 'February'; 'March'; 'April'; 'May'; 'June'; 'July'; 'August'; 'September'; 'October'; 'November'; 'December';' Sum'};
T = table(mth, averageMth(:,1), averageMth(:,2));
T.Properties.VariableNames = {'Month','ET8', 'ET14'}
writetable(T, strcat(figdir, 'ET_MonthlyFluxes_HRU8_14_mmday.csv'))

% for wet and dry season insetad
wetseason = [10,11,12,1,2,3,4];
dryseason= [5:9];
peakdry = 4;
peakwet = 1;

averagewet = mean(averageMth(wetseason, :));
averagedry = mean(averageMth(dryseason, :));
peakwetseason = averageMth(peakwet, :);
peakdryseason = averageMth(peakdry, :);
x = round([averagewet; averagedry; peakwetseason;peakdryseason], 3);
label = {'wet season';'dry season';'peak wet';' peak dry'};
T = table(label, x(:,1), x(:,2));
T.Properties.VariableNames = {'Month','ET8', 'ET14'}
writetable(T, strcat(figdir, 'ET_WetDryFluxes_HRU8_14.csv'))


%% a supplementray figure of the snowpack
hru_elev = [5276 5495 5020 4793 5251 4893 4290 4045 5152 5138 4762 4358 4778 4450 4976 4666 4299 4780 4298 ];
hru_elev([2,3,4,10,16])

load('CRHM\output\Cuchi_20230823.mat', 'SWE' ,'time')
figure
subplot(2,1,1)
plot(time, SWE(:, [2,3,4,10,16]))
legend ('2 (5495m)','3 (5020m)','4 (4793m)','10 (5138m)','16 (4666m)', 'orientation','horizontal', 'location','northoutside')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')])
ylim([0 700])
ylabel ('SWE (mm w.e.)')
text (datetime('10-Jul-2014'), 620, '(a)')

subplot(2,1,2)
plot(time, SWE(:, [2,3,4,10,16]))
xlim ([datetime('01-Nov-2017') datetime('01-Jun-2018')])
ylim([0 60])
text (datetime('5-Nov-2017'),55, '(b)')
ylabel ('SWE (mm w.e.)')

figname ='SWEexample';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%% Basin facts
hru_elev = [5276 5495 5020 4793 5251 4893 4290 4045 5152 5138 4762 4358 4778 4450 4976 4666 4299 4780 4298 ];
ratio_glacier = sum(hru_area([2,3,5,6,9]))./sum(hru_area)

%% Glacier melt values
load('E\CRHM\output\Cuchi_20230823.mat', 'ice', 'time', 'SWEmelt', 'icemelt', 'firnmelt', 'time')
plot(time, ice(:,[2,3,5,6,9]))
legend ('2','3','5','6','9')

t1 = find(time=='01-Jul-2014')
t2 = find(time=='01-Jul-2016')
hru_elev ([6,3,9])
melt_hru6 = (ice(t2, 6) - ice(t1,6))/1000
melt_hru3 = (ice(t2, 3) - ice(t1,3))/1000
melt_hru9 = (ice(t2, 9) - ice(t1,9))/1000

% find annual melt rate
dt = datevec(time)
dt_yr = find( dt(:,2)==01 & dt(:,3) ==1 & dt(:,4) ==0)
time(dt_yr)
ice_sub = ice(dt_yr, [6,3,9,5,2])
for i = 1:length(ice_sub)-1
melt_yr(i, :) = (ice_sub(i+1, :)-ice_sub(i, :))./1000
end 
melt_yr(length(ice_sub), :) = mean(melt_yr)
yr = {'2014';'2015';'2016';'2017';'2018';'2019';'average'}
T = table(yr,melt_yr(:,1),melt_yr(:,2),melt_yr(:,3),melt_yr(:,4),melt_yr(:,5))
T.Properties.VariableNames = {'Year','HRU6, 4893m','HRU3, 5020', 'HRU9, 5152m','HRU5, 5251m','HRU2, 5495m'}
writetable(T, strcat(figdir, 'AnnualIceMelt_perHRU.csv'))

% Find melt between 23 June- 10 July 2014
swemelt = SWEmelt./24;
icemelt = icemelt./24 ;
firnmelt = firnmelt./24;

t1 = find(time=='23-Jun-2014')
t2 = find(time=='10-Jul-2014')
melt_hru3 = (sum(icemelt(t1:t2, 3)) + sum(swemelt(t1:t2,3))+  sum(firnmelt(t1:t2,3))) /0.7/10
melt_hru6 = (sum(icemelt(t1:t2, 6)) + sum(swemelt(t1:t2,6))+  sum(firnmelt(t1:t2,6)))/0.7/10
melt_hru9 = (sum(icemelt(t1:t2, 9)) + sum(swemelt(t1:t2,9))+  sum(firnmelt(t1:t2,9)))/0.7/10
T = table(melt_hru3, melt_hru6, melt_hru9)
T.Properties.VariableNames = {'HRU3, 5020','HRU6, 4893m', 'HRU9, 5152m'}
writetable(T, strcat(figdir, 'IceMelt_perHRU_23Jun2014_10Jul2014.csv'))

