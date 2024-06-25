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

cd 'G:\11_CRHM_cuchi\'
figdir = 'G:\11_CRHM_cuchi\fig\modeleval\'; % where figures are saved
addpath 'G:\11_CRHM_cuchi\functions\' % local directory with CRHM model results

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
b = find(Qt =='01-Oct-2016'); % '07-Dec-2016');
Q(a:b)=nan;

%% Load modelled streamflow
load('CRHM\output\v9\Cuchi_20230823.mat', 'basinflow','basingw', 'time')
a = find(time =='26-June-2014');

% cut the first year (spin up)
basingw(1:a)=nan;
basinflow(1:a)=nan;

% Daily
T = timetable(time-days(0), basinflow/3600, basingw/3600);
TT = retime(T,'daily', 'mean');
bf_dd = TT.Var1;
gw_dd = TT.Var2;
t_dd = TT.Time;

% weekly
T = timetable(time-days(0), basinflow/3600, basingw/3600);
TT = retime(T,'weekly', 'mean');
bf_ww = TT.Var1;
gw_ww = TT.Var2;
t_ww = TT.Time;
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
plot(Qt, bf_d+gw_d); hold on
plot(Qt, Q)
plot(Qt, gw_d)
xlabel ('Time')
ylabel ('Streamflow (m3/s)')
legend ('streamflow - modelled','streamflow - measured', 'gw - modelled')
ylim([0 8])
% Calculate Performance Metrics
bf_mod  = bf_d+gw_d; % groundwater and surface water
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
%% 
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
load('G:\11_CRHM_cuchi\CRHM\output\v9\Cuchi_20230823.mat', 'SWEmelt', 'basinflow', 'basingw', 'firnmelt', 'icemelt', 'net_rain_org','hru_actet', 'time', 'hru_snow','hru_rain')
rainfallrunoff = net_rain_org;
swemelt = SWEmelt/24;% hourloy values
icemelt = (icemelt+firnmelt)/24; % hourloy values
hru_area = [1.528 1.322 1.237 0.839 7.366 4.056 3.492 3.663 3.201 4.633 9.491 2.721 2.418 2.008 6.249 10.62 1.903 0.14 0.45];
ratio_area = hru_area./sum(hru_area);
% flow in m3s


% Calculate for swemelt
[swemelt_m3s, swemelt_basinavg_mmwe_day] = calculate_melt(swemelt, hru_area);

% Calculate for icemelt
[icemelt_m3s, icemelt_basinavg_mmwe_day] = calculate_melt(icemelt, hru_area);

% Calculate for rainfall runoff
[rr_m3s, rr_basinavg_mmwe_day] = calculate_melt(rainfallrunoff, hru_area);

% Calculate for ET
[et_m3s, et_basinavg_mmwe_day] = calculate_melt(hru_actet, hru_area);

% for mm w.e. in basin
T =  timetable (time, swemelt_basinavg_mmwe_day, icemelt_basinavg_mmwe_day, rr_basinavg_mmwe_day, et_basinavg_mmwe_day);
TT = retime(T, 'daily','mean');
rrbasin_dd = TT. rr_basinavg_mmwe_day;
swebasin_dd = TT.swemelt_basinavg_mmwe_day;
icebasin_dd = TT.icemelt_basinavg_mmwe_day;
etbasin_dd = TT.et_basinavg_mmwe_day;
t_dd = TT.time;

T =  timetable (time, swemelt_basinavg_mmwe_day, icemelt_basinavg_mmwe_day, rr_basinavg_mmwe_day, et_basinavg_mmwe_day);
TT = retime(T, 'weekly','mean');
rrbasin_ww = TT. rr_basinavg_mmwe_day;
swebasin_ww = TT.swemelt_basinavg_mmwe_day;
icebasin_ww = TT.icemelt_basinavg_mmwe_day;
etbasin_ww = TT.et_basinavg_mmwe_day;
t_ww = TT.time;

% for m3/s
T =  timetable (time, swemelt_m3s, icemelt_m3s, rr_m3s, et_m3s, basinflow, basingw);
TT = retime(T, 'daily','mean');
rr_m3s_dd = TT. rr_m3s;
swe_m3s_dd = TT.swemelt_m3s;
ice_m3s_dd = TT.icemelt_m3s;
et_m3s_dd = TT.et_m3s;
flow_m3s_dd = TT.basinflow/3600 + TT.basingw/3600;
t_dd = TT.time;

T =  timetable (time, swemelt_m3s, icemelt_m3s, rr_m3s, et_m3s);
TT = retime(T, 'weekly','mean');
rr_m3s_ww = TT. rr_m3s;
swemelt_m3s_ww = TT.swemelt_m3s;
icemelt_m3s_ww = TT.icemelt_m3s;
et_m3s_ww = TT.et_m3s;
t_ww = TT.time;


figure
area(t_dd, [ smooth(ice_m3s_dd, 5),  smooth(swe_m3s_dd,5), smooth(rr_m3s_dd,5)], 'Linestyle', 'none'); hold on
plot(Qt, Q, 'k', 'linewidth', 0.8)
plot(t_dd,flow_m3s_dd, ':k',  'linewidth', 0.8)
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
load('G:\11_CRHM_cuchi\CRHM\output\v9\Cuchi_NoIce_20230823.mat', 'basinflow_noice','basingw_noice', 'time')
a = find(time =='26-June-2014');
% cut the first year 9spin upP
basingw_noice(1:a)=nan;
basinflow_noice(1:a)=nan;

t = timetable(time-days(2), basinflow_noice/3600, basingw_noice/3600);
tt = retime(t, 'weekly', 'mean');
bf_ww_noice = tt.Var1;
gw_ww_noice = tt.Var2;
t_ww = tt.Time;

t = timetable(time-days(2), basinflow_noice/3600, basingw_noice/3600);
tt = retime(t, 'daily', 'mean');
bf_dd_noice = tt.Var1;
gw_dd_noice = tt.Var2;
t_dd = tt.Time;

figure
plot(t_dd, bf_dd_noice); hold on
plot(t_dd, bf_dd)

a= gw_ww_noice./(bf_ww_noice+ gw_ww_noice);
b = gw_ww ./ (bf_ww+gw_ww);
figure
plot(a); hold on; plot (b)
%% Import the no groundwater simulation
load('CRHM\output\v9\Cuchi_NoGW_20230823.mat', 'basinflow_nogw', 'basingw_nogw', 'time')
a = find(time =='26-June-2014');
% cut the first year 9spin upP
basinflow_nogw(1:a)=nan;
basingw_nogw(1:a)=nan;

T = timetable(time-days(2), basinflow_nogw/3600, basingw_nogw/3600);
TT= retime(T, 'daily', 'mean');
bf_dd_nogw = TT.Var1;
gw_dd_nogw = TT.Var2;
t_dd = TT.Time;

T = timetable(time-days(2), basinflow_nogw/3600, basingw_nogw/3600);
TT= retime(T, 'weekly', 'mean');
bf_ww_nogw = TT.Var1;
gw_ww_nogw = TT.Var2;
t_ww = TT.Time;

figure
plot(t_ww, bf_ww_nogw); hold on
plot(t_ww, gw_ww_nogw)


%% lot streamflow components
lw =1.1;
fig = figure('units','inches','outerposition',[0 0 8 7]);
subplot(3,1,1)
a = area(t_dd, [gw_dd,bf_dd]); hold on
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
a = area(t_ww, [gw_ww_nogw,bf_ww_nogw]); hold on
a(1).FaceColor = [0.45 0.45 0.45];
a(2).FaceColor = [0 76 153]./255;
a(1).EdgeColor = [0.45 0.45 0.45];
a(2).EdgeColor = [0 76 153]./255;
plot(t_ww, gw_ww(1:end-1),'-k', 'Linewidth', lw,'color', [0.7 0.7 0.7]);
plot(t_ww, bf_ww(1:end-1)+gw_ww(1:end-1),'-k', 'Linewidth', lw, 'color', [102 178 255]./255);
legend('GW - low GW', 'Ov. + Vad. - low GW','GW (baseline)','Ov. + Vad. (baseline)', 'orientation','horizontal');
ylabel ({'Weekly Mean' ;'Streamflow, m^3/s'});
ylim([0 6]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
text (datetime('10-Jul-2014'), 5.6, '(b)');


subplot(3,1,3)
a = area(t_ww, [icebasin_ww(1:end-1), swebasin_ww(1:end-1), rrbasin_ww(1:end-1)]*10, 'Linestyle', 'none'); hold on
b = area(t_ww, [-etbasin_ww(1:end-1)]*10,  'Linestyle', 'none'); 
legend ('icemelt runoff','snowmelt runoff','rainfall runoff',  'ET' ,'orientation','horizontal');
ylabel ({'Weekly Basin-Averaged';'Hydrological Components (mm)'});
a(3).FaceColor = [0 76 153]./255;
a(2).FaceColor = [255 128 0]./255;
a(1).FaceColor = [153 153 255]./255;
b(1).FaceColor = [41 135 0]./255;
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim([-30 110]);
text (datetime('10-Jul-2014'), 99, '(c)');

figname ='F3_Basinflow_SurfaceSubsurface_CDA_withflowcomponent';
saveas (gcf, strcat( figdir, figname, '.pdf'));
saveas (gcf, strcat(figdir, figname, '.png'));
savefig(gcf, strcat(figdir, figname));

%% Just no ice
lw =1.1;
fig = figure('units','inches','outerposition',[0 0 8 4]);

a = area(t_ww, [gw_ww_noice,bf_ww_noice]); hold on
a(1).FaceColor = [0.45 0.45 0.45];
a(2).FaceColor = [0 76 153]./255;
a(1).EdgeColor = [0.45 0.45 0.45];
a(2).EdgeColor = [0 76 153]./255;
plot(t_ww, gw_ww(1:end-1),'-k', 'Linewidth', lw+0.3,'color', [0.7 0.7 0.7]);
plot(t_ww, bf_ww(1:end-1)+gw_ww(1:end-1),'-k', 'Linewidth', lw+0.3, 'color', [102 178 255]./255);
legend('GW - No Ice', 'Ov. + Vad. - No Ice', 'GW (baseline)','Ov. + Vad. (baseline)', 'orientation','horizontal');
ylabel ({'Weekly Mean' ;'Streamflow, m^3/s'});
ylim([0 6]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);

figname ='SW_GWsim_NoIce';
grid on
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))


%% lot streamflow components
lw =1.1;
fig = figure('units','inches','outerposition',[0 0 8 9]);
subplot(4,1,1)
a = area(t_dd, [gw_dd,bf_dd]); hold on
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
subplot(4,1,2)
a = area(t_ww, [gw_ww_noice,bf_ww_noice]); hold on
a(1).FaceColor = [0.45 0.45 0.45];
a(2).FaceColor = [0 76 153]./255;
a(1).EdgeColor = [0.45 0.45 0.45];
a(2).EdgeColor = [0 76 153]./255;
plot(t_ww, gw_ww(1:end-1),':k', 'Linewidth', lw+0.3,'color', [0.7 0.7 0.7]);
plot(t_ww, bf_ww(1:end-1)+gw_ww(1:end-1),':k', 'Linewidth', lw+0.3, 'color', [102 178 255]./255);
legend('GW - No Ice', 'Ov. + Vad. - No Ice', 'GW','Ov. + Vad.', 'orientation','horizontal');
ylabel ({'Weekly Mean' ;'Streamflow, m^3/s'});
ylim([0 6]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
text (datetime('10-Jul-2014'), 5.6, '(b)');

subplot(4,1,3)
a = area(t_ww, [gw_ww_nogw,bf_ww_nogw]); hold on
a(1).FaceColor = [0.45 0.45 0.45];
a(2).FaceColor = [0 76 153]./255;
a(1).EdgeColor = [0.45 0.45 0.45];
a(2).EdgeColor = [0 76 153]./255;
plot(t_ww, gw_ww,':k', 'Linewidth', lw+0.3,'color', [0.7 0.7 0.7]);
plot(t_ww, bf_ww+gw_ww,':k', 'Linewidth', lw+0.3, 'color', [102 178 255]./255);
legend('GW - low GW', 'Ov. + Vad. - low GW','GW','Ov. + Vad.', 'orientation','horizontal');
ylabel ({'Weekly Mean' ;'Streamflow, m^3/s'});
ylim([0 8]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
text (datetime('10-Jul-2014'), 7.3, '(c)');

subplot(4,1,4)
a = area(t_ww, [icebasin_ww(1:end-1), swebasin_ww(1:end-1), rrbasin_ww(1:end-1)]*10, 'Linestyle', 'none'); hold on
b = area(t_ww, [-etbasin_ww(1:end-1)]*10,  'Linestyle', 'none'); 
legend ('icemelt runoff','snowmelt runoff','rainfall runoff',  'ET' ,'orientation','horizontal');
ylabel ({'Weekly Basin-Averaged';'Hydrological Components (mm)'});
a(3).FaceColor = [0 76 153]./255;
a(2).FaceColor = [255 128 0]./255;
a(1).FaceColor = [153 153 255]./255;
b(1).FaceColor = [41 135 0]./255;
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
ylim([-30 110]);
text (datetime('10-Jul-2014'), 97, '(d)');

figname ='F3_Basinflow_SurfaceSubsurface_CDA_withflowcomponent_v2';
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
annual = [1:12];

averagewet = mean(averageMth(wetseason, :));
averagedry = mean(averageMth(dryseason, :));
peakwetseason = averageMth(peakwet, :);
peakdryseason = averageMth(peakdry, :);
avgannual=mean(averageMth(annual, :));

x = [averagewet; averagedry; peakwetseason;peakdryseason;avgannual];
label = {'wet season';'dry season';'peak wet';'peak dry'; 'annual'};
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

% can I make a figute of SW-GW sepearte as percent for with ice and no ice?

frac_sw = basinflow*100./(basinflow+basingw);
frac_gw =basingw*100./(basinflow+basingw);
frac_gw_noice = basingw_noice*100./(basinflow_noice+basingw_noice);
T = timetable(time, frac_sw, frac_gw, frac_gw_noice);
TT = retime(T,'weekly','mean');
tw = TT.time;
frac_sw=TT.frac_sw;
frac_gw=TT.frac_gw;
frac_gw_noice=TT.frac_gw_noice;

lw =1.1;
fig = figure('units','inches','outerposition',[0 0 8 4]);
a = area(tw, [frac_gw, frac_sw], 'edgecolor', 'none'); hold on
a(1).FaceColor = [0.7 0.7 0.7];
a(2).FaceColor = [102 178 255]./255;
a(1).EdgeColor = [0.7 0.7 0.7];
a(2).EdgeColor = [102 178 255]./255;
plot(tw,frac_gw_noice , '-k', 'Linewidth', lw); hold on
legend ('GW', 'Ov + Vad', 'No Ice', 'orientation','horizontal', 'location','best');
ylabel ({'Weekly Mean' ;'Groundwater Fraction'});
ylim([0 100])
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')]);
text (datetime('10-Jul-2014'), 95, '(a)');
figname ='FractionSWGW';
grid on
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

T = timetable(tw,frac_gw, frac_gw_noice); % streamflow, m3/s
TT = retime(T, 'monthly', 'mean');
months = month(TT.tw);
monthNames = datestr(TT.tw, 'mmm');
% Find the indices of January months
clear averageMth
for i = 1:12
mthIndices = find(months == i);
% Extract the monthly data
Data = TT(mthIndices, :);
% Calculate the average of data
averageMth (i,:) = varfun(@nanmean, Data(:, 1:2));
end 

% monthly values
averageMth = table2array(averageMth);
averageMth = [averageMth; mean(averageMth)]
mth = {'January'; 'February'; 'March'; 'April'; 'May'; 'June'; 'July'; 'August'; 'September'; 'October'; 'November'; 'December';' Mean'};
T = table(mth, averageMth(:,1),averageMth(:,2));
T.Properties.VariableNames = {'Month','Frac SW ','Frac SW no ice'}
writetable(T, strcat(figdir, 'MonthlyFrac_IceNoIce.csv'))

% for wet and dry season insetad
wetseason = [10,11,12,1,2,3,4];
dryseason= [5:9];
peakdry = [9:12];
peakwet = 12;
annual = [1:12];

averagewet = mean(averageMth(wetseason, :));
averagedry = mean(averageMth(dryseason, :));
peakwetseason = averageMth(peakwet, :);
peakdryseason =mean( averageMth(peakdry, :));
avgannual=mean(averageMth(annual, :));

x = [averagewet; averagedry; peakwetseason;peakdryseason;avgannual];
label = {'wet season';'dry season';'peak wet';'peak dry'; 'annual'};
T = table(label, x(:,1), x(:,2));
T.Properties.VariableNames = {'Month','SurfaceWater - Ice ','SW - No Ice'}
writetable(T, strcat(figdir, 'WetDryRatioSW_IceNoIce.csv'))



%% look at how much leakage is occuring
load('CRHM\output\v9\Cuchi_20230823.mat', 'basingw','basinflow')
load('CRHM\output\v9\Cuchi_NoLeak_20230823.mat', 'basingw_noleak', 'basinflow_noleak','time')
load('CRHM\output\v9\Cuchi_NoIce_20230823.mat', 'basinflow_noice','basingw_noice', 'time')

streamflow_noleak = basinflow_noleak+ basingw_noleak;
streamflow_withleak = basinflow + basingw;

figure;
plot(time, streamflow_noleak/3600); hold on
plot(time, streamflow_withleak/3600)
plot(Qt, Q)
streamflow_noice =  basinflow_noice + basingw_noice;
leakage = sum(streamflow_noleak)*100./sum(streamflow_withleak)-100
icecontrib = sum(streamflow_noice)*100./sum(streamflow_withleak)-100

x = [leakage, icecontrib];
T = table(x(:,1), x(:,2));
T.Properties.VariableNames = {'Leakage','Ice Contribution'}
writetable(T, strcat(figdir, 'Leakage_IceContribution.csv'))

streamflow_noleak = basinflow_noleak+ basingw_noleak;
streamflow_withleak = basinflow + basingw;
streamflow_noice =  basinflow_noice+ basingw_noice;

leakage = sum(streamflow_noleak)*100./sum(streamflow_withleak)-100
icecontrib = sum(streamflow_noice)*100./sum(streamflow_withleak)-100

plot(time, basingw); hold on
plot(time, basingw_noleak)
leakage_gw = sum(basingw_noleak)*100./sum(basingw)-100
icecontrib_gw = sum(basingw_noice)*100./sum(basingw)-100

x = [leakage_gw, icecontrib_gw];

T = table(x(:,1), x(:,2));
T.Properties.VariableNames = {'Leakage','Ice Contribution'}
writetable(T, strcat(figdir, 'LeakageGW_IceContribution.csv'))

%% Calculate percentages

% table of all the montlhy values: surface water, gw, icemelt, snowmelt,
% rainfall runoff, et
T = timetable(t_dd, bf_dd/3600, gw_dd/3600,  rrbasin_dd, swebasin_dd, icebasin_dd, etbasin_dd); % streamflow, m3/sTT = retime(T, 'monthly', 'sum');
TT = retime(T, 'monthly', 'sum');
% Extract the month from the time variable
months = month(TT.t_dd);
monthNames = datestr(TT.t_dd, 'mmm');
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
T = timetable(t_dd, bf_dd/3600, gw_dd/3600,  rrbasin_dd, swebasin_dd, icebasin_dd, etbasin_dd); % streamflow, m3/sTT = retime(T, 'monthly', 'sum');
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
T = timetable(time, hru_actet(:, [5,6,8,16]));
TT = retime(T, 'daily', 'sum');
averageET = mean(table2array(TT));
lw = 0.8
ET= table2array(TT);
timed = TT.time;
figure
plot(timed, smooth(ET(:,1), 2), 'linewidth', lw);hold on
plot(timed, smooth(ET(:,2), 2), 'linewidth', lw);
plot(timed, smooth(ET(:,3), 2), 'linewidth', lw);
plot(timed, smooth(ET(:,4), 2), 'linewidth', lw);
legend ('Glacier Acc. - 5 ','Glac. Abl. -6','Valley -8','Hillside -14')
xlim ([datetime('01-Nov-2017') datetime('01-Nov-2019')])
ylabel('Evapotranspiration (mm day^{-1})')
figname ='DailyET';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))


T = timetable(time, hru_actet(:, [5,6,8,17]));
TT = retime(T, 'daily', 'sum');
TT = retime(TT, 'weekly', 'mean');

averageET = mean(table2array(TT));


ET= table2array(TT);
timed = TT.time;
figure
plot(timed, ET(:,1), 'linewidth', lw);hold on
plot(timed, ET(:,2), 'linewidth', lw);
plot(timed, ET(:,3), 'linewidth', lw);
plot(timed, ET(:,4), 'linewidth', lw);

legend ('Glacier Acc. - 5 ','Glac. Abl. -6','Valley -8','Hillside -14')
xlim ([datetime('01-Nov-2017') datetime('01-Nov-2019')])
ylabel('Evapotranspiration (mm day^{-1})')
figname ='WeeklyET';
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

% monthly values
averageMth = table2array(averageMth);
averageMth = [averageMth; mean(averageMth);sum(averageMth)]
mth = {'January'; 'February'; 'March'; 'April'; 'May'; 'June'; 'July'; 'August'; 'September'; 'October'; 'November'; 'December';'Average';'Sum'};
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
hru_elev([1,4,10,15])

load('CRHM\output\v9\Cuchi_20230823.mat', 'SWE' ,'time')
figure
subplot(2,1,1)
plot(time, SWE(:, [1,4,10,15]), 'linewidth', 0.8)
legend ('1 (5276m)','4 (4793m)','10 (5138m)','15 (4976m)', 'orientation','horizontal', 'location','northoutside')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')])
ylim([0 700])
ylabel ('Off-Glacier SWE (mm w.e.)')
text (datetime('10-Jul-2014'), 620, '(a)')

subplot(2,1,2)
plot(time, SWE(:, [2,3,5,6, 9]), 'linewidth', 0.8)
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')])
ylim([0 900])
legend ('2 (5495m)','3 (5020m)','5 (5251m)','6 (4893m)','9 (5152m) ', 'orientation','horizontal', 'location','northoutside')
text (datetime('5-Nov-2017'),55, '(b)')
ylabel ('On-Glacier SWE (mm w.e.)')

figname ='SWEexample';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%% a supplementray figure of the snowpack on glacier
hru_elev = [5276 5495 5020 4793 5251 4893 4290 4045 5152 5138 4762 4358 4778 4450 4976 4666 4299 4780 4298 ];
hru_elev([2,3,5,6,9])

load('CRHM\output\v9\Cuchi_20230823.mat', 'SWE' ,'time')
figure
subplot(2,1,1)
plot(time, SWE(:, [2,3,5,6,9]))
legend ('2 (5495m)','3 (5020m)','5 (5251m)','6 (4893m)','9 (5152m)', 'orientation','horizontal', 'location','northoutside')
xlim ([datetime('26-Jun-2014') datetime('01-Apr-2020')])
ylim([0 700])
ylabel ('SWE (mm w.e.)')
text (datetime('10-Jul-2014'), 620, '(a)')

subplot(2,1,2)
plot(time, SWE(:, [2,3,4, 5,6,9]))
xlim ([datetime('01-Nov-2017') datetime('01-Jun-2018')])
ylim([0 60])
text (datetime('5-Nov-2017'),55, '(b)')
ylabel ('SWE (mm w.e.)')

figname ='SWEexample_onglacier';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
%% Basin facts
hru_elev = [5276 5495 5020 4793 5251 4893 4290 4045 5152 5138 4762 4358 4778 4450 4976 4666 4299 4780 4298 ];
ratio_glacier = sum(hru_area([2,3,5,6,9]))./sum(hru_area)

%% Glacier melt values
load('CRHM\output\v9\Cuchi_20230823.mat', 'ice', 'time', 'SWEmelt', 'icemelt', 'firnmelt', 'time')
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
melt_hru3 = (sum(icemelt(t1:t2, 3)) + sum(swemelt(t1:t2,3))+  sum(firnmelt(t1:t2,3))) /0.7/10;
melt_hru6 = (sum(icemelt(t1:t2, 6)) + sum(swemelt(t1:t2,6))+  sum(firnmelt(t1:t2,6)))/0.7/10;
melt_hru9 = (sum(icemelt(t1:t2, 9)) + sum(swemelt(t1:t2,9))+  sum(firnmelt(t1:t2,9)))/0.7/10;
T = table(melt_hru3, melt_hru6, melt_hru9)
T.Properties.VariableNames = {'HRU3, 5020','HRU6, 4893m', 'HRU9, 5152m'}
writetable(T, strcat(figdir, 'IceMelt_perHRU_23Jun2014_10Jul2014.csv'))




