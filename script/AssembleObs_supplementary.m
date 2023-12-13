% Pre-processing meteorological data to create observation file for
% Quilcayhuanca valley CRHM model
% By Caroline Aubry-Wake
% May 2023

% this script create a complete observation file form available data in the
% Cordillera Blanca and exports it as a CRHM observation file. It also
% proprocess the Casa de agua streamflow to make it as an obervation file
% for CRHM. The processing of the MET data is shown in the supplementary
% file 1.

% the data analysed here is downloaded from 
%Mateo, E. I., B. G. Mark, R. Å. Hellström, M. Baraer, J. M. McKenzie, T. Condom, A. C. Rapre, G. Gonzales, J. Q. Gómez, R. C. C. Encarnación (2022). High temporal resolution hydrometeorological data collected in the tropical Cordillera Blanca, Peru (2004-2020), HydroShare, https://doi.org/10.4211/hs.35a670e6c5824ff89b3b74fe45ca90e0

%% set up
figdir = 'F:\11_CRHM_cuchi\fig\obs\'
folderPath = 'F:\11_CRHM_cuchi\data\Wx\';  % find the met files
files = dir([folderPath '*.csv']);  % Get a list of all CSV files in the folder

numFiles = numel(files); % Get the number of CSV files
clear data
% 4 filea: Cuchillacocha, Caasa de Agua, Llaganuco Wx and LLanganuco Port

% Import each file and apply the precipitation undercatch correction. 
for i = 1:numFiles
 filepath = [folderPath files(i).name]; % Get the full file path
 d = readtable(filepath); % Read the table from the file

% Identify columns that are numeric, datetime, or duration types
numericCols = varfun(@isnumeric, d, 'OutputFormat', 'uniform');
datetimeCols = varfun(@isdatetime, d, 'OutputFormat', 'uniform');
durationCols = varfun(@isduration, d, 'OutputFormat', 'uniform');

% Find column indices to keep (remove clumn with other data type
colsToKeep = numericCols | datetimeCols | durationCols;
% Keep only desired columns
d = d(:, colsToKeep);

% Set to hourly values (mean for all, sum for precip)
dt = table2timetable(d);
dd_m = retime(dt, 'hourly', 'mean') ;
dd_s = retime(dt, 'hourly', 'sum') ;

% Get the variable names from the table
variableNames = dd_s.Properties.VariableNames;

% Find the header that contains "precp" (case-insensitive)
matchingColumnIndex = -1;
for ii = 1:numel(variableNames)
    if contains(variableNames{ii}, 'Precip', 'IgnoreCase', true)
        matchingColumnIndex = ii;
        break;
    end
end

% Retrieve the column with the header 'Precip'
dataToKeep = dd_s(:, matchingColumnIndex);
dd_combined = dd_m;
dd_combined(:, matchingColumnIndex) = dataToKeep ;

  % applying undercatch correction to the precipitation data using the winf
  % speed
   a = 0.728;
   b = 0.230;
   c = 0.336;
% find the wind variable
UmatchingColumnIndex = -1;
    for ii = 1:numel(variableNames)
        if contains(variableNames{ii}, 'WindS', 'IgnoreCase', true)
            UmatchingColumnIndex = ii;
            break;
        end
    end

   U = table2array(dd_combined(:, UmatchingColumnIndex));
   CE = (a) * exp(-b*U) + c;
   Pcorr = table2array(dd_combined(:, matchingColumnIndex)).*(1./CE);
   dd_combined = addvars(dd_combined, Pcorr, 'After', size(dd_combined, 2));

   data{i} = dd_combined; % Store the modified table in the cell array
end

% Divide the data into the four segments
CDA = data{1};
Cuchi  = data{2};
LlanWx = data{3};
LlanPort = data{4};

% Plot the correction precip for each dataset
figure
subplot(2,1,1)
plot(CDA.Datetime, CDA.Pcorr);
ylabel ({'Hourly Precipitation,';' Corrected (mm)'})
hold on

subplot(2,1,2)
plot(CDA.Datetime, cumsum(CDA.Precipitation_mm_)); hold on
plot(CDA.Datetime, cumsum(CDA.Pcorr));
ylabel ('Cumul. Precipitation (mm)')
legend ('Originial', 'Corrected', 'location', 'best')

figname ='CDA_precipitation_hourly_cumul';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
% 



%% Infill with cuchi
T = retime(CDA, Cuchi.Datetime, 'fillwithmissing')
CDA_values = T(:, 9).Variables;
Cuchi_values = Cuchi(:, 8).Variables;
CDA_time = T.Datetime;
Cuchi_time = Cuchi.Datetime;

% Find indices where both columns have a value
t = CDA_time
d  =[CDA_values, Cuchi_values];

%% Do a regresison correction;
valid_indices = find(~isnan(d(:,2))  & ~isnan(d(:,1)));
subset_d1 = d(valid_indices, 1);
subset_d2 = d(valid_indices, 2);

% Perform linear regression
coefficients = polyfit(subset_d1, subset_d2, 1);
yfit = polyval(coefficients, d(:,2))

figure
hold on
subplot(2,1,1)
plot(t, yfit); hold on; plot(t, d(:,1)); plot(t, d(:,2))
legend ('Fitted','CDA','Cuchi', 'orientation','horizontal','location','best')
ylabel({'P (mm)'});
subplot(2,1,2)
plot(t, yfit); hold on; plot(t, d(:,1)); plot(t, d(:,2))
xlim([datetime('10-Feb-2018') datetime('15-Feb-2018')])
legend ('Fitted','CDA','Cuchi', 'orientation','horizontal','location','best')
ylabel({'P (mm)'});

figname ='CDA_precip_RegressedCuchi';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all


r = corrcoef(d(:,1), yfit, 'rows', 'pairwise')
r_corr = r(2);
RMSE = rmse(d(:,1), yfit, 'omitnan')
r = corrcoef(d(:,1), d(:,2), 'rows', 'pairwise')
r_or = r(2);
RMSE_or = rmse(d(:,1), d(:,2), 'omitnan')


cuchi_filled = yfit;

figure;
scatter(d(:, 1), d(:,2), 'r.');hold on
scatter(d(:, 1), yfit, 'k.'); 
grid on
lsline
rl = refline(1,0)
rl.Color ='k'
ylabel('Cuchillacocha')
xlabel('CDA')
legend( 'Original','Corrected', 'best fit, or','best fit, corr','1:1', 'location', 'best')
text(-4, 22, strcat('y = ', num2str(round(coefficients(1), 2)), '+' , strcat(num2str(round(coefficients(2), 2)))))
text (-4, 20, strcat( 'RMSEor=', num2str(round(RMSE_or,2)), ', RMSEcorr = ', num2str(round(RMSE,2))))
text (-4, 18, strcat( 'R =', num2str(round(r_or,2)), ', Rcorr = ', num2str(round(r_corr,2))))

grid on
lsline
rl = refline(1,0)
rl.Color ='k'
figname = 'CDAvsCuchi_linearregression_p'
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

close all

% fill with missing values
CDA_filled = retime(CDA,Cuchi.Datetime, 'fillwithmissing')
CDA_filled = CDA_filled(:, 9).Variables;
CDA_filled_0 = CDA_filled;
CDA_filled(isnan(CDA_filled)) = Cuchi_values(isnan(CDA_filled))*0.85;

a = find(Cuchi_time == '26-Jun-2014')
CDAp_filled = CDA_filled(a:end);
CDA_filled_0 = CDA_filled_0(a:end);

CDA_time = Cuchi_time (a:end);


figure;
plot(CDA_time, CDAp_filled); hold on
plot(CDA_time, CDA_filled_0)
legend ('Corrected Cuchi','CDA','orientation','horizontal','location','best')
ylabel({'Precip (mm)'});
ylim([0 12])
figname ='CDA_precip_Filled';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%%
% calculate annual sum
Tt = timetable(CDA_time, CDAp_filled)
T = retime(Tt, 'yearly', 'sum')
figure
subplot(2,1,1)

X = T.CDA_time(2:end-1);
Y = T.CDAp_filled(2:end-1);
bar(X, Y)
labels = arrayfun(@(value) num2str(round(value)), Y, 'UniformOutput', false);
text(X,Y,labels,'HorizontalAlignment','center','VerticalAlignment','bottom') 
      % clears X axis data
ylabel('Annual Precipitation (mm)' )
ylim([0 1100])
subplot(2,1,2)
Tt = timetable(CDA_time, CDAp_filled)
T = retime(Tt, 'monthly', 'sum')
X = T.CDA_time;
Y = T.CDAp_filled;
bar(X, Y)

figname ='CDA_annualPrecip_bar';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

x = timetable(CDA_time,  CDAp_filled)
filename = 'F:\11_CRHM_cuchi\data\processed\InfilledP_20130705_20200406.csv'
writetimetable(x, filename);

close all
%% Air temperature
figure
hold on
subplot(2,1,1)
plot (LlanWx.Datetime, LlanWx.Temperature__C_); hold on
plot (LlanPort.Datetime, LlanPort.AirTemperature__C_)
plot (Cuchi.Datetime, Cuchi.Temperature__C_)
plot (CDA.Datetime, CDA.Temperature__C_)
ylim([-5 23])
legend ('LlanWx','LlanPort','Cuchi','CDA', 'orientation','horizontal','location','best')
subplot(2,1,2)
plot (LlanWx.Datetime, LlanWx.Temperature__C_); hold on
plot (LlanPort.Datetime, LlanPort.AirTemperature__C_)
plot (Cuchi.Datetime, Cuchi.Temperature__C_)
plot (CDA.Datetime, CDA.Temperature__C_)
ylim([-4 20])
xlim([datetime('10-Feb-2018') datetime('15-Feb-2018')])
legend ('LlanWx','LlanPort','Cuchi','CDA', 'orientation','horizontal','location','best')

figname ='CDA_temperature_allstations_hourly_2006_2020';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all


%% Regresison with Cuchillacocha
T = retime(CDA, Cuchi.Datetime, 'fillwithmissing')
CDA_values = T(:, 6).Variables;
Cuchi_values = Cuchi(:, 2).Variables;
CDA_time = T.Datetime;
Cuchi_time = Cuchi.Datetime;

% Find indices where both columns have a value
t = CDA_time
d  =[CDA_values, Cuchi_values];

%% Do a regresison correction;
valid_indices = find(~isnan(d(:,2))  & ~isnan(d(:,1)));
subset_d1 = d(valid_indices, 1);
subset_d2 = d(valid_indices, 2);

% Perform linear regression
coefficients = polyfit(subset_d2, subset_d1, 2);
yfit = polyval(coefficients, d(:,2))

figure
hold on
subplot(2,1,1)
plot(t, yfit); hold on; plot(t, d(:,1)); plot(t, d(:,2))
ylim([-5 23])
legend ('Fitted','CDA','Cuchi', 'orientation','horizontal','location','best')
ylabel({'Air Temperature ({\circ}C)'});
title ('Air Temperature, Cuchillacocha and CDA')
subplot(2,1,2)
plot(t, yfit); hold on; plot(t, d(:,1)); plot(t, d(:,2))
ylim([-4 20])
xlim([datetime('10-Feb-2018') datetime('15-Feb-2018')])
legend ('Fitted','CDA','Cuchi', 'orientation','horizontal','location','best')
ylabel({'Air Temperature ({\circ}C)'});

figname ='CDA_temperature_RegressedCuchi_2006_2020';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all


r = corrcoef(d(:,1), yfit, 'rows', 'pairwise')
r_corr = r(2);
RMSE = rmse(d(:,1), yfit, 'omitnan')
r = corrcoef(d(:,1), d(:,2), 'rows', 'pairwise')
r_or = r(2);
RMSE_or = rmse(d(:,1), d(:,2), 'omitnan')


cuchi_filled = yfit;

figure;
scatter(d(:, 1), d(:,2), 'r.');hold on
scatter(d(:, 1), yfit, 'k.'); 
grid on
lsline
rl = refline(1,0)
rl.Color ='k'
ylabel('Cuchillacocha')
xlabel('CDA')
legend( 'Original','Corrected', 'best fit, or','best fit, corr','1:1', 'location', 'best')
text(-4, 22, strcat('y = ', num2str(round(coefficients(1), 2)), '+' , strcat(num2str(round(coefficients(2), 2)))))
text (-4, 20, strcat( 'RMSEor=', num2str(round(RMSE_or,2)), ', RMSEcorr = ', num2str(round(RMSE,2))))
text (-4, 18, strcat( 'R =', num2str(round(r_or,2)), ', Rcorr = ', num2str(round(r_corr,2))))

grid on
lsline
rl = refline(1,0)
rl.Color ='k'
figname = 'CDAvsCuchi_linearregression_airtemp'
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

close all
%% Regression correctio with Llanguno
T = retime(CDA, Cuchi.Datetime, 'fillwithmissing')
CDA_values = T(:, 6).Variables;
T = retime(LlanWx, Cuchi.Datetime, 'fillwithmissing')
Llan_values = T.Temperature__C_;
CDA_time = T.Datetime;
Llan_time = Cuchi.Datetime;


% Find indices where both columns have a value
t = CDA_time
d  =[CDA_values, Llan_values];

valid_indices = ~isnan(d(:, 1)) & ~isnan(d(:, 2));

%% Infill CDA with Llan
%% Do a regresison correction;
valid_indices = find(~isnan(d(:,2))  & ~isnan(d(:,1)));
subset_d1 = d(valid_indices, 1);
subset_d2 = d(valid_indices, 2);

% Perform linear regression
coefficients = polyfit(subset_d2, subset_d1, 1);
yfit = polyval(coefficients, d(:,2))

r = corrcoef(d(:,1), yfit, 'rows', 'pairwise')
r_corr = r(2);
RMSE = rmse(d(:,1), yfit, 'omitnan')
r = corrcoef(d(:,1), d(:,2), 'rows', 'pairwise')
r_or = r(2);
RMSE_or = rmse(d(:,1), d(:,2), 'omitnan')


figure
hold on
subplot(2,1,1)
plot(t, yfit); hold on; plot(t, d(:,1)); plot(t, d(:,2))
ylim([-5 23])
legend ('Fitted','CDA','LlanWx', 'orientation','horizontal','location','best')
ylabel({'Air Temperature ({\circ}C)'});
title ('Air Temperature, Llanganuco and CDA')
subplot(2,1,2)
plot(t, yfit); hold on; plot(t, d(:,1)); plot(t, d(:,2))
ylim([-4 20])
xlim([datetime('10-Feb-2018') datetime('15-Feb-2018')])
legend ('Fitted','CDA','Llan Wx', 'orientation','horizontal','location','best')
ylabel({'Air Temperature ({\circ}C)'});

figname ='CDA_temperature_RegressedLlan_2006_2020';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

llan_filled = yfit; 

figure;
scatter(d(:, 1), d(:,2), 'r.');hold on
scatter(d(:, 1), yfit, 'k.'); 
grid on
lsline
rl = refline(1,0)
rl.Color ='k'
ylabel('Llan')
xlabel('CDA')
legend( 'Original','Regression', 'best fit, or','best fit, corr','1:1', 'location', 'best')
grid on
lsline
rl = refline(1,0)
rl.Color ='k'
text(-4, 22, strcat('y = ', num2str(round(coefficients(1), 2)), '+' , strcat(num2str(round(coefficients(2), 2)))))
text (-4, 20, strcat( 'RMSEor=', num2str(round(RMSE_or,2)), ', RMSEcorr = ', num2str(round(RMSE,2))))
text (-4, 18, strcat( 'R =', num2str(round(r_or,2)), ', Rcorr = ', num2str(round(r_corr,2))))


figname = 'CDAvsLlan_linearregression_airtemp'
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

close all
%% Infill air temp at CDA with llan thn Cuchi
CDA_filled = retime(CDA,Cuchi.Datetime, 'fillwithmissing')
CDA_values = CDA_filled(:, 6).Variables;
CDA_filled = CDA_filled(:, 6).Variables;
CDA_filled(isnan(CDA_filled)) = llan_filled(isnan(CDA_filled));
CDA_filled(isnan(CDA_filled)) = cuchi_filled(isnan(CDA_filled));

a = find(Cuchi_time == '26-Jun-2014')
CDAta_filled = CDA_filled(a:end);
CDA_time = Cuchi_time (a:end);


figure;
plot(CDA_time, CDAta_filled); hold on
legend ('Filled CDA','orientation','horizontal','location','best')
ylabel({'Air Temperature ({\circ}C)'});
title ('Air Temperature, Filled CDA')
figname ='CDA_temperature_Filled';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

% export that time series:
x = timetable(CDA_time, CDAta_filled)
filename = 'F:\11_CRHM_cuchi\data\processed\InfilledTa_20130705_20200406.csv'
writetimetable(x, filename);

t_obs = CDA_time;


%%%%%%%%% Wind speed
%% Infill u
 
figure
hold on
subplot(2,1,1)
plot (CDA.Datetime, CDA.WindSpeed_m_s_)
hold on
plot (Cuchi.Datetime, Cuchi.WindSpeed_m_s_)
legend ('CDA', 'Cuchi', 'orientation','horizontal','location','best')
ylabel({'Wind Speed (m s^{-1})'});
subplot(2,1,2)
plot (CDA.Datetime, CDA.WindSpeed_m_s_)
hold on
plot (Cuchi.Datetime, Cuchi.WindSpeed_m_s_)
xlim([datetime('10-Feb-2018') datetime('15-Feb-2018')])
legend ('CDA', 'Cuchi', 'orientation','horizontal','location','best')
ylabel({'Wind Speed (m s^{-1})'});

figname ='CDA_windspeed';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

T = retime(CDA, Cuchi_time,"fillwithmissing");
CDA_values = T(:, 4).Variables;
Cuchi_values =Cuchi(:, 4).Variables;
CDA_time = T.Datetime;
Cuchi_time = Cuchi.Datetime;

Cuchi_values(Cuchi_values>4.5)=nan

t = CDA_time
d  =[CDA_values, Cuchi_values];


% Find indices where both columns have a value
valid_indices = ~isnan(d(:, 1)) & ~isnan(d(:, 2));
valid_indices = find(valid_indices)

coefficients = polyfitZero(d(valid_indices, 2), d(valid_indices, 1),1)
yfit = polyval(coefficients, Cuchi_values) % this is also looking super weong

r = corrcoef(d(:,1), yfit, 'rows', 'pairwise')
r_corr = r(2);
RMSE = rmse(d(:,1), yfit, 'omitnan')
r = corrcoef(d(:,1), d(:,2), 'rows', 'pairwise')
r_or = r(2);
RMSE_or = rmse(d(:,1), d(:,2), 'omitnan')


figure;
scatter(d(:, 1), d(:,2), 'r.');hold on
scatter(d(:, 1), yfit, 'k.'); 
grid on
lsline
rl = refline(1,0)
rl.Color ='k'
ylabel('Cuchi')
xlabel('CDA')
legend( 'Original','Regression', 'best fit, or','best fit, corr','1:1', 'location', 'best')
grid on
text(0.5, 14, strcat('y = ', num2str(round(coefficients(1), 2)), '+' , strcat(num2str(round(coefficients(2), 2)))))
text (0.5, 13, strcat( 'RMSEor=', num2str(round(RMSE_or,2)), ', RMSEcorr = ', num2str(round(RMSE,2))))
text (0.5, 12, strcat( 'R =', num2str(round(r_or,2)), ', Rcorr = ', num2str(round(r_corr,2))))

figname ='CDA_windspeed_scatterplot';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all


figure
hold on
subplot(2,1,1)
plot (CDA.Datetime, CDA.WindSpeed_m_s_)
hold on
plot (Cuchi.Datetime, Cuchi.WindSpeed_m_s_)
plot (Cuchi.Datetime, yfit)
legend ('CDA', 'Cuchi', 'Regressed','orientation','horizontal','location','best')
ylabel({'Wind Speed (m s^{-1})'});
subplot(2,1,2)
plot (CDA.Datetime, CDA.WindSpeed_m_s_)
hold on
plot (Cuchi.Datetime, Cuchi.WindSpeed_m_s_)
plot (Cuchi.Datetime, yfit)
xlim([datetime('10-Feb-2018') datetime('15-Feb-2018')])
legend ('CDA', 'Cuchi', 'Regressed', 'orientation','horizontal','location','best')
ylabel({'Wind Speed (m s^{-1})'});

figname ='CDA_windspeed_regressed';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

% tehre is a seaosnal poattern at Cuchi dthat we dont get at CDA, so
% instead, I will just rep[eat the precivous year
CDA_values = T(:, 4).Variables;
Cuchi_values =Cuchi(:, 4).Variables;
CDA_time = T.Datetime;
Cuchi_time = Cuchi.Datetime;

a = find(CDA_time == datetime('02-Jul-2019 12:00'))
aa = find(CDA_time == datetime('02-Jul-2015 12:00'))
bb = find(CDA_time == datetime('06-Apr-2016 07:00'))
clear replacement tobereplaced
replacement  = CDA_values(aa:bb);
tobereplaced = CDA_values(a:end);
CDA_values(a:end)=CDA_values(aa:bb)

% plot final results
figure
hold on
plot (Cuchi.Datetime, CDA_values); hold on
plot (CDA.Datetime, CDA.WindSpeed_m_s_)
hold on
plot (Cuchi.Datetime, Cuchi.WindSpeed_m_s_)
legend ('Filled with 2015-2016','CDA', 'Cuchi', 'orientation','horizontal','location','best')
ylabel({'Wind Speed (m s^{-1})'});

figname ='CDA_windspeed_filledwith20152016';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

% export result

a = find(Cuchi_time == '26-Jun-2014')
CDAu_filled = CDA_values(a:end);
CDA_time = Cuchi_time (a:end);

x = timetable(CDA_time, CDAu_filled)
filename = 'F:\11_CRHM_cuchi\data\processed\InfilledU_20130705_20200406.csv'
writetimetable(x, filename);

%% shortwave
%%%%%%%%% 
% Infill 
close all
figure
hold on
subplot(2,1,1)
plot (LlanWx.Datetime, LlanWx.Solar_Wm_2_);hold on
plot (Cuchi.Datetime, Cuchi.Solar_Wm_2_)
plot (CDA.Datetime, CDA.Solar_Wm_2_)
xlim([Cuchi.Datetime(1) Cuchi.Datetime(end)])
legend ('LLanWx', 'Cuchi', 'CDA', 'orientation','horizontal','location','best')
ylabel({'SWin (W m^{-2})'});
subplot(2,1,2)
plot (LlanWx.Datetime, LlanWx.Solar_Wm_2_);hold on
plot (Cuchi.Datetime, Cuchi.Solar_Wm_2_)
plot (CDA.Datetime, CDA.Solar_Wm_2_)
xlim([datetime('10-Sep-2014') datetime('15-Sep-2014')])
legend ('LLanWx', 'Cuchi', 'CDA', 'orientation','horizontal','location','best')
ylabel({'SWin (W m^{-2})'});

figname ='CDA_SW';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

%% Compare SW Cuchi and CDA

T = retime(CDA, Cuchi_time,"fillwithmissing");
CDA_values = T(:, 2).Variables;
Cuchi_values =Cuchi(:, 7).Variables;
CDA_time = T.Datetime;
Cuchi_time = Cuchi.Datetime;
T = retime(LlanWx, Cuchi_time,"fillwithmissing");
Llan_values = T(:, 7).Variables;
Llan_time = T.Datetime;

% Set CDA values that are bad to nan
t1 = '15-Sep-2015'
a = find(Cuchi_time==t1)
CDA_values(a:end)= nan;

%% Regresison 
t = CDA_time
d  =[CDA_values, Cuchi_values];


% Find indices where both columns have a value
valid_indices = ~isnan(d(:, 1)) & ~isnan(d(:, 2));
valid_indices = find(valid_indices)

coefficients = polyfitZero(d(valid_indices, 1),d(valid_indices, 2),1)
yfit = polyval(coefficients, Cuchi_values) % this is also looking super weong

CDA_filled = yfit;

r = corrcoef(d(:,1), yfit, 'rows', 'pairwise')
r_corr = r(2);
RMSE = rmse(d(:,1), yfit, 'omitnan')
r = corrcoef(d(:,1), d(:,2), 'rows', 'pairwise')
r_or = r(2);
RMSE_or = rmse(d(:,1), d(:,2), 'omitnan')


figure;
scatter(d(:, 1), d(:,2), 'r.');hold on
scatter(d(:, 1), yfit, 'k.'); 
grid on
lsline
rl = refline(1,0)
rl.Color ='k'
ylabel('Cuchi')
xlabel('CDA')
legend( 'Original','Regression', 'best fit, or','best fit, corr','1:1', 'location', 'best')
grid on
text(800, 200, strcat('y = ', num2str(round(coefficients(1), 2)), 'x'))
text (800, 150,  strcat( 'RMSEor=', num2str(round(RMSE_or,2)), ', RMSEcorr = ', num2str(round(RMSE,2))))
text (800, 100,  strcat( 'R =', num2str(round(r_or,2)), ', Rcorr = ', num2str(round(r_corr,2))))

figname ='CDA_SW_scatterplot_cuchi';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all


figure
hold on
subplot(2,1,1)
plot (CDA_time, CDA_values)
hold on
plot (Cuchi.Datetime, Cuchi.Solar_Wm_2_)
plot (Cuchi.Datetime, yfit)
legend ('CDA', 'Cuchi', 'Regressed','orientation','horizontal','location','best')
ylabel({'SWin (W m^{-2})'});
subplot(2,1,2)
plot (CDA_time, CDA_values)
hold on
plot (Cuchi.Datetime, Cuchi.Solar_Wm_2_)
plot (Cuchi.Datetime, yfit)
legend ('CDA', 'Cuchi', 'Regressed','orientation','horizontal','location','best')
ylabel({'SWin (W m^{-2})'});
xlim([datetime('10-Oct-2014') datetime('15-Oct-2014')])

figname ='CDA_SW_regressed_cuchi';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all


% Now infil to llanganuco
t = CDA_time
d  =[CDA_filled, Llan_values];


% Find indices where both columns have a value
valid_indices = ~isnan(d(:, 1)) & ~isnan(d(:, 2));
valid_indices = find(valid_indices)

coefficients = polyfitZero(d(valid_indices, 2),d(valid_indices, 1),1)
coefficients = [0.9,0]
yfit = polyval(coefficients, Llan_values) % this is also looking super weong

r = corrcoef(d(:,1), yfit, 'rows', 'pairwise')
r_corr = r(2);
RMSE = rmse(d(:,1), yfit, 'omitnan')
r = corrcoef(d(:,1), d(:,2), 'rows', 'pairwise')
r_or = r(2);
RMSE_or = rmse(d(:,1), d(:,2), 'omitnan')


figure;
scatter(d(:, 1), d(:,2), 'r.');hold on
scatter(d(:, 1), yfit, 'k.'); 
grid on
lsline
rl = refline(1,0)
rl.Color ='k'
ylabel('Llan')
xlabel('CDA filled')
legend( 'Original','Regression', 'best fit, or','best fit, corr','1:1', 'location', 'best')
grid on
text(5, 1750, strcat('y = ', num2str(round(coefficients(1), 2)), 'x'))
text (5, 1650,  strcat( 'RMSEor=', num2str(round(RMSE_or,2)), ', RMSEcorr = ', num2str(round(RMSE,2))))
text (5, 1550,  strcat( 'R =', num2str(round(r_or,2)), ', Rcorr = ', num2str(round(r_corr,2))))

figname ='CDA_SW_scatterplot_llan';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all


figure
hold on
subplot(2,1,1)
plot (CDA_time, CDA_filled)
hold on
plot(CDA_time, yfit)
legend ('CDA filled', 'Llan regressed','orientation','horizontal','location','best')
ylabel({'SWin (W m^{-2})'});
subplot(2,1,2)
plot (CDA_time, CDA_filled)
hold on
plot(CDA_time, yfit)
legend ('CDA filled', 'Llan regressed', 'orientation','horizontal','location','best')
ylabel({'SWin (W m^{-2})'});
xlim([datetime('10-Oct-2014') datetime('15-Oct-2014')])

figname ='CDA_SW_regressed_llan';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all


% export result
CDA_filled(isnan(CDA_filled)) = yfit(isnan(CDA_filled));

figure
plot (CDA_time, CDA_filled)

% one more gap to fill

a = find(CDA_time == datetime('22-Jun-2014 17:00'))
b = find(CDA_time == datetime('26-Jun-2014 9:00'))
aa = find(CDA_time == datetime('18-Jun-2014 17:00'))
bb = find(CDA_time == datetime('22-Jun-2014 9:00'))
clear replacement tobereplaced
replacement  = CDA_filled(aa:bb);
tobereplaced = CDA_filled(a:b);
CDA_filled(a:b)=CDA_filled(aa:bb)
a = find(isnan(CDA_filled))

CDA_time(a(1))
CDA_time(a(end))
a = find(CDA_time == datetime('3-Jul-2017 14:00'))
b = find(CDA_time == datetime('5-Jul-2017 1:00'))
aa = find(CDA_time == datetime('1-Jun-2017 14:00'))
bb = find(CDA_time == datetime('3-Jun-2017 1:00'))
clear replacement tobereplaced
replacement  = CDA_filled(aa:bb);
tobereplaced = CDA_filled(a:b);
CDA_filled(a:b)=CDA_filled(aa:bb)
a = find(isnan(CDA_filled))
CDA_time(a(1))
CDA_time(a(end))

a = find(CDA_time == datetime('3-Jul-2019 14:00'))
b = find(CDA_time == datetime('4-Jul-2019 23:00'))
aa = find(CDA_time == datetime('1-Jul-2019 14:00'))
bb = find(CDA_time == datetime('2-Jul-2019 23:00'))
clear replacement tobereplaced
replacement  = CDA_filled(aa:bb);
tobereplaced = CDA_filled(a:b);
CDA_filled(a:b)=CDA_filled(aa:bb)


a = find(Cuchi_time == '26-Jun-2014')
CDAsw_filled = CDA_filled(a:end);
CDA_time = Cuchi_time (a:end);

figure
plot(CDA_time, CDAsw_filled)
legend ('CDA filled','orientation','horizontal','location','best')
ylabel({'SWin (W m^{-2})'});

figname ='CDA_SW_filled';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all
x = timetable(CDA_time, CDAsw_filled)
filename = 'F:\11_CRHM_cuchi\data\processed\InfilledSW_20130705_20200406.csv'
writetimetable(x, filename);


%% RH
figure
hold on
subplot(2,1,1)
plot(LlanWx.Datetime, LlanWx.RH___); hold on
plot (Cuchi.Datetime, Cuchi.RH___)
plot (CDA.Datetime, CDA.RH___)
xlim([Cuchi.Datetime(1) Cuchi.Datetime(end)])
legend ('LlanWx','Cuchi','CDA',  'orientation','horizontal','location','best')
ylabel({'RH (%)'});
subplot(2,1,2)
plot(LlanWx.Datetime, LlanWx.RH___); hold on
plot (Cuchi.Datetime, Cuchi.RH___)
plot (CDA.Datetime, CDA.RH___)
legend ('LlanWx','Cuchi','CDA',  'orientation','horizontal','location','best')
ylabel({'RH (%)'});
xlim([datetime('10-Sep-2017') datetime('15-Sep-2017')])
figname ='CDA_RH';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%% Compare RH Cuchi and CDA
% Fille values 1:1
T = retime(CDA, Cuchi_time,"fillwithmissing");
CDA_values = T(:, 7).Variables;
Cuchi_values =Cuchi(:, 3).Variables;

T= retime(LlanWx,Cuchi_time,"fillwithmissing");
Llan_values = T(:,6).Variables;
CDA_time = T.Datetime;
Cuchi_time = Cuchi.Datetime;

figure
plot(Cuchi_time, CDA_values); hold on
plot(Cuchi_time,Cuchi_values)
plot(Cuchi_time,Llan_values)



CDA_filled =  CDA_values;
% export result
CDA_filled(isnan(CDA_filled)) = Cuchi_values(isnan(CDA_filled));
CDA_filled(isnan(CDA_filled)) = Llan_values(isnan(CDA_filled));
a = find(isnan(CDA_filled))


a = find(Cuchi_time == '26-Jun-2014')
CDArh_filled = CDA_filled(a:end);
CDA_time = Cuchi_time (a:end);


figure
plot(CDA_time, CDArh_filled); hold on
ylabel('RH (%)')
figname ='CDA_RH_filled';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

x = timetable(CDA_time, CDArh_filled)
filename = 'F:\11_CRHM_cuchi\data\processed\InfilledRH_20130705_20200406.csv'
writetimetable(x, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% Collate all in an obs file
T = [CDAta_filled, CDArh_filled, CDAsw_filled ,CDAu_filled,CDAp_filled];
T = timetable(CDA_time, T);
T = retime(T, CDA_time, 'linear');
tt = T.CDA_time;
xx = table2array(T);

t = datevec(T.CDA_time);
T = table2array(T);

% Add a spin up year
% add year 2
Tadd = T(8761:17520, :);
T2 = [Tadd; T];

Tadd = T(1:8760, :);
T2 = [Tadd; T];

tadd = t(1:8760, :);
tadd (:,1)= tadd(:,1)-1;
t2 = [tadd; t]

Tt = timetable(datetime(t2), T2);
Tt = retime(Tt, datetime(t2), 'linear')

t = datevec(Tt.Time);
T = table2array(Tt);

    % create the obs matrix and save it in a mat file
t = t(:, 1:5);

obs = [t T]; % compiled time and data together, and save it in a matlab format
obs = obs(2:end, :);

fn = strcat('data\processed\CuchiObs_2014_2020.mat');
save (fn, 'obs');  
 % create the obs text file
headerlines = {'Obs file, WRF PGW';
              't 1 (C)';
              'rh 1 (hpa)';
              'Qsi 1 (W/m2)';
              'u 1 (m/s)';
              'p 1 (mm)';
              '$$ Missing ' ;
              '#####      t.1  rh.1   qsi.1 u.1    p.1'}
fp = strcat('data\processed\CuchiObs_2014_2020.obs');
fid = fopen(fp, 'wt');
for l = 1:numel(headerlines)
   fprintf(fid, '%s\n', headerlines{l});
end
fclose(fid);
dlmwrite(fp, obs , '-append', 'delimiter', '\t'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare steramflow obs
d = readtable('F:\11_CRHM_cuchi\data\Q\cda_lev-q-t.csv'); % Read the table from the file
d = table2timetable(d)
d = retime(d, 'daily', 'mean')
Q = d.Discharge_cms_* 3600;
Q(isnan(Q))= 0;

    % create the obs matrix and save it in a mat file
t = datevec(d.Datetime);
t = t(:, 1:5);
t(:,4)=1;

Qobs= [t Q]; % compiled time and data together, and save it in a matlab format


fn = strcat('F:\11_CRHM_cuchi\data\processed\CuchiQ_2009_2019_hrly.mat');
save (fn, 'Qobs');  
 % create the obs text file
headerlines = {'Obs file, WRF PGW';
              'a 1 (C)';
              '$$ Missing ' ;
              '#####      a.1 '}
fp = strcat('F:\11_CRHM_cuchi\data\processed\CuchiQ_2009_2019_hrly.obs');
fid = fopen(fp, 'wt');
for l = 1:numel(headerlines)
   fprintf(fid, '%s\n', headerlines{l});
end
fclose(fid);
dlmwrite(fp, Qobs , '-append', 'delimiter', '\t'); 

%% Prepare steramflow obs
d = readtable('F:\11_CRHM_cuchi\data\Q\cda_lev-q-t.csv'); % Read the table from the file
d = table2timetable(d)
d = retime(d, 'daily', 'mean')
Q = d.Discharge_cms_* 86400;
Q(isnan(Q))= 0;

    % create the obs matrix and save it in a mat file
t = datevec(d.Datetime);
t = t(:, 1:5);
t(:,4)=1;

Qobs= [t Q]; % compiled time and data together, and save it in a matlab format


fn = strcat('F:\11_CRHM_cuchi\data\processed\CuchiQ_2009_2019.mat');
save (fn, 'Qobs');  
 % create the obs text file
headerlines = {'Obs file, WRF PGW';
              'b 1 (C)';
              '$$ Missing ' ;
              '#####      b.1 '}
fp = strcat('F:\11_CRHM_cuchi\data\processed\CDAQ_2009_2019.obs');
fid = fopen(fp, 'wt');
for l = 1:numel(headerlines)
   fprintf(fid, '%s\n', headerlines{l});
end
fclose(fid);
dlmwrite(fp, Qobs , '-append', 'delimiter', '\t'); 

