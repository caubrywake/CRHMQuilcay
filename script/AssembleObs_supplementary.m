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
figdir = 'G:\11_CRHM_cuchi\fig\obs\v2\'
folderPath = 'G:\11_CRHM_cuchi\data\Wx\';  % find the met files
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

% Adjust the size of the figure
figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]

% First subplot
subplot(2,1,1)
plot(CDA.Datetime, CDA.Pcorr);
ylabel({'Hourly Precipitation,'; 'Corrected (mm)'})
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)
hold on

% Adjust x-axis ticks to every 6 months
xt = datetime(min(CDA.Datetime), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(CDA.Datetime), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
hold on

% Second subplot
subplot(2,1,2)
plot(CDA.Datetime, cumsum(CDA.Precipitation_mm_)); hold on
plot(CDA.Datetime, cumsum(CDA.Pcorr));
ylabel('Cumul. Precipitation (mm)')
legend('Original', 'Corrected', 'location', 'southeast')
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)
hold on

% Adjust x-axis ticks to every 6 months
xt = datetime(min(CDA.Datetime), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(CDA.Datetime), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
hold on

% Save the figure
figname = 'F20_CDA_precipitation_hourly_cumul';
saveas(gcf, fullfile(figdir, strcat(figname, '.pdf')))
saveas(gcf, fullfile(figdir, strcat(figname, '.png')))
savefig(gcf, fullfile(figdir, figname))

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

figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]

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
grid on; box on
lsline
rl = refline(1,0)
rl.Color ='k'
ylabel('Cuchillacocha')
xlabel('CDA')
legend( 'Original','Corrected', 'best fit, or','best fit, corr','1:1', 'location', 'best')
text(-4, 22, strcat('y = ', num2str(round(coefficients(1), 2)), '+' , strcat(num2str(round(coefficients(2), 2)))))
text (-4, 20, strcat( 'RMSEor=', num2str(round(RMSE_or,2)), ', RMSEcorr = ', num2str(round(RMSE,2))))
text (-4, 18, strcat( 'R =', num2str(round(r_or,2)), ', Rcorr = ', num2str(round(r_corr,2))))

grid on; box on
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


figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
plot(CDA_time, CDAp_filled); hold on
plot(CDA_time, CDA_filled_0)
legend ('Filled from Cuchillacocha','CDA','orientation','horizontal','location','best')
ylabel({'Hourly Precipitation (mm)'});
ylim([0 12])
xt = datetime(min(CDA_time), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(CDA_time), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
hold on

figname ='F21_CDA_precip_Filled';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%% Elevation gradient for precipitation
% Elevation of the stations
elevation1 = 3924; % meters
elevation2 = 4642; % meters

% Calculate the elevation difference
elevation_diff = elevation2 - elevation1;

% Convert the data to a timetable
TT = timetable(t, d(:, 1), d(:, 2), 'VariableNames', {'Station1', 'Station2'});

% Calculate the monthly sums
monthly_sum = retime(TT, 'monthly', 'sum');

% Calculate the monthly precipitation difference
monthly_diff = monthly_sum.Station2 - monthly_sum.Station1;


% Calculate the precipitation gradient


monthly_gradient = monthly_diff / elevation_diff;
a = find(monthly_sum.Station2==0);
monthly_gradient(a) = nan;
a = find(monthly_sum.Station1==0);
monthly_gradient(a) = nan;
a = find(monthly_gradient<=0);
monthly_gradient(a) = nan;
% Add the precipitation gradient to the timetable
monthly_gradient_table = timetable(monthly_sum.t, monthly_gradient, ...
                                   'VariableNames', {'PrecipitationGradient'});

pfactor = 0.059;
monthly_sum.Station2corr = monthly_sum.Station1*(1+elevation_diff *pfactor/100);

figure('Position', [100, 100, 800, 400]);

% Subplot 1: Bar plot for monthly precipitation sums
subplot(2, 1, 1);
plot(monthly_gradient_table.Time, monthly_gradient_table.PrecipitationGradient, '-xk');
grid on;
box on;
ylabel({'Precipitation Gradient'; '(mm/m)'});
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)

% Set x-ticks to every 6 months
xt = datetime(min(monthly_gradient_table.Time), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(monthly_gradient_table.Time), 'Format', 'MMM-yyyy');
xticks(xt);
xtickformat('MMM-yyyy');

% Subplot 2: Line plot for monthly precipitation gradient
subplot(2, 1, 2);
bar(monthly_sum.t, [monthly_sum.Station1, monthly_sum.Station2, monthly_sum.Station2corr]);
ylabel({'Monthly Precipitation Sum';'(mm)'});
legend({'CDA', 'Cuchillacocha', 'CDA with elevation correction to Cuchillacocha'}, 'Location', 'best', 'orientation', 'horizontal');
grid on;
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)

% Set x-ticks to every 6 months
xticks(xt);
xtickformat('MMM-yyyy');


% Save the figure
figname = 'Monthly_Precipitation_Sums_and_Gradient';
saveas(gcf, fullfile(figdir, strcat(figname, '.pdf')));
saveas(gcf, fullfile(figdir, strcat(figname, '.png')));
savefig(gcf, fullfile(figdir, figname));



nanmean(monthly_gradient_table.PrecipitationGradient)
%%
% calculate annual sum
Tt = timetable(CDA_time, CDAp_filled)
T = retime(Tt, 'yearly', 'sum')

figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
subplot(2,1,1)

X = T.CDA_time(2:end-1);
Y = T.CDAp_filled(2:end-1);
bar(X, Y)
labels = arrayfun(@(value) num2str(round(value)), Y, 'UniformOutput', false);
text(X,Y,labels,'HorizontalAlignment','center','VerticalAlignment','bottom') 
      % clears X axis data
ylabel('Annual Precipitation (mm)' )
ylim([0 1100])
xlim ([X(1)-calmonths(6) X(end)+calmonths(6)])
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)


subplot(2,1,2)
Tt = timetable(CDA_time, CDAp_filled)
T = retime(Tt, 'monthly', 'sum')
X = T.CDA_time;
Y = T.CDAp_filled;
bar(X, Y)
xt = datetime(min(X), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(X), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
ylabel('Monthly Precipitation (mm)' )
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)
ylim([0 300])
figname ='F22_CDA_annualPrecip_bar';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

x = timetable(CDA_time,  CDAp_filled)
filename = 'G:\11_CRHM_cuchi\data\processed\InfilledP_20130705_20200406.csv'
writetimetable(x, filename);

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Air temperature
figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
hold on
subplot(2,1,1)
plot (LlanWx.Datetime, LlanWx.Temperature__C_); hold on
plot (LlanPort.Datetime, LlanPort.AirTemperature__C_)
plot (Cuchi.Datetime, Cuchi.Temperature__C_)
plot (CDA.Datetime, CDA.Temperature__C_)
ylim([-5 23])
legend ('LlanWx','LlanPort','Cuchi','CDA', 'orientation','horizontal','location','northoutside')
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)
xt = datetime(min(LlanWx.Datetime), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(Cuchi.Datetime), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
ylabel('Hourly Air Temperature (°C)')

subplot(2,1,2)
plot (LlanWx.Datetime, LlanWx.Temperature__C_); hold on
plot (LlanPort.Datetime, LlanPort.AirTemperature__C_)
plot (Cuchi.Datetime, Cuchi.Temperature__C_)
plot (CDA.Datetime, CDA.Temperature__C_)
ylim([-4 20])
xlim([datetime('10-Feb-2018') datetime('15-Feb-2018')])
legend ('LlanWx','LlanPort','Cuchi','CDA', 'orientation','horizontal','location','best')
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)
ylabel('Hourly Air Temperature (°C)')

figname ='F2_CDA_temperature_allstations_hourly_2006_2020';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
hold on
subplot(2,1,1)
plot(t, yfit); hold on; plot(t, d(:,1)); plot(t, d(:,2))
ylim([-5 23])
legend ('Cuchillacocha - corrected to CDA','CDA','Cuchillacocha', 'orientation','horizontal','location','best')
ylabel({'Hourly Air Temperature ({\circ}C)'});
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)
xt = datetime(min(Cuchi.Datetime), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(Cuchi.Datetime), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')

subplot(2,1,2)
plot(t, yfit); hold on; plot(t, d(:,1)); plot(t, d(:,2))
ylim([-4 20])
xlim([datetime('10-Feb-2018') datetime('15-Feb-2018')])
legend ('Cuchillacocha - corrected to CDA','CDA','Cuchillacocha', 'orientation','horizontal','location','best')
ylabel({'Hourly Air Temperature ({\circ}C)'});
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)
ylabel('Hourly Air Temperature (°C)')

figname ='F6_CDA_temperature_RegressedCuchi_2006_2020';
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

close all
figure;
hold on

% Scatter plots
scatter(d(:, 1), d(:,2), 'r.'); 
scatter(d(:, 1), yfit, 'k.');

% Add lsline for each scatter plot
hls1 = lsline;
hls1(1).Color = 'k';
hls1(2).Color = 'r';

% Add refline with dashed grey color
rl = refline(1,0);
rl.Color = [0.5, 0.5, 0.5]; % grey color
rl.LineStyle = '--'; % dashed line

% Labels and legend
ylabel('Cuchillacocha Air Temperature ({\circ}C)');
xlabel('CDA Air Temperature ({\circ}C)');
hleg = legend( 'CDA vs Cuchi., original', 'CDA vs Cuchi., after correction', 'best fit, corrected', 'best fit, original','1:1', 'location', 'southeast');
%set(hleg, 'Box', 'off', 'Color', 'none'); % Remove box and set transparent background

% Add text annotations
text(-4, 22, ['Cuchi. corrected = ' num2str(round(coefficients(1), 2)) ' * CDA + ' num2str(round(coefficients(2), 2))], 'FontSize', 10);
text(-4, 21, ['RMSE (CDA and Cuchi. original) = ' num2str(round(RMSE_or, 2)), '{\circ}C'], 'FontSize', 10);
text(-4, 20, ['RMSE (CDA and Cuchi. corrected) = ' num2str(round(RMSE, 2)), '{\circ}C'], 'FontSize', 10);
text(-4, 19, ['R (CDA and Cuchi. original) = ' num2str(round(r_or, 2))], 'FontSize', 10);
text(-4, 18, ['R (CDA and Cuchi. corrected) = ' num2str(round(r_corr, 2))], 'FontSize', 10);

% Grid
grid on; box on

figname = 'F5_CDAvsCuchi_linearregression_airtemp'
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


llan_filled = yfit; 

figure;
hold on

% Scatter plots
scatter(d(:, 1), d(:,2), 'r.'); 
scatter(d(:, 1), yfit, 'k.');


% Add lsline for each scatter plot
hls1 = lsline;
hls1(1).Color = 'k';
hls1(2).Color = 'r';

% Add refline with dashed grey color
rl = refline(1,0);
rl.Color = [0.5, 0.5, 0.5]; % grey color
rl.LineStyle = '--'; % dashed line

% Labels and legend
ylabel('Llanganuco Wx Air Temperature ({\circ}C)');
xlabel('CDA Air Temperature ({\circ}C)');
hleg = legend('CDA vs LlwanWx., original','CDA vs LlanWx., after correction', 'best fit, corrected',  'best fit, original', '1:1', 'location', 'southeast');
%set(hleg, 'Box', 'off', 'Color', 'none'); % Remove box and set transparent background

% Add text annotations
text(-4, 24, ['LlanWx corrected = ' num2str(round(coefficients(1), 2)) ' * CDA + ' num2str(round(coefficients(2), 2))], 'FontSize', 10);
text(-4, 23, ['RMSE (CDA and LlanWx original) = ' num2str(round(RMSE_or, 2)), '{\circ}C'], 'FontSize', 10);
text(-4, 22, ['RMSE (CDA and LlanWx corrected) = ' num2str(round(RMSE, 2)), '{\circ}C'], 'FontSize', 10);
text(-4, 21, ['R (CDA and LlanWx original) = ' num2str(round(r_or, 2))], 'FontSize', 10);
text(-4, 20, ['R (CDA and LlanWx corrected) = ' num2str(round(r_corr, 2))], 'FontSize', 10);

% Grid
grid on; box on

figname ='F3_CDA_temperature_RegressedLlan_2006_2020';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))



close all

figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
hold on
subplot(2,1,1)
plot(t, d(:,1));hold on; plot(t, d(:,2)); plot(t, yfit);
ylim([-5 23])
legend ('CDA','Llan','LlanWx- corrected to CDA', 'orientation','horizontal','location','best')
ylabel({'Hourly Air Temperature ({\circ}C)'});
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)
xt = datetime(min(t), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(t), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')

subplot(2,1,2)
plot(t, d(:,1));hold on; plot(t, d(:,2)); plot(t, yfit);
ylim([-4 20])
xlim([datetime('10-Feb-2018') datetime('15-Feb-2018')])
legend ('CDA','Llan','LlanWx- corrected to CDA', 'orientation','horizontal','location','best')
ylabel({'Hourly Air Temperature ({\circ}C)'});
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)
ylabel('Hourly Air Temperature (°C)')

figname ='F4_CDA_temperature_RegressedLlan_2006_2020';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

%% Infill air temp at CDA with llan thn Cuchi
% Infill air temp at CDA with llan then Cuchi
CDA_filled = retime(CDA, Cuchi.Datetime, 'fillwithmissing');
CDA_values = CDA_filled(:, 6).Variables;
CDA_filled = CDA_filled(:, 6).Variables;
CDA_filled(isnan(CDA_filled)) = llan_filled(isnan(CDA_filled));
CDA_filled(isnan(CDA_filled)) = cuchi_filled(isnan(CDA_filled));

a = find(Cuchi_time == '26-Jun-2014');
CDAta_filled = CDA_filled(a:end);
CDA_time = Cuchi_time(a:end);

% Calculate daily average, min, and max temperatures
CDA_table = timetable(CDA_time, CDAta_filled);
daily_stats = retime(CDA_table, 'daily', @mean);
daily_min = retime(CDA_table, 'daily', @min);
daily_max = retime(CDA_table, 'daily', @max);

% Create figure with 2 subplots
figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]

% First subplot: Infilled temperature
subplot(2, 1, 1);
plot(CDA_time, CDAta_filled);
legend('Filled CDA', 'Orientation', 'horizontal', 'Location', 'best');
ylabel('Hourly Air Temperature ({\circ}C)');
grid on; box on;
xt = datetime(min(CDA_time), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(CDA_time), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)

% Second subplot: Daily average, min, and max temperatures
subplot(2, 1, 2);
hold on;
plot(daily_stats.CDA_time, daily_stats.CDAta_filled, 'DisplayName', 'Average');
plot(daily_min.CDA_time, daily_min.CDAta_filled, 'DisplayName', 'Min');
plot(daily_max.CDA_time, daily_max.CDAta_filled,  'DisplayName', 'Max');
legend('show', 'Location', 'best', 'Orientation','horizontal');
ylabel('Daily Average Temperature ({\circ}C)');
grid on; box on;
xt = datetime(min(CDA_time), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(CDA_time), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)

% Save the figure
figname = 'F7_CDA_temperature_Filled';
saveas(gcf, fullfile(figdir, strcat(figname, '.pdf')));
saveas(gcf, fullfile(figdir, strcat(figname, '.png')));
savefig(gcf, fullfile(figdir, figname));
close all;

% export that time series:
x = timetable(CDA_time, CDAta_filled)
filename = 'G:\11_CRHM_cuchi\data\processed\InfilledTa_20130705_20200406.csv'
writetimetable(x, filename);

t_obs = CDA_time;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Wind speed
%% Infill u
 
figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
hold on
subplot(2,1,1)
plot (CDA.Datetime, CDA.WindSpeed_m_s_)
hold on
plot (Cuchi.Datetime, Cuchi.WindSpeed_m_s_)
legend ('CDA', 'Cuchillacocha', 'orientation','horizontal','location','best')
ylabel({'Hourly Wind Speed (m s^{-1})'});
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)
hold on

% Adjust x-axis ticks to every 6 months
xt = datetime(min(CDA.Datetime), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(Cuchi.Datetime), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
hold on

subplot(2,1,2)
plot (CDA.Datetime, CDA.WindSpeed_m_s_)
hold on
plot (Cuchi.Datetime, Cuchi.WindSpeed_m_s_)
xlim([datetime('10-Feb-2018') datetime('15-Feb-2018')])
legend ('CDA', 'Cuchillacocha', 'orientation','horizontal','location','best')
ylabel({'Hourly Wind Speed (m s^{-1})'});
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)
hold on


figname ='F16_CDA_windspeed';
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
hold on

% Scatter plots
scatter(d(:, 1), d(:,2), 'r.'); 
scatter(d(:, 1), yfit, 'k.');


% Add lsline for each scatter plot
hls1 = lsline;

hls1(2).Color = 'r';
hls1(1).Color = 'k';
% Add refline with dashed grey color
rl = refline(1,0);
rl.Color = [0.5, 0.5, 0.5]; % grey color
rl.LineStyle = '--'; % dashed line

% Labels and legend
ylabel('Cuchillacocha Wind Speed (m s^{-1})');
xlabel('CDA Wind Speed (m s^{-1})');
hleg = legend('CDA vs Cuchi., original','CDA vs Cuchi., after correction',   'best fit, corrected', 'best fit, original','1:1', 'location', 'southeast');
%set(hleg, 'Box', 'off', 'Color', 'none'); % Remove box and set transparent background
xlim ([0 9])
ylim([0 9])

% Add text annotations
text(0.2, 8.8, ['Cuchi. corrected = ' num2str(round(coefficients(1), 2)) ' * CDA + ' num2str(round(coefficients(2), 2))], 'FontSize', 10);
text(0.2, 8.4, ['RMSE (CDA and Cuchi. original) = ' num2str(round(RMSE_or, 2)), 'm s^{-1}'], 'FontSize', 10);
text(0.2, 8, ['RMSE (CDA and Cuchi. corrected) = ' num2str(round(RMSE, 2)), 'm s^{-1}'], 'FontSize', 10);
text(0.2, 7.6, ['R (CDA and Cuchi. original) = ' num2str(round(r_or, 2))], 'FontSize', 10);
text(0.2, 7.2, ['R (CDA and Cuchi. corrected) = ' num2str(round(r_corr, 2))], 'FontSize', 10);

% Grid
grid on; box on
figname ='F17_CDA_windspeed_scatterplot';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all




%%
close all
figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
hold on
subplot(2,1,1)
plot (CDA.Datetime, CDA.WindSpeed_m_s_)
hold on
plot (Cuchi.Datetime, yfit)
plot (Cuchi.Datetime, Cuchi.WindSpeed_m_s_)

legend ('CDA', 'Cuchillacocha - corrected to CDA','Cuchillacocha', 'orientation','horizontal','location','best')
ylabel({'Hourly Wind Speed (m s^{-1})'});
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)
hold on

% Adjust x-axis ticks to every 6 months
xt = datetime(min(CDA.Datetime), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(Cuchi.Datetime), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
hold on


subplot(2,1,2)
plot (CDA.Datetime, CDA.WindSpeed_m_s_)
hold on
plot (Cuchi.Datetime, yfit)
plot (Cuchi.Datetime, Cuchi.WindSpeed_m_s_)
legend ('CDA', 'Cuchillacocha - corrected to CDA','Cuchillacocha', 'orientation','horizontal','location','best')
xlim([datetime('10-Feb-2018') datetime('15-Feb-2018')])
ylabel({'Hourly Wind Speed (m s^{-1})'});
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)


figname ='F18_CDA_windspeed_regressed';
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
figure('Position', [100, 100, 800, 300]); % Adjust the figure size if needed
hold on
plot (Cuchi.Datetime, CDA_values); hold on
plot (CDA.Datetime, CDA.WindSpeed_m_s_)
legend ('Filled with 2015-2016','CDA', 'orientation','horizontal','location','best')
ylabel({'Hourly Wind Speed (m s^{-1})'});
xt = datetime(min(CDA.Datetime), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(Cuchi.Datetime), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
box on
figname ='F19_CDA_windspeed_filledwith20152016';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

% export result

a = find(Cuchi_time == '26-Jun-2014')
CDAu_filled = CDA_values(a:end);
CDA_time = Cuchi_time (a:end);

x = timetable(CDA_time, CDAu_filled)
filename = 'G:\11_CRHM_cuchi\data\processed\InfilledU_20130705_20200406.csv'
writetimetable(x, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shortwave


%%%%%%%%% 
% Infill 
close all
figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
hold on
subplot(2,1,1)
plot (LlanWx.Datetime, LlanWx.Solar_Wm_2_);hold on
plot (Cuchi.Datetime, Cuchi.Solar_Wm_2_)
plot (CDA.Datetime, CDA.Solar_Wm_2_)
xlim([Cuchi.Datetime(1) Cuchi.Datetime(end)])
legend ('LLanWx', 'Cuchillacocha', 'CDA', 'orientation','horizontal','location','best')
ylabel({'Hourly SWin (W m^{-2})'});
box on;
xt = datetime(min(CDA.Datetime), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(Cuchi.Datetime), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)

subplot(2,1,2)
plot (LlanWx.Datetime, LlanWx.Solar_Wm_2_);hold on
plot (Cuchi.Datetime, Cuchi.Solar_Wm_2_)
plot (CDA.Datetime, CDA.Solar_Wm_2_)
xlim([datetime('10-Sep-2014') datetime('15-Sep-2014')])
legend ('LLanWx', 'Cuchillacocha', 'CDA', 'orientation','horizontal','location','best')
ylabel({'Hourly SWin (W m^{-2})'});
box on; 
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)

figname ='F10_CDA_SW2';
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


close all

figure;
hold on

% Scatter plots
scatter(d(:, 1), d(:,2), 'r.'); 
scatter(d(:, 1), yfit, 'k.');


% Add lsline for each scatter plot
hls1 = lsline;

hls1(2).Color = 'r';
hls1(1).Color = 'k';
% Add refline with dashed grey color
rl = refline(1,0);
rl.Color = [0.5, 0.5, 0.5]; % grey color
rl.LineStyle = '--'; % dashed line

% Labels and legend
ylabel('Cuchillacocha SWin  (W m^{-2})');
xlabel('CDA SWin  (W m^{-2})');
hleg = legend('CDA vs Cuchi., original','CDA vs Cuchi., after correction',   'best fit, corrected', 'best fit, original','1:1', 'location', 'southeast');
%set(hleg, 'Box', 'off', 'Color', 'none'); % Remove box and set transparent background

% Add text annotations
text(10, 1350, ['Cuchi. corrected = ' num2str(round(coefficients(1), 2)) ' * CDA + ' num2str(round(coefficients(2), 2))], 'FontSize', 10);
text(10, 1300, ['RMSE (CDA and Cuchi. original) = ' num2str(round(RMSE_or)), ' W m^{-2}'], 'FontSize', 10);
text(10, 1250, ['RMSE (CDA and Cuchi. corrected) = ' num2str(round(RMSE)), ' W m^{-2}'], 'FontSize', 10);
text(10, 1200, ['R (CDA and Cuchi. original) = ' num2str(round(r_or, 2))], 'FontSize', 10);
text(10, 1150, ['R (CDA and Cuchi. corrected) = ' num2str(round(r_corr, 2))], 'FontSize', 10);

% Grid
grid on; box on
figname ='F11_CDA_SW_scatterplot_cuchi';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

%% Time series infilled

figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
hold on
subplot(2,1,1)
plot (CDA_time, CDA_values)
hold on
plot (Cuchi.Datetime, Cuchi.Solar_Wm_2_)
plot (Cuchi.Datetime, yfit)
legend ('CDA', 'Cuchi', 'Cuchi. - corrected to CDA','orientation','horizontal','location','best')
ylabel({'Hourly SWin (W m^{-2})'});
box on
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)

% Adjust x-axis ticks to every 6 months
xt = datetime(min(CDA.Datetime), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(Cuchi.Datetime), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
hold on


subplot(2,1,2)
plot (CDA_time, CDA_values)
hold on
plot (Cuchi.Datetime, Cuchi.Solar_Wm_2_)
plot (Cuchi.Datetime, yfit)
legend ('CDA', 'Cuchi', 'Cuchi. - corrected to CDA','orientation','horizontal','location','best')
ylabel({'Hourly SWin (W m^{-2})'});
xlim([datetime('10-Oct-2014') datetime('15-Oct-2014')])
box on
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)

figname ='F12_CDA_SW_regressed_cuchi';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Scatter plots
scatter(d(:, 1), d(:,2), 'r.'); hold on
scatter(d(:, 1), yfit, 'k.');


% Add lsline for each scatter plot
hls1 = lsline;

hls1(2).Color = 'r';
hls1(1).Color = 'k';
% Add refline with dashed grey color
rl = refline(1,0);
rl.Color = [0.5, 0.5, 0.5]; % grey color
rl.LineStyle = '--'; % dashed line
xlim([0 1800]) 
ylim([0 1800])

% Labels and legend
ylabel('LlanWxSWin  (W m^{-2})');
xlabel('CDA SWin  (W m^{-2})');
hleg = legend('CDA vs Llan., original','CDA vs Llan., after correction',   'best fit, corrected', 'best fit, original','1:1', 'location', 'southeast');
%set(hleg, 'Box', 'off', 'Color', 'none'); % Remove box and set transparent background

% Add text annotations
text(10, 1750, ['Llan. corrected = ' num2str(round(coefficients(1), 2)) ' * CDA + ' num2str(round(coefficients(2), 2))], 'FontSize', 10);
text(10, 1680, ['RMSE (CDA and Llan. original) = ' num2str(round(RMSE_or)), ' W m^{-2}'], 'FontSize', 10);
text(10, 1610, ['RMSE (CDA and Llan. corrected) = ' num2str(round(RMSE)), ' W m^{-2}'], 'FontSize', 10);
text(10, 1540, ['R (CDA and Llan. original) = ' num2str(round(r_or, 2))], 'FontSize', 10);
text(10, 1470, ['R (CDA and Llan. corrected) = ' num2str(round(r_corr, 2))], 'FontSize', 10);

% Grid
grid on; box on

figname ='F13_CDA_SW_scatterplot_llan';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all


figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
hold on
subplot(2,1,1)
plot (CDA_time, CDA_filled)
hold on
plot(CDA_time, yfit)
legend ('CDA', 'LlanWx- corrected to CDA','orientation','horizontal','location','north')
ylabel({'Hourly SWin (W m^{-2})'});
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)
hold on

% Adjust x-axis ticks to every 6 months
xt = datetime(min(CDA.Datetime), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(CDA.Datetime), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
hold on


subplot(2,1,2)
plot (CDA_time, CDA_filled)
hold on
plot(CDA_time, yfit)
legend ('CDA', 'LlanWx- corrected to CDA', 'orientation','horizontal','location','best')
ylabel({'SWin (W m^{-2})'});
xlim([datetime('10-Oct-2014') datetime('15-Oct-2014')])
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)

figname ='F14_CDA_SW_regressed_llan';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all


% export result
CDA_filled(isnan(CDA_filled)) = yfit(isnan(CDA_filled));

figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
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

figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
plot(CDA_time, CDAsw_filled)
legend ('CDA filled','orientation','horizontal','location','best')
ylabel({'SWin (W m^{-2})'});


% Calculate daily maximum SWin
CDA_table = timetable(CDA_time, CDAsw_filled);
daily_max_SWin = retime(CDA_table, 'daily', @max);

% Create figure with 2 subplots
figure('Position', [100, 100, 800, 300]); % Adjust the figure size if needed

% Subplot 1: Plot of filled CDA SWin
plot(CDA_time, CDAsw_filled);
legend('CDA filled', 'Orientation', 'horizontal', 'Location', 'best');
ylabel('Hourly SWin (W m^{-2})');
grid on;
% Adjust x-axis ticks to every 6 months
xt = datetime(min(CDA_time), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(CDA_time), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
hold on


figname ='F15_CDA_SW_filled';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))
close all


x = timetable(CDA_time, CDAsw_filled)
filename = 'G:\11_CRHM_cuchi\data\processed\InfilledSW_20130705_20200406.csv'
writetimetable(x, filename);


%% RH
figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]

hold on
subplot(2,1,1)
plot(LlanWx.Datetime, LlanWx.RH___); hold on
plot (Cuchi.Datetime, Cuchi.RH___)
plot (CDA.Datetime, CDA.RH___)
xlim([Cuchi.Datetime(1) Cuchi.Datetime(end)])
legend ('LlanWx','Cuchillacocha','CDA',  'orientation','horizontal','location','best')
ylabel({'Hourly RH (%)'});
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)
hold on

% Adjust x-axis ticks to every 6 months
xt = datetime(min(CDA.Datetime), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(Cuchi.Datetime), 'Format', 'MMM-yyyy');
xticks(xt)
datetick('x', 'mmm-yyyy', 'keepticks')
hold on

subplot(2,1,2)
plot(LlanWx.Datetime, LlanWx.RH___); hold on
plot (Cuchi.Datetime, Cuchi.RH___)
plot (CDA.Datetime, CDA.RH___)
legend ('LlanWx','Cuchillacocha','CDA',  'orientation','horizontal','location','best')
ylabel({'Hourly RH (%)'});
xlim([datetime('10-Sep-2017') datetime('15-Sep-2017')])
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)


figname ='F8_CDA_RH2';
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

figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]
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


% Calculate daily statistics for RH
CDA_table = timetable(CDA_time, CDArh_filled);
daily_avg_RH = retime(CDA_table, 'daily', @mean);
daily_min_RH = retime(CDA_table, 'daily', @min);
daily_max_RH = retime(CDA_table, 'daily', @max);

% Create figure with 2 subplots
figure('Position', [100, 100, 800, 400]); % [left, bottom, width, height]

% Subplot 1: Plot of filled CDA RH
subplot(2, 1, 1);
plot(CDA_time, CDArh_filled);
ylabel('Hourly RH (%)');
grid on;
hold on;
text(0.01, 0.9, '(a)', 'Units', 'normalized', 'FontSize', 10)

% Adjust x-axis ticks to every 6 months
xt = datetime(min(CDA_time), 'Format', 'MMM-yyyy'):calmonths(6):datetime(max(CDA_time), 'Format', 'MMM-yyyy');
xticks(xt);
datetick('x', 'mmm-yyyy', 'keepticks');

% Subplot 2: Plot of daily average, min, max RH
subplot(2, 1, 2);
hold on;
plot(daily_avg_RH.CDA_time, daily_avg_RH.CDArh_filled, 'DisplayName', 'Avg');
plot(daily_min_RH.CDA_time, daily_min_RH.CDArh_filled, 'DisplayName', 'Min');
plot(daily_max_RH.CDA_time, daily_max_RH.CDArh_filled, 'DisplayName', 'Max');
legend('show', 'Location', 'best', 'orientation','horizontal');
ylabel('Daily RH (%)');
text(0.01, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 10)
box on
grid on;

% Save the figure
figname = 'F9_CDA_RH_filled';
saveas(gcf, fullfile(figdir, strcat(figname, '.pdf')));
saveas(gcf, fullfile(figdir, strcat(figname, '.png')));
savefig(gcf, fullfile(figdir, figname));
close all;



x = timetable(CDA_time, CDArh_filled)
filename = 'G:\11_CRHM_cuchi\data\processed\InfilledRH_20130705_20200406.csv'
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare steramflow obs
d = readtable('G:\11_CRHM_cuchi\data\Q\cda_lev-q-t.csv'); % Read the table from the file
d = table2timetable(d)
d = retime(d, 'daily', 'mean')
Q = d.Discharge_cms_* 3600;
Q(isnan(Q))= 0;

    % create the obs matrix and save it in a mat file
t = datevec(d.Datetime);
t = t(:, 1:5);
t(:,4)=1;

Qobs= [t Q]; % compiled time and data together, and save it in a matlab format


fn = strcat('G:\11_CRHM_cuchi\data\processed\CuchiQ_2009_2019_hrly.mat');
save (fn, 'Qobs');  
 % create the obs text file
headerlines = {'Obs file, WRF PGW';
              'a 1 (C)';
              '$$ Missing ' ;
              '#####      a.1 '}
fp = strcat('G:\11_CRHM_cuchi\data\processed\CuchiQ_2009_2019_hrly.obs');
fid = fopen(fp, 'wt');
for l = 1:numel(headerlines)
   fprintf(fid, '%s\n', headerlines{l});
end
fclose(fid);
dlmwrite(fp, Qobs , '-append', 'delimiter', '\t'); 

%% Prepare steramflow obs
d = readtable('G:\11_CRHM_cuchi\data\Q\cda_lev-q-t.csv'); % Read the table from the file
d = table2timetable(d)
d = retime(d, 'daily', 'mean')
Q = d.Discharge_cms_* 86400;
Q(isnan(Q))= 0;

    % create the obs matrix and save it in a mat file
t = datevec(d.Datetime);
t = t(:, 1:5);
t(:,4)=1;

Qobs= [t Q]; % compiled time and data together, and save it in a matlab format


fn = strcat('G:\11_CRHM_cuchi\data\processed\CuchiQ_2009_2019.mat');
save (fn, 'Qobs');  
 % create the obs text file
headerlines = {'Obs file, WRF PGW';
              'b 1 (C)';
              '$$ Missing ' ;
              '#####      b.1 '}
fp = strcat('G:\11_CRHM_cuchi\data\processed\CDAQ_2009_2019.obs');
fid = fopen(fp, 'wt');
for l = 1:numel(headerlines)
   fprintf(fid, '%s\n', headerlines{l});
end
fclose(fid);
dlmwrite(fp, Qobs , '-append', 'delimiter', '\t'); 

