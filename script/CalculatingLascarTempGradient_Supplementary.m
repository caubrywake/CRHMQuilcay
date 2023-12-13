%% Temperature Gradient, Quilcayhuanca

% the data analysed here is downloaded from 
%Mateo, E. I., B. G. Mark, R. Å. Hellström, M. Baraer, J. M. McKenzie, T. Condom, A. C. Rapre, G. Gonzales, J. Q. Gómez, R. C. C. Encarnación (2022). High temporal resolution hydrometeorological data collected in the tropical Cordillera Blanca, Peru (2004-2020), HydroShare, https://doi.org/10.4211/hs.35a670e6c5824ff89b3b74fe45ca90e0

% This is used to determine the temperature gradient as shown in the
% supplementray 1. 
close all
clear all
figdir = 'D:\11_CRHM_cuchi\fig\obs\'

folderPath = 'D:\11_CRHM_cuchi\data\lascar\';  % find the met files
files = dir([folderPath '*.csv']);  % Get a list of all CSV files in the folder

numFiles = numel(files); % Get the number of CSV files
clear data
% 4 filea: Cuchillacocha, Caasa de Agua, Llaganuco Wx and LLanganuco Port
label={'3955','3955','4122','4122','4355','4355','4561'}
% craete a time array
t1 = datetime('01-Jul-2006');
t2 = datetime('01-Jul-2019');
tt = t1:calmonths(1):t2;
% Import each file and apply the precipitation undercatch correction. 
for i = 1:numFiles
 filepath = [folderPath files(i).name]; % Get the full file path
 d = readtable(filepath); % Read the table from the file
d = table2timetable(d(:, 1:2));

dd_m = retime(d, tt, 'mean') ;
plot(dd_m.Datetime, dd_m.Temperature__C_); hold on
 data(:, i) = dd_m.Temperature__C_; % Store the modified table in the cell array
end

% use only (1,3,5,8)
plot(tt, data(:,[1,3,5,8]))
legend ('4355 m','3955 m','4122 m','4561 m')
ylabel({'Monthly Air Temperature ({\circ}C)'});
figname ='Lascar_monthlyTemp';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

elev = [4355,3955,4122,4561]


for i = 1:length(data)
coefficients = polyfit(elev*1000, data(i,[1,3,5,8]), 1);
gradient(i) = coefficients(1);
if isnan(gradient(i))
   continue
else 
scatter(elev, data(i,[1,3,5,8]), 'k.');  hold on
end 
end
lsline
ylabel({'Monthly Air Temperature ({\circ}C)'});
xlabel('Elevations (m.a.s.l.)')
figname ='Lascar_TempGradient';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

figure 
plot(tt, gradient)
months = month(tt)
for i = 1:12
a  =find(months==i);
b =gradient(a);
b(isnan(b))=[]
rates{i}=b;
meangrad(i) = nanmean(gradient(a))
meanstd (i)= nanstd(gradient(a))
end

x = [1:12]
y = meangrad;
errlow = meanstd;
errhigh=errlow;
figure

er = errorbar(x,y,errlow);    
er.Color = [0 0 0 ];                            
ylabel({'Monthly Air Temperature Gradient({\circ}C m{^-1})'});
xlabel('Months')
xlim([0.5 12.5])
figname ='Lascar_TempGradient_monthly';
saveas (gcf, strcat( figdir, figname, '.pdf'))
saveas (gcf, strcat(figdir, figname, '.png'))
savefig(gcf, strcat(figdir, figname))

mt ={'January'; 'February'; 'March'; 'April'; 'May'; 'June'; 'July'; 'August'; 'September'; 'October'; 'November'; 'December'}


savedir = 'D:\11_CRHM_cuchi\data\processed\'
t = table(mt, meangrad')
t.Properties.VariableNames= {'Month', 'TempGradident_C_m'};
writetable(t, strcat(savedir, 'Lascar_MonthlyTempGradident.csv'));
