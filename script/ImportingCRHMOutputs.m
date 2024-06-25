%% Import CRHM outputs and save as matlab format for the main simulations
% edited by Caroline Aubry-Wake, 2024-05-30
% contact at caroline.aubrywake@gmail.com

% This script imports the CRHM outputs (either .txt or .csv) and organises
% them and converts to a matlab format, which is used in further analysis.

%% Set-up
close all
clear all
cd 'G:\11_CRHM_cuchi\CRHM\output\v9\' % set to working folder
addpath 'G:\11_CRHM_cuchi\script\' % local directory with CRHM model results

%% Regular simulation
fileList = dir('R*.txt'); % get all the CRHM outputs starting with PeytoCUR_1OBS

for i = 1:numel(fileList)
fn = fileList(i).name; % file to import
H = importdata(fn,' ',2); %  % import headers 
D = importdata(fn) ; % import data
headers = regexp(H(1, 1), '\s+', 'split'); % split headers
headers = string(vertcat(headers{:})); % split headers
idxvar = [1, find(contains(headers,'(1)'))]; % select the number of variable in the file
% time is always the first one, followed by 2 or 3 variables
numvar = numel(idxvar); % number of variables

for ii= 1:numvar % for each variable, select all the column with data corresponding to that name
    varname =char(headers(idxvar(ii)));
    varname = strcat(varname(1:end-3)); % remove hru number from name
    if varname == 't' % excpetion is time - it does not have a number 
        varname = 'time';
    end 
Index = strfind(headers, varname);
Index = find(not(cellfun('isempty', Index)));
varname = strcat(varname);
assignin('base',varname,D.data(:, Index));
end 

end 
%convert the time to a datetime format
time= datetime(datevec(time+ 693960));
time = dateshift(time,'start','hour', 'nearest');
clear D H fileList fn headers i idxvar ii Index numvar pattern splitcells str varname  
save ('Cuchi_20230823.mat')

%% No ice simulation
cd 'G:\11_CRHM_cuchi\CRHM\output\v9\' % set to working folder

close all
clear all
fileList = dir('NoI*.txt'); % get all the CRHM outputs starting with PeytoCUR_1OBS

for i = 1:numel(fileList)
fn = fileList(i).name; % file to import
H = importdata(fn,' ',2); %  % import headers 
D = importdata(fn) ; % import data
headers = regexp(H(1, 1), '\s+', 'split'); % split headers
headers = string(vertcat(headers{:})); % split headers
idxvar = [1, find(contains(headers,'(1)'))]; % select the number of variable in the file
% time is always the first one, followed by 2 or 3 variables
numvar = numel(idxvar); % number of variables

for ii= 1:numvar % for each variable, select all the column with data corresponding to that name
    varname =char(headers(idxvar(ii)));
    varname = strcat(varname(1:end-3)); % remove hru number from name
    if varname == 't' % excpetion is time - it does not have a number 
        varname = 'time';
    end 
Index = strfind(headers, varname);
Index = find(not(cellfun('isempty', Index)));
varname = strcat(varname);
assignin('base',varname,D.data(:, Index));
end 

end 
basinflow_noice = basinflow;
basingw_noice =basingw;
%convert the time to a datetime format
time= datetime(datevec(time+ 693960));
time = dateshift(time,'start','hour', 'nearest');
clear D H fileList fn headers i idxvar ii Index numvar pattern splitcells str varname   basinflow basingw
save ('Cuchi_NoIce_20230823.mat')

%% No leakage simulation
cd 'G:\11_CRHM_cuchi\CRHM\output\v9\' % set to working folder

close all
clear all
fileList = dir('NoL*.txt'); % get all the CRHM outputs starting with PeytoCUR_1OBS

for i = 1:numel(fileList)
fn = fileList(i).name; % file to import
H = importdata(fn,' ',2); %  % import headers 
D = importdata(fn) ; % import data
headers = regexp(H(1, 1), '\s+', 'split'); % split headers
headers = string(vertcat(headers{:})); % split headers
idxvar = [1, find(contains(headers,'(1)'))]; % select the number of variable in the file
% time is always the first one, followed by 2 or 3 variables
numvar = numel(idxvar); % number of variables

for ii= 1:numvar % for each variable, select all the column with data corresponding to that name
    varname =char(headers(idxvar(ii)));
    varname = strcat(varname(1:end-3)); % remove hru number from name
    if varname == 't' % excpetion is time - it does not have a number 
        varname = 'time';
    end 
Index = strfind(headers, varname);
Index = find(not(cellfun('isempty', Index)));
varname = strcat(varname);
assignin('base',varname,D.data(:, Index));
end 

end 
basinflow_noleak = basinflow;
basingw_noleak =basingw;
%convert the time to a datetime format
time= datetime(datevec(time+ 693960));
time = dateshift(time,'start','hour', 'nearest');
clear D H fileList fn headers i idxvar ii Index numvar pattern splitcells str varname  basinflow basingw 
save ('Cuchi_NoLeak_20230823.mat')

%% No GW
cd 'G:\11_CRHM_cuchi\CRHM\output\v9\' % set to working folder

close all
clear all
fileList = dir('NoG*.txt'); % get all the CRHM outputs starting with PeytoCUR_1OBS

for i = 1:numel(fileList)
fn = fileList(i).name; % file to import
H = importdata(fn,' ',2); %  % import headers 
D = importdata(fn) ; % import data
headers = regexp(H(1, 1), '\s+', 'split'); % split headers
headers = string(vertcat(headers{:})); % split headers
idxvar = [1, find(contains(headers,'(1)'))]; % select the number of variable in the file
% time is always the first one, followed by 2 or 3 variables
numvar = numel(idxvar); % number of variables

for ii= 1:numvar % for each variable, select all the column with data corresponding to that name
    varname =char(headers(idxvar(ii)));
    varname = strcat(varname(1:end-3)); % remove hru number from name
    if varname == 't' % excpetion is time - it does not have a number 
        varname = 'time';
    end 
Index = strfind(headers, varname);
Index = find(not(cellfun('isempty', Index)));
varname = strcat(varname);
assignin('base',varname,D.data(:, Index));
end 

end 
basinflow_nogw = basinflow;
basingw_nogw =basingw;
%convert the time to a datetime format
time= datetime(datevec(time+ 693960));
time = dateshift(time,'start','hour', 'nearest');
clear D H fileList fn headers i idxvar ii Index numvar pattern splitcells str varname  basinflow basingw 
save ('Cuchi_NoGW_20230823.mat')

% clear all
% d = readtable('NoGW_RunFlow.csv');
% tt_nogw = d.Datetime;
% gw_nogw  = d.Basingw_m3hr_nogw;
% bf_nogw = d.Basinflow_m3hr_nogw;
% clear d
% save ('Cuchi_NoGW_20230823.mat')
