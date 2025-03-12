%% HARMONIC ANALYSIS with UTide
close all
clear
clc

load USAtfwind.mat
load CRMS0479.mat

USAtf_wind_data.dt  = USAtf_wind_data.date + timeofday(USAtf_wind_data.time);
USAtf_wind_data.dt.Format = 'dd.MM.yyyy HH:mm';
USAtf_wind_data = removevars(USAtf_wind_data, ["date","time"]);
TT = table2timetable(USAtf_wind_data);
TT2 = retime(TT,'hourly');

TT2 = TT2(ismember(TT2.dt,time_UTC),:);
wl_AHD = wl_AHD.*0.3048; % convert ft to m
wl_AHD = wl_AHD(ismember(time_UTC,TT2.dt));
%% Rank wind speeds
Q = quantile(TT2.wind_speed,3);
windrank = categorical(height(TT2),1);

for i = 1:length(TT2.wind_speed)
    if TT2.wind_speed(i)<=Q(1)
        windrank(i) = "low";
    elseif (TT2.wind_speed(i)>Q(1)) && (TT2.wind_speed(i)<=Q(2))
        windrank(i) = "mid";
    elseif (TT2.wind_speed(i)>Q(2)) && (TT2.wind_speed(i)<=Q(3))
        windrank(i) = "mid high";
    else
        windrank(i) = "high";
    end
end

windrank = windrank';
%% Create scalar time variable in MATLAB datenum format
wl_AHD(wl_AHD>20) = NaN;
T = timetable2table(TT2);
T = [T table(wl_AHD) table(windrank)];
T = removevars(T, "wind_speed");
T = rmmissing(T);

sz = size(T.wind_dir);
T.cardinal = deg2cardinal(T.wind_dir)';

ind = find(ismember(T.windrank=='high','mid high' & T.wind_dir=='S'));
ix = find(~ismember(T.windrank=='high','mid high' & T.wind_dir=='S'));
H = T(ind,1:3);
O = T(ix,1:3);
Z = T;

% Pre-allocate cell array for saving output variables
k = 2;          % define number of stations in cell array = 1
tpt_Hobart2021 = cell(k,1);

for i = 1:3
    if i == 1
        T = H;% high speed, northerly winds
    elseif i == 2% no high speed or northerly winds
        T = O;
    else 
        T = Z;% all wind speeds and directions
    end

    % Write station variables to cell array
    tpt_Hobart2021{i}.name = name;
    tpt_Hobart2021{i}.lat = lat;
    tpt_Hobart2021{i}.lon = lon;

    % PREPARE TIME SERIES
    % Define input variables
    time_raw = T.(1);
    wl_raw = T.(3);
    interval = 1; % Default output sampling interval is 1 hour
    % Raw data contains duplicates, gaps, and changes in sampling frequency and timing
    % Execute TP_REGDATA function to regulate data and fill gaps with NaNs
    [time,wl] = tp_regdata(time_raw,wl_raw,interval);
    
    % Remove unwanted variables from Workspace
    % Define year (OR write a loop to build a table with multiple years)
    tpt_Hobart2021{i}.year = 2024;

    % Write regulated time series to cell array
    VarNames = {'time_UTC','wl_AHD'};
    tpt_Hobart2021{i}.regdata = table(time,wl,'VariableNames',VarNames);

    % Create scalar time variable in MATLAB datenum format
    time_scalar = convertTo(time,"datenum");

    coef = ut_solv(time_scalar,wl,[],lat,'auto');
    % Fit tide predictions with UTide
    [wl_pred,~] = ut_reconstr(time_scalar,coef);
    % Create BASIC TIDE variable:
    % Extract 8th largest constituents by amplitude (excl long period)
    cnstit = coef.name(coef.aux.frq > 12/(24*365));       % choose constituents with frequency > 0.0014 cycles per hour (12 per year)
    cnstit = cnstit(1:8);
    % Reconstruct tidal predictions from basic constituents
    coef8 = ut_solv(time_scalar,wl,[],lat,cnstit);
    [wl_basic,~] = ut_reconstr(time_scalar,coef8);
    % Define Z0: mean sea level
    Z0 = coef.mean;
    M2 = coef.A(coef.name == "M2");                  % define M2 amplitude
    S2 = coef.A(coef.name == "S2");                  % define S2 amplitude
    K1 = coef.A(coef.name == "K1");                  % define K1 amplitude
    O1 = coef.A(coef.name == "O1");                  % define O1 amplitude
    F = (K1 + O1)/(M2 + S2);                         % calculate Form Factor (Defant, 1958)

    % Save variables to cell array
    tpt_Hobart2021{i}.regdata.wl_pred = wl_pred;
    tpt_Hobart2021{i}.coef = coef;

end

%% Store tidal variables
A = [(tpt_Hobart2021{1, 1}.coef.A)';(tpt_Hobart2021{1,1}.coef.A_ci)'];
T = array2table(A,"VariableNames",tpt_Hobart2021{1,1}.coef.name,'RowNames',{'Coefs','CI'});

%% Save the results as a CSV
writetable(T, 'htidecoefs.csv');
