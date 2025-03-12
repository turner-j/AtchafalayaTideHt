function [time,wl] = tp_regdata(time_raw,wl_raw,interval)
% tp_regdata.m
% Interpolates tide gauge time series to specified regular intervals 
% 
% NOTES:
%   - Duplicates are removed
%   - Gaps in the data greater than 3 hours are filled with NaNs
%   - The interpolated time series begins at the start of the first hour
%   - For further information refer to the Tide Peaks Toolbox User Manual
%
% Syntax:  [time,wl] = tp_regdata(time_raw,wl_raw,interval)
%
% Inputs:
%    time_raw - Input time vector in datetime format
%    wl_raw   - Input water levels (in metres)
%    interval - Specified output frequency in fraction of hours  
%
% Outputs:
%    time     - Regular interval time vector in datetime format
%    wl       - Interpolated water levels in metres
%
% BEFORE RUNNING THIS FUNCTION:
%    1. Create a new 'time_raw' variable MATLAB datetime "dd/MM/yyyy HH:mm" format
%       e.g. time_raw = datetime(time_raw,'ConvertFrom',datenum,'TimeZone','local','Format','dd/MM/yyyy HH:mm')
%    2. Create a new 'wl_raw' variable in metres with the required vertical reference
%    3. Replace any error flag values in "wl_raw" with NaNs
%    4. Divide the time series into annual lots January to December
%       - Where available, include 6-12 hours before 1-Jan and after 31-Dec
%         (ensures peaks occurring near New Years Eve are counted)
%
% Author: Karen Palmer, University of Tasmania
% Created: 1 October 2021 | Last revision: 9 February 2023

%------------- BEGIN CODE --------------
% Remove duplicates
idx = find(diff(time_raw)==0);
time_raw(idx) = [];
wl_raw(idx) = [];
% Insert a single Nan into wl_raw vector 1 hour into gaps > 3 hours
% When interpolation is executed wl values in the gaps will be NaNs
idx = find(diff(time_raw)>hours(3));                % index gaps longer than 3 hours
time2 = time_raw;                                   % make a copy of time_raw
wl2 = wl_raw;                                       % make a copy of wl_raw
timefill = NaT(size(idx));                          % pre-allocate time filler vector
% Insert extra time 1 interval into each gap
for i=1:numel(idx)
    tfill = time2(idx(i))+ interval/24;             % create extra time value 1 hour into gap
    timefill(i) = tfill;
end
time2 = [time2;timefill];
[timesort,torder] = sort(time2);
wlfill = NaN(size(idx));                            % create NaNs one per gap
wl2 = [wl2;wlfill];                                 % concatenate WLs with gap fillers
wlsort = wl2(torder);                               % sort in order of time series

t1 = time_raw(1);                                   % identify first time in data
t1 = dateshift(t1,'start','hour');                  % shift start time to beginning of first hour
time = t1:hours(interval):time_raw(end);            % create regular time vector at specified interval
time = time';                                       % transpose time to column vector
wl = interp1(timesort,wlsort,time);                 % interpolate edited time series to hourly vector

%------------- END OF CODE --------------
