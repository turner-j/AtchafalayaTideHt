function [basicHW,basicLW,predHW,predLW,obsHW,obsLW] = tp_peaks(time,wl_basic,wl_pred,wl)
% tp_peaks.m
% Identifies high water (HW) peaks and low water (LW) troughs in tide gauge data
%
% NOTES:
%   - Requires MATLAB Signal Processing Toolbox
%   - HW/LW outputs are tables with the following headers
%        HW: [timeHW | HW] and LW: [timeLW | LW]
%   - High frequency peaks/troughs (e.g., wind waves and seiches) are excluded
%     by the minimum separation between peaks of 4 hours (MinPeakDistance = 4/24)
%   - The minimum height difference between consecutive peak and trough is
%     is set to 10 cm (MinPeakProminence = 0.1)
%   - For further information refer to the Tide Peaks Toolbox User Manual
%
% Syntax:  [basicHW,basicLW,predHW,predLW,obsHW,obsLW] = tp_peaks(time,wl_basic,wl_pred,wl)
%
% Inputs:
%    time     - Regular time vector in datetime format
%    wl_basic - Basic water levels modelled from 8 most influential constituents
%    wl_pred  - Predicted water levels modelled from all constituents
%    wl       - Observed water levels in metres
%
% Outputs:
%    basicHW  - Table of 'basic' high water (HW) time (datetime) and heights (m)
%    basicLW  - Table of 'basic' low water (LW) time (datetime) and heights (m)
%    predHW   - Table of predicted high water (HW) time (datetime) and heights (m)
%    predLW   - Table of predicted low water (LW) time (datetime) and heights (m)
%    obsHW    - Table of observed high water (HW) time (datetime) and heights (m)
%    obsLW    - Table of observed low water (LW) time (datetime) and heights (m)
%
% BEFORE EXECUTING THIS FUNCTION (see m_tpdemo.m for template):
%   1. Execute 'tp_regdata.m' to interpolate raw data to regular intervals,
%      remove duplicates, and fill gaps with NaNs
%   2. Execute UTide functions to compute harmonic coefficients and
%      model basic and predicted tide
%
% Author: Karen Palmer, University of Tasmania
% email address: karen.palmer@utas.edu.au  
% Created: 1 October 2021 | Last revision: 9 February 2023

%------------- BEGIN CODE --------------

% Execute MATLAB findpeaks function to locate tidal maximas in basic tide model
[HW,timeHW] = findpeaks(wl_basic,time,...
    'MinPeakDistance',4/24,...              % MinPeakDistance = 4 hrs apart
    'MinPeakProminence',0.1);               % MinPeakProminence = 0.1 metres
% Locate tidal minimas by findpeaks function on inverted model (-mod)
[LW,timeLW] = findpeaks(-wl_basic,time,...
    'MinPeakDistance',4/24,...              % MinPeakDistance = 4 hrs apart
    'MinPeakProminence',0.1);               % MinPeakProminence = 0.1 metres
LW = -LW;
% Copy basic peaks to new variables
bHW = HW; btimeHW = timeHW;
bLW = LW; btimeLW = timeLW;
% Output to tables
basicHW = table(timeHW,HW);
basicLW = table(timeLW,LW);
% Pre-allocate tables for predictions and observations
predHW = basicHW; predLW = basicLW;
obsHW = basicHW; obsLW = basicLW;
% Use findpeaks to locate tidal maximas in predicted tide
[HW,timeHW] = findpeaks(wl_pred,time,...
    'MinPeakDistance',0,...                 % MinPeakDistance = none
    'MinPeakProminence',0);                 % MinPeakProminence = none
% Locate tidal minimas by findpeaks command on inverted predictions (-wl_fit)
[LW,timeLW] = findpeaks(-wl_pred,time,...
    'MinPeakDistance',0,...                 % MinPeakDistance = none
    'MinPeakProminence',0);                 % MinPeakProminence = none
LW = -LW;
%
% Match basic model peaks to predicted maxima/minima
for i=1:numel(bHW)
    idx = find(timeHW > btimeHW(i)-4/24 & timeHW < btimeHW(i)+4/24);
    if ~isempty(idx)
    pk = HW(idx);
    hw = max(pk);
    t = timeHW(idx);
    t = t(pk == hw);
    % save predicted HW/timeHW variables to table
    predHW.HW(i) = hw;
    predHW.timeHW(i) = t(1);
    else
        predHW.HW(i) = NaN;
        predHW.timeHW(i) = NaT;
    end
end
for i=1:numel(bLW)
    idx = find(timeLW > btimeLW(i)-4/24 & timeLW < btimeLW(i)+4/24);
    if ~isempty(idx)
    pk = LW(idx);
    lw = min(pk);
    t = timeLW(idx);
    t = t(pk == lw);
    % overwrite HW/timeHW variables with observed values
    predLW.LW(i) = lw;
    predLW.timeLW(i) = t(1);
    else
        predLW.LW(i) = NaN;
        predLW.timeLW(i) = NaT;
    end
end
%
% Use findpeaks command to locate tidal maximas in observations
[HW,timeHW] = findpeaks(wl,time,...
    'MinPeakDistance',0,...                 % MinPeakDistance = none
    'MinPeakProminence',0);                 % MinPeakProminence = none
% Locate tidal minimas by findpeaks command on inverted observations (-mod)
[LW,timeLW] = findpeaks(-wl,time,...
    'MinPeakDistance',0,...                 % MinPeakDistance = none
    'MinPeakProminence',0);                 % MinPeakProminence = none
LW = -LW;
% Match basic model peaks to observed maxima/minima
for i=1:numel(bHW)
    idx = find(timeHW > btimeHW(i)-4/24 & timeHW < btimeHW(i)+4/24);
    if ~isempty(idx)
        pk = HW(idx);
        hw = max(pk);
        t = timeHW(idx);
        t = t(pk == hw);
        % overwrite HW/timeHW variables with observed values
        obsHW.HW(i) = hw;
        obsHW.timeHW(i) = t(1);
    else
        obsHW.HW(i) = NaN;
        obsHW.timeHW(i) = NaT;
    end
end
for i=1:numel(bLW)
    idx = find(timeLW > btimeLW(i)-4/24 & timeLW < btimeLW(i)+4/24);
    if ~isempty(idx)
        pk = LW(idx);
        lw = min(pk);
        t = timeLW(idx);
        t = t(pk == lw);
        % overwrite HW/timeHW variables with observed values
        obsLW.LW(i) = lw;
        obsLW.timeLW(i) = t(1);
    else
        obsLW.LW(i) = NaN;
        obsLW.timeLW(i) = NaT;
    end
end
%
%------------- END OF CODE --------------
