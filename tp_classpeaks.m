function [HHW,LLW,HWS,LWS] = tp_classpeaks(basicHW,basicLW,HW,LW)
% tp_classpeaks
% Identifies high water peaks and low water troughs in tide data
% 
% NOTES:
%   - Classifies higher high water (HHW) and lower low water (LLW) from identified HW/LW
%   - Requires MATLAB Signal Processing Toolbox and UTide Toolbox (Codiga, 2011)
%   - HHW/HW/LW/LLW outputs are tables with "Time" and "HHW/HW/LW/LLW" variables respectively
%   - Values are given as Time (datetime) and Height (m)
%   - Diurnal peaks are identified from the maximum peak and minimum trough in
%     each sequence of two semidiurnal cycles, starting from the first
%     crossing of Z0. Predicted HW/LW are identified from the 'Basic' tide
%     fit using only the 8 most influential constituents (wl_basic)
%     (this reduces misclassifications in complex tidal regimes)
%   - For further information refer to the Tide Peaks Toolbox User Manual
%
% Syntax:  [HHW,LLW,HWS,LWS] = tp_classpeaks(basicHW,basicLW,HW,LW)
%
% Inputs as produced by the tp_peaks.m function:
%    basicHW  - High Water times and levels from the basic tide
%    basicLW  - Low Water times and levels from the basic tide
%    HW       - High Water peaks from the selected time series
%    LW       - Low Water troughs from the selected time series
%
% Outputs:
%    HHW      - Time and height of diurnal higher high water peaks
%    LLW      - Time and height of diurnal lower low water troughs
%    HWS      - Time and height of spring high water peaks (fortnightly)
%    LWS      - Time and height of spring low water troughs (fortnightly)
%
% BEFORE EXECUTING THIS FUNCTION (see m_tpdemo.m for template):
%   1. Execute 'tp_regdata.m' to interpolate raw data to regular intervals,
%      remove duplicates, and fill gaps with NaNs
%   2. Execute UTide functions to compute harmonic coefficients and
%      model basic and predicted tide
%   3. Execute 'tp_peaks.m' to identify HW/LW peaks and troughs
%
% Author: Karen Palmer, University of Tasmania
% email address: karen.palmer@utas.edu.au  
% Created: 1 October 2021 | Last revision: 28 April 2023

%------------- BEGIN CODE --------------

% Higher High Water | HHW
% Pre-allocate HHW & LLW vectors
fmt = 'dd/MM/yyyy hh:mm';
timeHHW = NaT(size(basicHW.timeHW),'Format',fmt);
HHW = NaN(size(basicHW.HW));
timeLLW = NaT(size(basicLW.timeLW),'Format',fmt);
LLW = NaN(size(basicLW.LW));
%
% Classify HHW peaks in basic tide
for a=1:2:numel(basicHW.HW)-2
    hw = [basicHW.HW(a);basicHW.HW(a+1);basicHW.HW(a+2)];                 % compare two consecutive values
    t = [basicHW.timeHW(a);basicHW.timeHW(a+1);basicHW.timeHW(a+2)];      % index respective times
    if t(2)-t(1)<16/24                                              % for indexed times < 18 hours apart
        hhw = max(hw(1:2));                                         % find the highest HW
        t = t(hw==hhw);                                             % find when HHW occurs
        HHW(a) = hhw;                                               % save HHW to matrix
        timeHHW(a) = t(1);                                          % if two peaks same height save first one
    elseif t(2)-t(1)>=16/24 && t(2)-t(1)<32/24 && t(3)-t(2)>=16/24 && t(3)-t(2)<32/24
        hhw = hw(1:2);                                              % save 1st two values
        t = t(1:2);                                                 % save time
        HHW(a:a+1) = hhw;                                           % save HHW variable to matrix
        timeHHW(a:a+1) = t;                                         % save time of HHW to matrix
    elseif t(2)-t(1)>=16/24 && t(2)-t(1)<32/24 && t(3)-t(2)<16/24
        lo = min(hw(1:3));                                          % compare values for lowest HW peak
        if lo == hw(3)                                              % if the 3rd peak is lowest
            hhw = hw(1:2);                                          % save 1st two values
            t = t(1:2);                                             % save 1st two times
            HHW(a:a+1) = hhw;                                       % write to matrix
            timeHHW(a:a+1) = t;                                     % write to matrix
        else
            HHW(a) = hw(1);                                         % save only the 1st HW peak
            timeHHW(a) = t(1);                                      % save 1st HW peak time
        end
    else
        HHW(a) = NaN;                                               % if no values found save NaN
        timeHHW(a) = NaT;                                           % if no values found save NaT
    end
end
% Remove leftover NaT/NaNs
timeHHW(isnat(timeHHW)) = [];
HHW(isnan(HHW)) = [];
% Copy basic peaks to new variables
bHHW = HHW; btimeHHW = timeHHW;
%
% Locate each basic HHW peak in HW peaks
for i=1:numel(HHW)
    idx = find(HW.timeHW > btimeHHW(i)-8/24 & HW.timeHW < btimeHHW(i)+8/24);
    if ~isempty(idx)
        pk = HW.HW(idx);
        hhw = max(pk);
        t = HW.timeHW(idx);
        t = t(pk == hhw);
        % overwrite HHW variables with new values
        HHW(i) = hhw;
        timeHHW(i) = t(1);
    else
        HHW(i) = NaN;
        timeHHW(i) = NaT;
    end
end
% Remove leftover NaT/NaNs
timeHHW(isnat(timeHHW)) = [];
HHW(isnan(HHW)) = [];

% High Water Springs | HWS
% Execute MATLAB findpeaks function to locate HWS in basic tide
if numel(bHHW)>3
    [hws,time_hws] = findpeaks(bHHW,btimeHHW,...
        'MinPeakDistance',8);               % MinPeakDistance = 8 days apart
    % Pre-allocate HWS array
    HWS = NaN(size(hws));
    timeHWS = NaT(size(time_hws),'format',fmt);
    %
    % Locate each basic HWS peak in HHW peaks
    for i=1:numel(HWS)
        idx = find(timeHHW > time_hws(i)-2 & timeHHW < time_hws(i)+2);
        if ~isempty(idx)
            pk = HHW(idx);
            ihws = max(pk);
            t = timeHHW(idx);
            t = t(pk == ihws);
            % overwrite HWS variables with new values
            HWS(i) = ihws;
            timeHWS(i) = t(1);
        else
            HWS(i) = NaN;
            timeHWS(i) = NaT;
        end
    end
    % Remove leftover NaT/NaNs
    timeHWS(isnat(timeHWS)) = [];
    HWS(isnan(HWS)) = [];
else
    timeHWS = NaT;
    HWS = NaN;
end

% Classify LLW peaks in basic tide
for a=1:2:numel(basicLW.LW)-2
    lw = [basicLW.LW(a);basicLW.LW(a+1);basicLW.LW(a+2)];                 % compare two consecutive values
    t = [basicLW.timeLW(a);basicLW.timeLW(a+1);basicLW.timeLW(a+2)];      % index respective times
    if t(2)-t(1)<16/24                                              % for values < 18 hours apart
        llw = min(lw(1:2));                                         % find the lowest LW
        t = t(lw==llw);                                             % find when LLW occurs
        LLW(a) = llw;
        timeLLW(a) = t(1);                                          % if two peaks of same height found choose first one
    elseif t(2)-t(1)>=16/24 && t(2)-t(1)<32/24 && t(3)-t(2)>=16/24 && t(3)-t(2)<32/24
        llw = lw(1:2);                                              % save value
        t = t(1:2);                                                 % save time
        LLW(a:a+1) = llw;                                           % save LLW variable to matrix
        timeLLW(a:a+1) = t;                                         % save time of LLW to matrix
    elseif t(2)-t(1)>=16/24 && t(2)-t(1)<32/24 && t(3)-t(2)<16/24
        hi = max(lw(1:3));                                          % compare values for lowest HW peak
        if hi == lw(3)                                              % if the 3rd peak is highest
            llw = lw(1:2);                                          % save 1st two values
            t = t(1:2);                                             % save time
            LLW(a:a+1) = llw;                                       % save 1st two values
            timeLLW(a:a+1) = t;                                     % save 1st two times
        else
            LLW(a) = lw(1);                                         % save HHW to matrix
            timeLLW(a) = t(1);                                      % save time to matrix
        end
    else
        LLW(a) = NaN;                                               % if no values found save NaN
        timeLLW(a) = NaT;                                           % if no values found save NaT
    end
end
%
% Remove leftover NaT/NaNs
timeLLW(isnat(timeLLW)) = [];
LLW(isnan(LLW)) = [];
% Copy basic LLW troughs to new variables
bLLW = LLW; btimeLLW = timeLLW;
%
% Locate each basic LLW trough in LW trough
for i=1:numel(LLW)
    idx = find(LW.timeLW > btimeLLW(i)-8/24 & LW.timeLW < btimeLLW(i)+8/24);
    if ~isempty(idx)
        pk = LW.LW(idx);
        llw = min(pk);
        t = LW.timeLW(idx);
        t = t(pk == llw);
        % overwrite LLW variables with observed values
        LLW(i) = llw;
        timeLLW(i) = t(1);
    else
        LLW(i) = NaN;
        timeLLW(i) = NaT;
    end
end
% Remove leftover NaT/NaNs
timeLLW(isnat(timeLLW)) = [];
LLW(isnan(LLW)) = [];

% Low Water Springs | LWS
% Execute MATLAB findpeaks function to locate LWS in basic tide
if numel(bLLW)>3
    [lws,time_lws] = findpeaks(-bLLW,btimeLLW,...
        'MinPeakDistance',8);               % MinPeakDistance = 8 days apart
    lws = -lws;
    % Pre-allocate LWS array
    LWS = NaN(size(lws));
    timeLWS = NaT(size(time_lws),'format',fmt);
    % Locate each basic LWS peak in LLW peaks
    for i=1:numel(LWS)
        idx = find(timeLLW > time_lws(i)-2 & timeLLW < time_lws(i)+2);
        if ~isempty(idx)
            pk = LLW(idx);
            ilws = min(pk);
            t = timeLLW(idx);
            t = t(pk == ilws);
            % overwrite LWS variables with new values
            LWS(i) = ilws;
            timeLWS(i) = t(1);
        else
            LWS(i) = NaN;
            timeLWS(i) = NaT;
        end
    end
    % Remove leftover NaT/NaNs
    timeLWS(isnat(timeLWS)) = [];
    LWS(isnan(LWS)) = [];
else
    timeLWS = NaT;
    LWS = NaN;
end
%
% Write to output tables
HWS = table(timeHWS,HWS);
HHW = table(timeHHW,HHW);
LLW = table(timeLLW,LLW);
LWS = table(timeLWS,LWS);
end
%
%------------- END OF CODE --------------
