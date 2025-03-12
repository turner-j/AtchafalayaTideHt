function [surgeHW,surgeLW] = tp_skewsurge(predHW,predLW,obsHW,obsLW)
% tp_skewsurge - Computes heights and timing offsets between predicted and
% observed High Waters (HW) and Low Waters (LW
%
%   - surgeHW/surgeLW outputs are tables with the following headers
%        surgeHW: [predTIME | pred | obsTIME | obs | shift | surge]
%        surgeLW: [predTIME | pred | obsTIME | obs | shift | surge]
%   - For further information refer to the Tide Peaks Toolbox User Manual
%
% Syntax:  [surgeHW,surgeLW] = tp_skewsurge(predHW,predLW,obsHW,obsLW)
%
% Inputs:
%    predHW  - Table of times (datetime) and heights (m) of predicted HWs
%    predLW  - Table of times (datetime) and heights (m) of predicted HWs
%    obsHW   - Table of times (datetime) and heights (m) of predicted HWs
%    obsLW   - Table of times (datetime) and heights (m) of predicted HWs
%              Input tables have format: [timeHW | HW]
%
% Outputs:
%    surgeHW - Table of predicted/observed times and heights of HWs and the
%              differences between them
%    surgeLW - Table of predicted/observed times and heights of LWs and the
%              differences between them
%
% BEFORE EXECUTING THIS FUNCTION (see m_tpdemo.m for template):
%   1. Execute 'tp_regdata.m' to interpolate raw data to regular intervals,
%      remove duplicates, and fill gaps with NaNs
%   2. Execute UTide functions to compute harmonic coefficients and
%      model basic and predicted tide
%   3. Execute 'tp_peaks.m' function to identify HW and LW peaks and
%      troughs for both predicted and observed water level time series
%
% Author: Karen Palmer, University of Tasmania
% email address: karen.palmer@utas.edu.au  
% Created: 2 December 2021 | Last revision: 9 February 2023

%------------- BEGIN CODE --------------
%
% Define output table variable names
colNames = {'predTIME','pred','obsTIME','obs','shift','surge'};
% High water skew surge
surgeHW = table('Size',[height(obsHW) 6],...
    'VariableTypes',["datetime","double","datetime","double","duration","double"],...
    'VariableNames',colNames);
for i=1:height(predHW)
    timeobs = obsHW.timeHW(i); timepred = predHW.timeHW(i);
    obs = obsHW.HW(i); pred = predHW.HW(i);
    surgeHW.predTIME(i) = timepred;
    surgeHW.pred(i) = pred;
    surgeHW.obsTIME(i) = timeobs;
    surgeHW.obs(i) = obs;
    surgeHW.shift(i) = diff([timeobs timepred]);      % calculate phase shift in hours
    surgeHW.surge(i) = obs - pred;                    % calculate height difference between peaks
end
% Low water skew surge
surgeLW = table('Size',[height(obsLW) 6],...
    'VariableTypes',["datetime","double","datetime","double","duration","double"],...
    'VariableNames',colNames);
for i=1:height(predLW)
    timeobs = obsLW.timeLW(i); timepred = predLW.timeLW(i);
    obs = obsLW.LW(i); pred = predLW.LW(i);
    surgeLW.predTIME(i) = timepred;
    surgeLW.pred(i) = pred;
    surgeLW.obsTIME(i) = timeobs;
    surgeLW.obs(i) = obs;
    surgeLW.shift(i) = diff([timeobs timepred]);      % calculate phase shift in hours
    surgeLW.surge(i) = obs - pred;                    % calculate height difference between troughs
end
end
%
%------------- END OF CODE --------------
