function [MAX,MHWS,MHHW,MHW,MSL,MTL,MLW,MLLW,MLWS,MIN,Gt,Mn] = tp_pkdatums(time,wl,HWS,HHW,HW,LW,LLW,LWS)
% tp_pkdatums.m
% Takes a series of identified HHW/HW/LW/LLW peaks/troughs and averages
% 
% NOTES:
%   - Expected time series has a period of one calendar year Jan-Dec *
%   - Extracts the maximum (MAX) and minimum (MIN) levels in the period
%   - Computes the mean and mean absolute deviation (MAD) of MHWS/MHHW/MHW/MSL/MTL/MLW/MLLW/MLWS
%   - Computes the Greater (Gt) and Mean (Mn) Tidal Range
%   - Output datum variables are tables for each datum and respective deviation
%   - For further information refer to the Tide Peaks Toolbox User Manual
%
% Syntax:  [MAX,MHWS,MHHW,MHW,MSL,MTL,MLW,MLLW,MLWS,MIN,Gt,Mn] = tp_pkdatums(time,wl,HWS,HHW,HW,LW,LLW,LWS)
%
% Inputs:
%    time     - Time vector in datetime format, local standard time
%    wl       - Water levels in metres
%    interval - Time series interval in hours (e.g., 15 mins = 0.25)
%    HHW      - Time and height of diurnal higher high water peaks
%    HW       - Time and height of all high water peaks
%    LW       - Time and height of all low water troughs
%    LLW      - Time and height of diurnal lower low water peaks
%
% Outputs:
%    Year     - The datum year
%    MAX      - Maximum water level in period
%    MHWS     - Mean High Water Springs
%    MHHW     - Mean Higher High Water
%    MHW      - Mean High Water
%    MSL      - Mean Sea Level
%    MTL      - Mean Tide Level
%    MLW      - Mean Low Water
%    MLLW     - Mean Lower Low Water
%    MLWS     - Mean Low Water Springs
%    MIN      - Minimum water level in period
%    Gt       - Greater Tidal Range MHHW - MLLW
%    Mn       - Mean Tidal Range MHW - MLW
%
%
% Author: Karen Palmer, University of Tasmania
% email address: karen.palmer@utas.edu.au  
% Created: 4 October 2021 | Last revision: 29 August 2022

%------------- BEGIN CODE --------------
%
% Find the MAX and MIN water levels and time they occur
MAX = max(wl);
MAXtime = time(wl == MAX);
MAXtime = MAXtime(1);   % in case of two equal maximum values take the first
MAX = table(MAXtime,MAX);
MIN = min(wl);
MINtime = time(wl == MIN);
MINtime = MINtime(1);   % in case of two equal maximum values take the first
MIN = table(MINtime,MIN);
% Calculate Mean Sea Level
MSL = mean(wl,'omitnan');
MSLmad = mad(wl);
MSL = table(MSL,MSLmad);
% Calculate the means of classified peaks/troughs
MHWS = mean(HWS.HWS,'omitnan');
MHWSmad = mad(HWS.HWS);
MHWS = table(MHWS,MHWSmad);
MHHW = mean(HHW.HHW,'omitnan');
MHHWmad = mad(HHW.HHW);
MHHW = table(MHHW,MHHWmad);
MHW = mean(HW.HW,'omitnan');
MHWmad = mad(HW.HW);
MHW = table(MHW,MHWmad);
MLW = mean(LW.LW,'omitnan');
MLWmad = mad(LW.LW);
MLW = table(MLW,MLWmad);
MLLW = mean(LLW.LLW,'omitnan');
MLLWmad = mad(LLW.LLW);
MLLW = table(MLLW,MLLWmad);
MLWS = mean(LWS.LWS,'omitnan');
MLWSmad = mad(LWS.LWS);
MLWS = table(MLWS,MLWSmad);
% Compute the Mean Tide Level
MTL = (MHW.MHW + MLW.MLW)/2;
% Compute the Greater Tidal Range = MHHW - MLLW
Gt = MHHW.MHHW - MLLW.MLLW;
% Compute the Mean Tidal Range = MHW - MLW
Mn = MHW.MHW - MLW.MLW;
%
%------------- END OF CODE --------------
