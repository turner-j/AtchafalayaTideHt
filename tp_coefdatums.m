function [hhwss,mhws,mhhw,mhw,msl,mtl,mlw,mllw,mlws,F] = tp_coefdatums(coef)
% tp_coefdatums.m
% Takes the 'coef' structure output from UTide harmonic analysis and
% computes constituent-based tidal datums from standard formulae
% 
% NOTES:
%   - Equations from ICSM Australian Tides Manual SP9 v.5.0, 2018 and
%     Woodworth & Pugh, 2014. 'Sea Level Science'. Cambridge Press.
%   - Expected constituents are computed from one calendar year Jan-Dec *
%   - Additional datums can be computed from user specified equations
%   - Also computes the Mean Tide Level (MTL) and Tidal Form Factor (F)
%   - For further information refer to the Tide Peaks Toolbox User Manual
%
% Syntax:  [hhwss,mhws,mhhw,mhw,msl,mtl,mlw,mllw,mlws,F] = tp_coefdatums(coef)
%
% Inputs:
%    coef  - Tidal harmonic coefficient results structure from UTide
%
% Outputs:
%    hhwss - High High Water Solstices Springs
%    mhws  - Mean High Water Springs
%    mhhw  - Mean Higher High Water
%    mhw   - Mean High Water
%    msl   - Mean Sea Level
%    mtl   - Mean Tide Level
%    mlw   - Mean Low Water
%    mllw  - Mean Lower Low Water
%    mlws  - Mean Low Water Springs
%    F     - Form Factor
%
% BEFORE EXECUTING THIS FUNCTION (see m_tpdemo.m for template):
%   1. Execute 'tp_regdata.m' to interpolate raw data to regular intervals,
%      remove duplicates, and fill gaps with NaNs
%   2. Execute UTide functions to compute harmonic coefficients and
%      model basic and predicted tide
%
% Author: Karen Palmer, University of Tasmania
% email address: karen.palmer@utas.edu.au  
% Created: 3 December 2021 | Last revision: 9 February 2023

%------------- BEGIN CODE --------------
% Define constituent variables
Z0 = coef.mean;
M2 = coef.A(ismember(coef.name,'M2'));                  % amplitude of M2
S2 = coef.A(ismember(coef.name,'S2'));                  % amplitude of S2
O1 = coef.A(ismember(coef.name,'O1'));                  % amplitude of O1
K1 = coef.A(ismember(coef.name,'K1'));                  % amplitude of K1

% Compute constituent tidal datums from standard formulae
hhwss = Z0 + M2 + S2 + 1.4*(K1 + O1);
mhws = Z0 + (M2 + S2);
mhw = Z0 + M2;
mlw = Z0 - M2;
mlws = Z0 - (M2 + S2);
msl = Z0;
mtl = (mhw-mlw)/2;
mhhw = Z0 + (M2 + K1 + O1);
mllw = Z0 - (M2 + K1 + O1);
F = (K1 + O1)/(M2 + S2);                                % calculate Form Factor

%------------- END OF CODE --------------
