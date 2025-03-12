%% TEMPLATE Tide Peaks Toolbox execution script
% Executes Tide Peaks Toolbox functions and demonstrates using sample data
% Sample data contained in directory "user\...\tidepeaks\HobartTG2021"
% For further information refer to the Tide Peaks Toolbox User Manual
%
% INPUTS:
%   - Minimum 1-hourly interval tide gauge time series [datetime] [water level (m)]
%   - Station name [string format]
%   - Station latitude in decimal degrees
%   - Station longitude in decimal degrees
%
% OUTPUTS:
%   - Figures illustrating each of the 8 process steps
%   - Cell array containing output variables
%
% REQUIREMENTS:
%   - UTide Unified Tidal Analysis and Prediction Functions' by Daniel Codiga (2017)
%     https://au.mathworks.com/matlabcentral/fileexchange/46523-utide-unified-tidal-analysis-and-prediction-functions
%     NOTE: Ensure the UTide folder is included in search path
%   - MATLAB Signal Processing Toolbox
%   - MATLAB Statistics Toolbox
% VERSION:
%   - The Tide Peaks Toolbox was created and tested in MATLAB R2020
%
% Author: Karen Palmer, University of Tasmania
% Created: 9 February 2023 | Last revision: 11 June 2023
%
% Attribution: If you use this toolbox please cite the following paper:
% Palmer, K., Hunter, J. R., Watson, C.S., Hague, B.H., & Power, H. E., 2023.
% An improved method for computing tidal datums. Coastal Engineering.
% DOI: https://doi.org/10.1016/j.coastaleng.2023.104354 
%
% Copyright 2023 Karen Palmer
% Open source initiative 3-clause BSD License applies
% https://opensource.org/license/bsd-3-clause/

%% Get started
close all
clear
clc

% Define colours for figures
blue = [31,86,163]/255;
darkred = [0.6353,0.0784,0.1843];
grey = [0.5 0.5 0.5];

%% 1. TIDE GAUGE DATA
% SAMPLE DATA from Hobart tide gauge provided by the Bureau of Meteorology
  % Time in Universal Coordinated Time (UTC)
  % Water levels in metres above Australian Height Datum (AHD-Tas83)
  % NOTE: time series contains gaps, duplicates, and irregular timing
% Load BOM Hobart tide gauge data file for 2021 from Toolbox directory
% load HobartTG2021.mat
load CRMS0479.mat
% Includes the following variables
% lat = 29.52; % Latitude
% lon = -91.45; % Longitude
% name = 'CRMS0479'; % Station name
% time_UTC = Datemmddyyyy + timeofday(Timehhmmss);
% time_UTC.Format = 'dd.MM.yyyy HH:mm';
%  Time UTC: time_UTC = [87529x1 datetime] in approx 10-minute intervals
% Water levels: wl_AHD = [87529x1 double] in metres above AHD-Tas83
% wl_AHD = wl_AHD./0.3048;
time_UTC(wl_AHD>20)=[];
wl_AHD(wl_AHD>20) = [];

% Pre-allocate cell array for saving output variables
k = 1;          % define number of stations in cell array = 1
tpt_Hobart2021 = cell(k,1);
% Write station variables to cell array
tpt_Hobart2021{k}.name = name;
tpt_Hobart2021{k}.lat = lat;
tpt_Hobart2021{k}.lon = lon;

% Preliminary plot: Observed water level time series 14 days
xstart = time_UTC(1)-7;
xend = time_UTC(end)+7;
WLmin = min(wl_AHD);
WLmax = max(wl_AHD);
r = range(wl_AHD);
ymin = WLmin;
ymax = WLmax;

fig1 = figure;
fig1.Units = 'centimeters';
fig1.Position = [2 5 20 12];
    hold on
    plot(time_UTC,wl_AHD,'color',blue,'linewidth',0.5);
    yline(0,'k');
    xlim([xstart xend]);
    ylim([ymin ymax]);
    ylabel('Metres above AHD-Tas83');
    xticks(time_UTC(1)+1:60:time_UTC(end));
    xtickformat('MMM-yyyy')
    % datetick('x','mmm-yyyy','keeplimits','keepticks');
    title('1. Tide gauge data');
    legend('Observed water level');
    box on
    hold off

% export_fig('Hobart_1-data.png','-r300','-dpng','-transparent','-painters')
%% 2. PREPARE TIME SERIES
% Define input variables
time_raw = time_UTC;
wl_raw = wl_AHD;
interval = 1;           % Default output sampling interval is 1 hour
% Raw data contains duplicates, gaps, and changes in sampling frequency and timing
% Execute TP_REGDATA function to regulate data and fill gaps with NaNs
[time,wl] = tp_regdata(time_raw,wl_raw,interval);
% Remove unwanted variables from Workspace
% clear time_raw wl_raw time_UTC wl_AHD
% Define year (OR write a loop to build a table with multiple years)
year = 2024;
tpt_Hobart2021{k}.year = year;

% Write regulated time series to cell array
VarNames = {'time_UTC','wl_AHD'};
tpt_Hobart2021{k}.regdata = table(time,wl,'VariableNames',VarNames);

% Plot 14-day sample with data points
xstart = datetime(year,1,1,0,0,0);
xend = xstart + 14;

fig2 = figure;
fig2.Units = 'centimeters';
fig2.Position = [3 5 20 12];
    hold on
    plot(time,wl,'color',blue,'linewidth',1.5);
    plot(time,wl,'ks','markerfacecolor','k','markersize',2.5);
    yline(0,'k');
    xlim([xstart xend]);
    ylim([ymin ymax]);
    xticks(xstart:3:xend);
    xtickformat('dd-MMM')
    % datetick('x','dd-mmm','keeplimits','keepticks');
    ylabel('Water level (mAHD)');
    xlabel(year);
    title('2. Prepare time series');
    legend('Observed water level','Sample point');
    box on
    hold off
    
% export_fig('Hobart_2-regdata.png','-r300','-dpng','-transparent','-painters')
    
%% 3. HARMONIC ANALYSIS with UTide
% Create scalar time variable in MATLAB datenum format
time_scalar = datenum(time);
% Solve harmonic coefficients with UTide
% Proceed if > 6mths data present (default)
if sum(~isnan(wl)) >  (1/interval)*24*180
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
    % Define tidal Form Factor (Pugh & Woodworth, 2014)
    Ftxt = ["Tidal form is semidiurnal";...
        "Tidal form is mixed, mainly semidiurnal";...
        "Tidal form is mixed, mainly diurnal";...
        "Tidal form is diurnal"];
        M2 = coef.A(coef.name == "M2");                  % define M2 amplitude
        S2 = coef.A(coef.name == "S2");                  % define S2 amplitude
        K1 = coef.A(coef.name == "K1");                  % define K1 amplitude
        O1 = coef.A(coef.name == "O1");                  % define O1 amplitude
        F = (K1 + O1)/(M2 + S2);                         % calculate Form Factor (Defant, 1958)
        % Save text describing semidiurnal/mixed/diurnal tidal regime
        if F>=0 && F<0.25
            Ftxt = sprintf('%s F=%3.2f',Ftxt(1),F);
        elseif F>=0.25 && F<1.5
            Ftxt = sprintf('%s F=%3.2f',Ftxt(2),F);
        elseif F>=1.5 && F<3.0
            Ftxt = sprintf('%s F=%3.2f',Ftxt(3),F);
        elseif F>=3.0
            Ftxt = sprintf('%s F=%3.2f',Ftxt(4),F);
        else
            Ftxt = "Tidal form is not defined";
        end
        % display tidal form factor text to command window
        sprintf('%s: %s',name,Ftxt)
else
    % if < 6mths data present save NaNs
    coef = {};
    wl_pred = NaN(size(wl));
    wl_basic = NaN(size(wl));
    Ftxt = "Tidal form is not defined";
end

% Save variables to cell array
tpt_Hobart2021{k}.regdata.wl_pred = wl_pred;
tpt_Hobart2021{k}.coef = coef;
tpt_Hobart2021{k}.Ftxt = Ftxt;

% Plot observed, predicted, and basic tide over 7 days
xstart = datetime(year,1,24,0,0,0);
xend = xstart+7;
fig3 = figure;
fig3.Units = 'centimeters';
fig3.Position = [4 5 20 12];
    hold on
    plot(time,wl,'color',blue,'linewidth',1.5);
    plot(time,wl_pred,'color',darkred,'linewidth',1.5);
    plot(time,wl_basic,'color',grey,'linewidth',1.5);
    yline(Z0,'color',darkred,'linestyle','--');
    yline(0,'k');
    text(xstart+0.05,Z0,'MSL','color','k','verticalalignment','bottom')
    % Add text for tidal form description
    text(xstart+0.2,ymax-r/8,Ftxt,'horizontalalignment','left');
    xlim([xstart xend]);
    ylim([ymin ymax]);
    xticks(xstart:1:xend);
    xtickformat('dd-MMM')
    % datetick('x','dd-mmm','keeplimits','keepticks');
    ylabel('Water level (mAHD)');
    xlabel(year);
    title('3. Harmonic analysis');
    legend('Observed','Predicted','Basic');
    box on
    hold off

% export_fig('Hobart_3-harmonic.png','-r300','-dpng','-transparent','-painters')
    
%% 4. CONSTITUENT TIDAL DATUMS
% Execute TP_COEFDATUMS function for constituent-based 2021 tidal datums
if ~isnan(coef.mean)
    [hhwss,mhws,mhhw,mhw,msl,mtl,mlw,mllw,mlws,F] = tp_coefdatums(coef);
    coefdatums = [year,hhwss,mhws,mhhw,mhw,msl,mtl,mlw,mllw,mlws,F];
else
    coefdatums(y,:) = [year,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end
% Define table output column names
VarNames = {'Year','hhwss','mhws','mhhw','mhw','msl','mtl','mlw','mllw','mlws','F'};
% Write constituent tidal datums to table       
coefdatums = array2table(coefdatums,'VariableNames',VarNames);
% Display constituent tidal datums to command window
disp(coefdatums)

% Save constituent-based tidal datums to cell array
tpt_Hobart2021{k}.coefdatums = coefdatums;

% Plot tidal datums computed from Z0/M2/S2/K1/O1 constituents
fig4 = figure;
fig4.Units = 'centimeters';
fig4.Position = [29 5 3 12];
    hold on
    line([1 1],[ymin ymax],'color','k','linewidth',0.1);
    line([2 2],[ymin ymax],'color','k','linewidth',0.1);
    line([0.5 3],[msl msl],'color','k',...
        'linestyle','--','linewidth',0.5);
    % Constituent tidal datums for {Year}
    text(0,coefdatums.mhws,'MHWS','fontsize',8);
    plot(1,coefdatums.mhws,'kd',...
        'markersize',8,'markerfacecolor','k','linewidth',1);
    text(1.15,coefdatums.mhws,sprintf('%3.2f',coefdatums.mhws),...
        'fontsize',9,'horizontalalignment','left');
    text(0,coefdatums.mhhw,'MHHW','fontsize',8);
    plot(1,coefdatums.mhhw,'kv',...
        'markersize',6,'markerfacecolor','k','linewidth',1);
    text(1.15,coefdatums.mhhw,sprintf('%3.2f',coefdatums.mhhw),...
        'fontsize',9,'horizontalalignment','left');
    text(0,coefdatums.mhw,'MHW','fontsize',8);
    plot(1,coefdatums.mhw,'kx',...
        'markersize',8,'linewidth',2);
    text(1.15,coefdatums.mhw,sprintf('%3.2f',coefdatums.mhw),...
        'fontsize',9,'horizontalalignment','left');
    text(0,coefdatums.msl,'MSL','fontsize',8);
    plot(1,coefdatums.msl,'ko',...
        'markersize',6,'markerfacecolor','k','linewidth',1);
    text(1.15,coefdatums.msl,sprintf('%3.2f',coefdatums.msl),...
        'fontsize',9,'horizontalalignment','left');
    text(0,coefdatums.mlw,'MLW','fontsize',8);
    plot(1,coefdatums.mlw,'k+',...
        'markersize',8,'linewidth',2);
    text(1.15,coefdatums.mlw,sprintf('%3.2f',coefdatums.mlw),...
        'fontsize',9,'horizontalalignment','left');
    text(0,coefdatums.mllw,'MLLW','fontsize',8);
    plot(1,coefdatums.mllw,'k^',...
        'markersize',6,'markerfacecolor','k','linewidth',1);
    text(1.15,coefdatums.mllw,sprintf('%3.2f',coefdatums.mllw),...
        'fontsize',9,'horizontalalignment','left');
    text(0,coefdatums.mlws,'MLWS','fontsize',8);
    plot(1,coefdatums.mlws,'ks',...
        'markersize',9,'markerfacecolor','k','linewidth',1);
    text(1.15,coefdatums.mlws,sprintf('%3.2f',coefdatums.mlws),...
        'fontsize',9,'horizontalalignment','left');

    xlim([-0.1 1.7]);
    xticklabels([]);
    xticks([]);
    xlabel(year);
    ylim([ymin ymax]);
    yticklabels([]);
    yticks([]);
    title(sprintf('4. Constituent\ntidal datums'));
    box on
    hold off
    
% export_fig('Hobart_4-constitdatum.png','-r300','-dpng','-transparent','-painters')

%% 5. PEAK/TROUGH IDENTIFICATION
% Execute TP_PEAKS function to identify HW/LW
% Outputs are tables with format: [time | height(m)]
[basicHW,basicLW,predHW,predLW,obsHW,obsLW] = tp_peaks(time,wl_basic,wl_pred,wl);
% Remove peaks identified outside of the present calendar year
ystart = datetime(year,1,1,0,0,0);
yend = datetime(year+1,1,1,0,0,0);
predHW = predHW(~(obsHW.timeHW<ystart | obsHW.timeHW>yend),:);
predLW = predLW(~(obsLW.timeLW<ystart | obsLW.timeLW>yend),:);
obsHW = obsHW(~(obsHW.timeHW<ystart | obsHW.timeHW>yend),:);
obsLW = obsLW(~(obsLW.timeLW<ystart | obsLW.timeLW>yend),:);
        
% Save predicted and observed HW/LW tables to cell array
tpt_Hobart2021{k}.predHW = predHW;
tpt_Hobart2021{k}.predLW = predLW;
tpt_Hobart2021{k}.obsHW = obsHW;
tpt_Hobart2021{k}.obsLW = obsLW;

% Plot identified HWs and LWs
xstart = datetime(year,1,1,0,0,0);
xend = datetime(year,2,1,0,0,0);
fig5 = figure;
fig5.Units = 'centimeters';
fig5.Position = [5 5 22 12];
    hold on
    o = plot(time,wl,'color',blue,'linewidth',0.5);
    p = plot(time,wl_pred,'k-','linewidth',1.5);
    yline(0,'k');
    h = plot(predHW.timeHW,predHW.HW,'kx',...
        'markersize',8,'linewidth',2);
    l = plot(predLW.timeLW,predLW.LW,'k+',...
        'markersize',8,'linewidth',2);
    xlim([xstart xend]);
    xlabel(year);
    ylim([ymin ymax]);
    ylabel('Water level (mAHD)');
    xticks(xstart:7:xend);
    xtickformat('dd-MMM')
    % datetick('x','dd-mmm','keeplimits','keepticks');
    title('5. Peak/trough identification');
    legend([o p h l],'Observed','Predicted','HW','LW',...
        'location','eastoutside');
    box on
    hold off
    
export_fig('Hobart_5-idpeaks.png','-r300','-dpng','-transparent','-painters')
%% 6. PEAK/TROUGH CLASSIFICATION
% Execute TP_CLASSPEAKS function to classify HWS/HHW/LLW/LWS
  % [HHW,LLW,HWS,LWS] = tp_classpeaks(basicHW,basicLW,predHW,predLW)
[predHHW,predLLW,predHWS,predLWS] = tp_classpeaks(basicHW,basicLW,predHW,predLW);
[obsHHW,obsLLW,obsHWS,obsLWS] = tp_classpeaks(basicHW,basicLW,obsHW,obsLW);

% Save predicted and observed HWS/HHW/LLW/LWS tables to cell array
tpt_Hobart2021{k}.predHWS = predHWS;
tpt_Hobart2021{k}.predHHW = predHHW;
tpt_Hobart2021{k}.predLLW = predLLW;
tpt_Hobart2021{k}.predLWS = predLWS;
tpt_Hobart2021{k}.obsHWS = obsHWS;
tpt_Hobart2021{k}.obsHHW = obsHHW;
tpt_Hobart2021{k}.obsLLW = obsLLW;
tpt_Hobart2021{k}.obsLWS = obsLWS;

% Plot classified HWS/HHW/HW/LW/LLW/LWS
% NOTE: Offset HWS/HHW/LLW/LWS markers for clarity
xstart = datetime(year,1,1,0,0,0);
xend = datetime(year,2,1,0,0,0);

fig6 = figure;
fig6.Units = 'centimeters';
fig6.Position = [6 5 22 12];
    hold on
    o = plot(time,wl,'color',blue,'linewidth',0.5);
    p = plot(time,wl_pred,'k-','linewidth',1.5);
    yline(0,'k');
    h = plot(predHW.timeHW,predHW.HW,'kx',...
        'markersize',8,'linewidth',2);
    hh = plot(predHHW.timeHHW,predHHW.HHW+0.1,'kv',...
        'markersize',6,'markerfacecolor','k','linewidth',1);
    hs = plot(predHWS.timeHWS,predHWS.HWS+0.2,'kd',...
        'markersize',8,'markerfacecolor','k','linewidth',1);
    l = plot(predLW.timeLW,predLW.LW,'k+',...
        'markersize',8,'linewidth',2);
    ll = plot(predLLW.timeLLW,predLLW.LLW-0.1,'k^',...
        'markersize',6,'markerfacecolor','k','linewidth',1);
    ls = plot(predLWS.timeLWS,predLWS.LWS-0.2,'ks',...
        'markersize',9,'markerfacecolor','k','linewidth',1);
    xlim([xstart xend]);
    xlabel(year);
    ylabel('Water level (mAHD)');
    xticks(xstart:7:xend);
    xtickformat('dd-MMM')
    % datetick('x','dd-mmm','keeplimits','keepticks');
    ylim([ymin ymax]);
    title('6. Peak/trough classification');
    legend([o p hs hh h l ll ls],'Observed','Predicted',...
        'HWS','HHW','HW','LW','LLW','LWS',...
        'location','eastoutside');
    box on
    hold off
    
export_fig('Hobart_6-classpeaks.png','-r300','-dpng','-transparent','-painters')
%% 7. PEAK-BASED TIDAL DATUMS
% Execute TP_PKDATUMS function to compute peak-based tidal datums: PREDICTED
[MAX,MHWS,MHHW,MHW,MSL,MTL,MLW,MLLW,MLWS,MIN,Gt,Mn] = tp_pkdatums(time,wl_pred,predHWS,predHHW,predHW,predLW,predLLW,predLWS);
% Populate table with yearly means for predicted tide
preddatums = [year,MAX.MAX,MHWS.MHWS,MHWS.MHWSmad,...
    MHHW.MHHW,MHHW.MHHWmad,MHW.MHW,MHW.MHWmad,MSL.MSL,MSL.MSLmad,MTL,...
    MLW.MLW,MLW.MLWmad,MLLW.MLLW,MLLW.MLLWmad,...
    MLWS.MLWS,MLWS.MLWSmad,MIN.MIN,Gt,Mn];       
% Execute TP_PKDATUMS function to compute peak-based tidal datums: OBSERVED
[MAX,MHWS,MHHW,MHW,MSL,MTL,MLW,MLLW,MLWS,MIN,Gt,Mn] = tp_pkdatums(time,wl,obsHWS,obsHHW,obsHW,obsLW,obsLLW,obsLWS);
% Populate table with yearly means for observed water levels
obsdatums = [year,MAX.MAX,MHWS.MHWS,MHWS.MHWSmad,...
    MHHW.MHHW,MHHW.MHHWmad,MHW.MHW,MHW.MHWmad,MSL.MSL,MSL.MSLmad,MTL,...
    MLW.MLW,MLW.MLWmad,MLLW.MLLW,MLLW.MLLWmad,...
    MLWS.MLWS,MLWS.MLWSmad,MIN.MIN,Gt,Mn];

% Define table output column names
VarNames = {'Year','MAX','MHWS','MHWSmad','MHHW','MHHWmad','MHW','MHWmad','MSL','MSLmad','MTL','MLW','MLWmad','MLLW','MLLWmad','MLWS','MLWSmad','MIN','Gt','Mn'};
% Write peak-based tidal datums to tables
preddatums = array2table(preddatums,'VariableNames',VarNames);
obsdatums = array2table(obsdatums,'VariableNames',VarNames);
% Display peak/trough tidal datums to command window
disp(preddatums)
disp(obsdatums)

% Save tables to cell array
tpt_Hobart2021{k}.preddatums = preddatums;
tpt_Hobart2021{k}.obsdatums = obsdatums;

% Plot computed datums from predicted peaks and troughs
fig7 = figure;
fig7.Units = 'centimeters';
fig7.Position = [32 5 3 12];
    hold on
    line([1 1],[ymin ymax],'color','k','linewidth',0.1);
    line([2 2],[ymin ymax],'color','k','linewidth',0.1);
    line([0.5 3],[msl msl],'color','k',...
        'linestyle','--','linewidth',0.5);
    % PREDICTED tidal datums for {Year}
    text(0,preddatums.MHWS,'MHWS','fontsize',8);
    plot(1,preddatums.MHWS,'kd',...
        'markersize',8,'markerfacecolor','k','linewidth',1);
    text(1.15,preddatums.MHWS,sprintf('%3.2f',preddatums.MHWS),...
        'fontsize',9,'horizontalalignment','left');
    text(0,preddatums.MHHW,'MHHW','fontsize',8);
    plot(1,preddatums.MHHW,'kv',...
        'markersize',6,'markerfacecolor','k','linewidth',1);
    text(1.15,preddatums.MHHW,sprintf('%3.2f',preddatums.MHHW),...
        'fontsize',9,'horizontalalignment','left');
    text(0,preddatums.MHW,'MHW','fontsize',8);
    plot(1,preddatums.MHW,'kx',...
        'markersize',8,'linewidth',2);
    text(1.15,preddatums.MHW,sprintf('%3.2f',preddatums.MHW),...
        'fontsize',9,'horizontalalignment','left');
    text(0,preddatums.MSL,'MSL','fontsize',8);
    plot(1,preddatums.MSL,'ko',...
        'markersize',6,'markerfacecolor','k','linewidth',1);
    text(1.15,preddatums.MSL,sprintf('%3.2f',preddatums.MSL),...
        'fontsize',9,'horizontalalignment','left');
    text(0,preddatums.MLW,'MLW','fontsize',8);
    plot(1,preddatums.MLW,'k+',...
        'markersize',8,'linewidth',2);
    text(1.15,preddatums.MLW,sprintf('%3.2f',preddatums.MLW),...
        'fontsize',9,'horizontalalignment','left');
    text(0,preddatums.MLLW,'MLLW','fontsize',8);
    plot(1,preddatums.MLLW,'k^',...
        'markersize',6,'markerfacecolor','k','linewidth',1);
    text(1.15,preddatums.MLLW,sprintf('%3.2f',preddatums.MLLW),...
        'fontsize',9,'horizontalalignment','left');
    text(0,preddatums.MLWS,'MLWS','fontsize',8);
    plot(1,preddatums.MLWS,'ks',...
        'markersize',9,'markerfacecolor','k','linewidth',1);
    text(1.15,preddatums.MLWS,sprintf('%3.2f',preddatums.MLWS),...
        'fontsize',9,'horizontalalignment','left');
    xlim([-0.1 1.7]);
    xticklabels([]);
    xlabel(year);
    xticks([]);
    ylim([ymin ymax]);
    yticklabels([]);
    yticks([]);
    title(sprintf('7. Peak/trough\ntidal datums'));
    box on
    hold off
    
export_fig('Hobart_7-peakdatums.png','-r300','-dpng','-transparent','-painters')
%% 8. SKEW SURGES
% Execute TP_SKEWSURGE function
% Outputs are tables with format:
% [Time(pred) | Height(m,pred) | Time(obs) | Height(m,obs) | Time(hms,diff) | Height(m,diff)]
% NOTE: time difference increments depend on interval (default 1 hour)
% NOTE: Surge (height difference) is +ve for higher and -ve for lower than predicted
[surgeHW,surgeLW] = tp_skewsurge(predHW,predLW,obsHW,obsLW);

% Save skew surge variables to cell array
tpt_Hobart2021{k}.surgeHW = surgeHW;
tpt_Hobart2021{k}.surgeLW = surgeLW;

% Plot skew surge heights at HWs and LWs
xstart = datetime(year,1,24,0,0,0);
xend = xstart +7;

fig8 = figure;
fig8.Units = 'centimeters';
fig8.Position = [7 5 22 12];
    hold on
    o = plot(time,wl,'color',blue,'linewidth',1.5);
    p = plot(time,wl_pred,'color','k','linewidth',1.5);
    yline(0,'k');
    h1 = plot(predHW.timeHW,predHW.HW,'kx',...
        'markersize',6,'linewidth',2);
    h2 = plot(obsHW.timeHW,obsHW.HW,'x','color',blue,...
        'markersize',6,'linewidth',2);
    l1 = plot(predLW.timeLW,predLW.LW,'k+',...
        'markersize',7,'linewidth',1.5);
    l2 = plot(obsLW.timeLW,obsLW.LW,'+','color',blue,...
        'markersize',7,'linewidth',1.5);
    s = surgeHW(surgeHW.obsTIME > xstart & surgeHW.obsTIME < xend,:);
    for i=1:height(s)
        text(s.obsTIME(i)+1/12,s.obs(i),...
            sprintf('%4.3f',s.surge(i)),...
            'color',darkred,'fontsize',9);
    end
    s = surgeLW(surgeLW.obsTIME > xstart & surgeLW.obsTIME < xend,:);
    for j=1:height(s)
        text(s.obsTIME(j)+1/12,s.obs(j),...
            sprintf('%4.3f',s.surge(j)),...
            'color',darkred,'fontsize',9);
    end
    xlim([xstart xend]);
    ylim([ymin ymax]);
    xlabel(year);
    ylabel('Water level (mAHD)');
    xticks(xstart:1:xend);
    xtickformat('dd-MMM')
    % datetick('x','dd-mmm','keeplimits','keepticks');
    title('8. Skew surges');
    legend([o p h1 h2 l1 l2],'Observed','Predicted','Predicted HW','Observed HW',...
        'Predicted LW','Observed LW',...
        'location','eastoutside');
    box on
    hold off

export_fig('Hobart_8-skewsurge.png','-r300','-dpng','-transparent','-painters')