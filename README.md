# AtchafalayaTideHt

Predict water level over the short-term in Atchafalaya Delta, Louisiana.

## Description

Use hourly water level data from the [Coastwide Reference Monitoring System](https://cims.coastal.la.gov/DataDownload/DataDownload.aspx?type=hydro_hourly) to predict future water levels with a harmonic oscillation model. This script also uses wind speed and direction, which may have an influence on water level at some sites in coastal Louisiana.

## Getting Started

### Dependencies

* Install UTide (Unified Tidal Analysis and Prediction Functions) to get the functions ut_reconstr and ut_solv. Add to the MATLAB path by running the file (will produce error: not enough output arguments or similar, can ignore) or [addpath()](https://www.mathworks.com/help/matlab/ref/addpath.html). Download available in [MATLAB File Exchange.](https://www.mathworks.com/matlabcentral/fileexchange/46523-utide-unified-tidal-analysis-and-prediction-functions)

### Executing program

* Run entire script at once.
