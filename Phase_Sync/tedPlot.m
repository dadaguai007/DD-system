
% MLTED, ELTED, ZCTED, GTED, or MMTED.
%%
clc;close all;clear;
addpath("Phase_Sync\")
TED='MMTED';
plotTedGain(TED)
%%
clc;close all;clear;
addpath("Phase_Sync\")
rolloff = 0.2;
plotSCurve = 1;
TED='MMTED';
Kp_analytic_ml = calcTedKp(TED, rolloff, 'analytic', plotSCurve);