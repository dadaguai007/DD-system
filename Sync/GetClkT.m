function [ clkT stdOfSampledClk ] = GetClkT( clkRawData, osciClkT, clkTRough, maxSampleTimeError )
%GETCLKT Summary of this function goes here
%   Detailed explanation goes here
NClkRaw = numel(clkRawData);
% approximate number of points after resample
NClkResampRough = NClkRaw * osciClkT / clkTRough;
% difference between adjacent test clock periods
deltaClkT = clkTRough * maxSampleTimeError / NClkResampRough;
NHalfClkTVec = ceil(0.1/2/maxSampleTimeError);
% vector containing test clock periods
clkTVec = clkTRough + (-NHalfClkTVec:NHalfClkTVec).' * deltaClkT;

% time sequence before resample
tSeqBefSamp = (0:NClkRaw-1).' * osciClkT;
maxt = max(tSeqBefSamp);
stdOfSampledClkVec = zeros(size(clkTVec));
for clkTIndex = 1:numel(clkTVec)
    tSeqAftSamp = (0/5 * clkTVec(clkTIndex):clkTVec(clkTIndex):maxt).';
    resampledData = ResampleInterp( clkRawData, tSeqBefSamp, tSeqAftSamp);
    stdOfSampledClkVec(clkTIndex) = std(resampledData);
end

% Find clock period corresponding to minimum std after its resampling
[valueMinStd minStdIndex] = min(stdOfSampledClkVec);
if minStdIndex == 1 || minStdIndex == numel(stdOfSampledClkVec)
    error('Real clock period out of searching range!');
end
polyForStd = polyfit((-1:1).',stdOfSampledClkVec(minStdIndex-1:minStdIndex+1), 2);
minStdIndexMinor = roots(polyder(polyForStd));
clkT = polyval(polyfit((-1:1).',clkTVec(minStdIndex-1:minStdIndex+1), 2),minStdIndexMinor);
% calculate the corresponding minimum std
tSeqAftSamp = (0/5 * clkT:clkT:maxt).';
resampledData = ResampleInterp( clkRawData, tSeqBefSamp, tSeqAftSamp);
stdOfSampledClk = std(resampledData);