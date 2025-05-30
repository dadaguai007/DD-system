function [ resampledData ] = ResampleInterp( rawData, tSeqBefResamp,  tSeqAftResamp )
%RESAMPLEINTERP Summary of this function goes here
% Detailed explanation goes here
if max(tSeqAftResamp) > max(tSeqBefResamp) || min(tSeqAftResamp) < min(tSeqBefResamp)
    error('insufficient input rawdata');
end

M = 200000;      % maximum resampled points at one time
NAft = numel(tSeqAftResamp);
resampEndPoits = (1:M:NAft).';
if resampEndPoits(end) < NAft
    resampEndPoits = [resampEndPoits; NAft];
end
resampledData = zeros(size(tSeqAftResamp));
for resampIndex = 1:numel(resampEndPoits)-1
    aftResampIndexL = resampEndPoits(resampIndex);
    aftResampIndexR = resampEndPoits(resampIndex+1);
    [temp befResampIndexL] = min(tSeqBefResamp <= tSeqAftResamp(aftResampIndexL));
    befResampIndexL = befResampIndexL - 1;
    [temp befResampIndexR] = max(tSeqBefResamp >= tSeqAftResamp(aftResampIndexR));
    if tSeqBefResamp(befResampIndexL) > tSeqAftResamp(aftResampIndexL) ...
            || tSeqBefResamp(befResampIndexR) < tSeqAftResamp(aftResampIndexR)
        error('Time for resampling out of range!');
    end    
    resampledData(aftResampIndexL:aftResampIndexR) = interp1(tSeqBefResamp(befResampIndexL:befResampIndexR), ...
        rawData(befResampIndexL:befResampIndexR), tSeqAftResamp(aftResampIndexL:aftResampIndexR),'spline');
end
