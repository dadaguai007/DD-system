function estFreq = KaysFrequencyEstimator(vX,samplingFreq)
% Refer to S. Kay's paper: A fast and accurate single frequency estimator
% in IEEE Trans. on  Acoustic
% See also: https://dsp.stackexchange.com/questions/76644/
% simple-and-effective-method-to-estimate-the-frequency-of-a-single-sine-signal-in

numSamples = length(vX);
estFreq = 0;
for ii = 1:(numSamples - 1)
    estFreq = estFreq + angle(vX(ii)' * vX(ii + 1));
end

estFreq = estFreq / (2 * pi * (numSamples - 1));
estFreq = samplingFreq * estFreq;