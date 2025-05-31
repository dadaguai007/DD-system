% Test MMSE equalizer
disp(datestr(now));
NFFT = 1000; % num of FFT or num of symbols in each frame
num_taps = 10;
channnel=[0.74 -0.514 0.37 0.216 0.062];
snr_db=(0:2:20);
%Add CP
cyclicPrefixes = [1 3 5 10];
for i=1:length(cyclicPrefixes)
cyclicPrefix=cyclicPrefixes(i);
end
if strcmp(type,'FDE') % addition of the cyclic prefix to enable circular convolution operation
    trSymVec = [trSymVec(end-cyclicPrefix+1:end) trSymVec];
end
% %%%%%%%%%%%%% CHANNEL %%%%%%%%%%%%%%%%%%%%
                % conv 的去延时，从信号尾部进行排除
%             channel_result = conv(trSymVec, channel);
%             channel_result = channel_result(1:end-channel_length+1);
%             trSymVec = channel_result;
%             noise=1/sqrt(2)*[randn(1, length(trSymVec)) + 1j*randn(1,length(trSymVec))];
%             recSigVec=trSymVec+sigma_noise*noise;
if strcmp(type,'FDE')
    equalizer=MMSE_FDE(channel, snr_db,NFFT);
    recSigVec = recSigVec(1+cyclicPrefix: end); % deletion of cyclic prefix
    recSigVec = fft(recSigVec, NFFT); % taking FFT of incoming signal
    recSigVec = recSigVec.*equalizer; % equalizing with multiplication
    recSigVec = ifft(recSigVec, NFFT); % IFFT
elseif strcmp(type,'TDE')
    equalizer = MMSE_TDE(num_taps, channel); % get time domain MMSE 10-tap equalizer
    equalizer_result = conv(recSigVec, equalizer);
    equalizer_result = equalizer_result(1:end-equalizer_length+1); % get causal conv output
    recSigVec = equalizer_result;
end