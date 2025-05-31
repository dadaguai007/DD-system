% 频域MMSE均衡器 得到权重向量
function [w]=getMMSEFDEW(channel,NFFT,sigma_noise)
% 计算得到频域权重
    w = 1./(fft(channel, NFFT) + (sigma_noise)^2); % implements the formula dicussed in the report
end