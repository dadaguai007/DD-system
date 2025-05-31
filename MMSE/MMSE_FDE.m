% MMSE 频域均衡器

function w=MMSE_FDE(channel, snr_db,NFFT)
snr_pow = 10^(snr_db/10);
sigma_noise = 1/sqrt(snr_pow);
w=getMMSEFDEW(channel,NFFT,sigma_noise);
end
