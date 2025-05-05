% DRE test


BitRateDefault =124e9;
FrmLen = 2^17;
SampleRateDefault = 2 * BitRateDefault;
UpRate = SampleRateDefault/BitRateDefault;
AWGSamplingRate = 224e9;
TrainLen = 1024;
N_T = 128; 
RollOff =0.01;
BitPerSym=2;
DREmode=0;
%% upsampling and RC filtering
TxWfm = upsample(TxSym,UpRate);
RRCFilt = rcosfir(RollOff, N_T, UpRate, 1, 'sqrt');
TxWfm = conv(TxWfm, RRCFilt, 'same');
% 归一化
TxWfm = TxWfm/sqrt(mean(abs(TxWfm).^2));

% 保留滤波器主瓣，减少旁瓣影响。
% 降低后续计算复杂度（滤波器长度从257点缩减为41点）。
h_DRE= RRCFilt(N_T*UpRate+1-20:N_T*UpRate+1+20);

%% Resampling
% h_DRE = RRCFilt;
% 将滤波器采样率从原值（如248GS/s）降低至80%（198.4GS/s）。
%适配目标系统的等效噪声带宽（如DAC模拟带宽限制）。
h_DRE = resample(h_DRE,SampleRateDefault*0.8,SampleRateDefault);
TxWfm = resample(TxWfm,AWGSamplingRate,SampleRateDefault);
SampleRateDefault = AWGSamplingRate;
% DRE 使用
switch DREmode
    case 0
        TxSig = TxWfm;
    case 1
        % % 基础DRE模式
        QAM_order = 16;
        h_DRE= 1;
        [sigout] = Quant_DRE((TxWfm),h_DRE,QAM_order);
        TxSig = sigout.';
        scatterplot(TxSig)
    case 2
        % 频响优化模式
        QAM_order = 16;
%         h_DRE= h_DRE(N_T*UpRate+1-20:N_T*UpRate+1+20);
%         h_DRE= Response_Tdomain;
        PointNum=length(h_DRE);
        X_axis=linspace(-SampleRateDefault/1e9/2,SampleRateDefault/1e9/2,PointNum);
        Y_axis=10*log10(abs(fftshift((fft(h_DRE,PointNum)))));
        figure;
        plot(X_axis,Y_axis);
        title('h_DRE');
        [sigout] = Quant_DRE((TxWfm),h_DRE,QAM_order);
        TxSig = sigout.';
end