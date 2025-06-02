% DD system Train
clc;close all;clear;
addpath('Plot\')
addpath('Dsp\')
addpath('Sync\')
addpath('Phase_Sync\')
addpath('THP\')

%% 系统初始化
% 信号生成
ddGeneration;

% 时间轴和频率轴创建
[freq,time]   =  freq_time_set(Tx.TxPHY.NSym *Tx.TxPHY.sps,fs);

% 信号生成
[signal,pamsignal]=Tx.dataOutput();

% DSP参考信号
label=pamsignal;
% 参考星座图
[const,Ksym]=Tx.creatReferenceConstellation();

% 数据性能起始位置
skip = 0.2 * Tx.TxPHY.NSym; % 跳过初始瞬态（20%符号）


if 0
    % THP
    THP=THPClase(1,1);
    load('THP_w.mat');
    modN=16;
    % Pre
    [pre_equalised_with_mod,xInput]=THP.preTHP((pamsignal),tapsTHP,modN);
    %     % 归一化
    %     preSignalNorm = pnorm(pre_equalised_with_mod);
    % 脉冲成型 装载成型后的信号
    signal=Tx.applyShapingFilter((pre_equalised_with_mod),hsqrt);
    signal=signal.';
end


if 0
    % 预编码
    % PRS code
    % linear
    D=0.5;
    taps=4;
    [filteredPRS,prs_sig,lcoeff]=Tx.prsSignal(D,taps);

    % nonlinear
    taps1=5;
    taps2=3;
    taps3=1;
    [filteredPRS_Nonlinear,prs_sig_Nonlinear,w]=Tx.prsSignalNonlinear(D,taps1,taps2,taps3);

    % 色散预补偿 生成FIR
    FiberLen=param.Ltotal*1e3;
    DER1=1;
    signalPreEDCFIR=Tx.cdFIRApply(DER1,FiberLen);
    % 色散预补偿 GS 算法
    DER2=22;
    signalPreEDCGS=Tx.cdGSApply(DER2,FiberLen);
end
%% Tx and Rx Lo
% phase noise TX
lw   =  0e3;
phase_Noise  =  Tx.phaseNoise(signal,lw);

%% 器件频响建立

% obj.Implementation.responType='Bessel';
% f3dB=20e9;
% order=5;
% verbose=0;
% [filt,H]=Rx.createFrequencyResponse(freq,order,f3dB,verbose);
% out1 = firFilter(filt.h, txSig); % 第一种滤波器工作方式,时域
% out2=filter_work(txSig,H); % 第二种滤波器工作方式，频域


%% 发射机创建
amp_factor=0.28; % EA放大
sigTx=signal*amp_factor;

% mod index
m=Modulation_index(sigTx,Vpi,'pam');
fprintf(' the module index =%.3f \n', m);

%mzm
sigTxo = mzm(phase_Noise, sigTx, Vb,Vpi);

% 设置功率
sigTxo=Tx.setSignalPower(sigTxo,channelPowerType,Pout_dBm);
power=signalpower(sigTxo);
fprintf(' optical signal power: %.2f dBm\n', 10 * log10(power / 1e-3));


%% 添加信道延迟
% timeOffset = 25;       % 信道延迟（采样点数）
% % 延迟模块（模拟信道延迟）
% DELAY = dsp.Delay(timeOffset);
% % 信道延迟
% delaySig = step(DELAY, sigTxo);    % 添加固定延迟;直接添加零


% if 1
%     Tx.Button.delaySignal='symbol';
%     delay=20;
% else
%     Tx.Button.delaySignal='time';
%     delay=50e-12;
% end
% delaySig=Tx.delaySignal(sigTxo,delay);

%% 信号传输
sigRxo=ssfm(sigTxo,param);
power2=signalpower(sigRxo);
fprintf(' after ssfm signal power: %.2f dBm\n', 10 * log10(power2 / 1e-3));
% 频谱参考
mon_ESA(sigRxo,fs);
% 接收
sigRxE = pd(sigRxo, paramPD);
%% 应用集成的时钟偏移函数
% ppm=-0.5;
% jitter_rms=0; %1e-10
% txResamp1=Tx.addSamplingClockOffset(txSig,ppm,jitter_rms);

%% 匹配滤波
% 滤波器赋值
Tx.TxPHY.hsqrt=hsqrt;
matchOut=Tx.matchFiltering(sigRxE);

%% 定时恢复算法
% 参数初始化
Bn_Ts    = 0.01;       % 环路带宽×符号周期（归一化噪声带宽）
eta      = 1;          % 环路阻尼系数（控制收敛速度）
debug_tl_static  = 1; % 静态调试标志（1=显示最终星座图）
debug_tl_runtime = 0; % 运行时调试标志（1=显示同步过程示波器）

clockRecovery = DspSyncDecoding( ...
    fs,...         % 接收信号的采样率
    fb, ...        % 接收信号的波特率
    M,...         % 接收信号的格式
    fs/fb, ...     % 上采样率
    2*fs/fb,...         % 时钟信号的上采样率
    skip, ...      % 误码计算起始位置
    label,...      % 参考信号
    "MLTED");
clockRecovery.Implementation.eta       = eta;              % 环路阻尼因子（稳定性控制）
clockRecovery.Implementation.Bn_Ts     = Bn_Ts;           % 环路带宽 × 符号周期（控制同步速度）
clockRecovery.Implementation.Ex =1;
clockRecovery.Implementation.intpl=2;          % 插值方法：0=多相，1=线性，2=二次，3=三次
rollOff=Tx.Nr.psfRollOff;
rcDelay=Tx.Nr.psfLength ;
sps=Tx.TxPHY.sps;
% [ rxSync1 ] = symbolTimingSync(TED, intpl, sps, sig_RxE, matchOut, K1, K2, ...
%     const, Ksym, rollOff, rcDelay, debug_tl_static, debug_tl_runtime);
% 集成的算法
% %接收信号，为sig_RxE\txResamp1
% [rxSync,Kp] =clockRecovery.TED_recoverOptimalSamplingPoints(sig_RxE,matchOut,rollOff,rcDelay,const,Ksym);

% 搜寻时间信号的极值
% [rxSync,P,OptSampPhase,MaxCorrIndex]=clockRecovery.time_phase_Recovery(matchOut);
%% LUT 计算
if 0
    hch = [0.207, 0.815, 0.207];
    sigRx=Tx.channelApply(hch,signal);
    matchOutO=Tx.matchFiltering(sigRx);
    % downsample
    outSignal=downsample(matchOutO,sps);
    % decode
    clockRecovery.PAM_ExecuteDecoding(outSignal);
    % 记忆深度(考虑前3符号)
    mem_len=3;
    % LUT 查找表建立
    [LUT, idx_I, idx_Q] = cal_lut(signal, sigRx, const*Ksym, mem_len);
    % 应用LUT
    sigTxoLUT=sigTx-LUT;
    % 过信道
    preLUTRxo=Tx.channelApply(hch,sigTxoLUT);
    % match
    matchOutLUT=Tx.matchFiltering(preLUTRxo);
    % downsample
    outSignal=downsample(matchOutLUT,sps);
    % decode
    clockRecovery.PAM_ExecuteDecoding(outSignal);
end
%% Eq
EQ=struct();
EQ.u=0.001;
EQ.k1=31;
EQ.k2=15;
EQ.ref=16;
EQ.sps=2;
EQ.lamda=0.9999;
EQ.delta=0.01;

refEq=Tx.createReferenceSignal(label);
% FFE_LMS
[yEq,en,w] = FFE_LMS(EQ, matchOut.', refEq.');
% DFE_LMS
% [yEq,en,w] = DFE_LMS(EQ,  matchOut.', refEq.');
plotEquParam(matchOut,yEq,refEq,w,en)
if 0
    yEq = THP.modulo(yEq, modN);  % 接收端模运算恢复信号
end
% decode
[decodedDataEq,berEq]=clockRecovery.PAM_ExecuteDecoding(yEq);
% Tx.PAM_ExecuteDecoding(yEq,const);
if 0
    tapsTHP=w(EQ.k1+1:end);
    save('THP_w', 'tapsTHP');
end

if 0
    % pnc and pec
    alpha_pec = 0.42;
    alpha_pnc = -0.5;
    alpha_snd = 0.3;
    alpha_snd2 = -0.15;
    alpha_snd3 = 0;
    pec_out = my_pec(yEq, alpha_pec, const);
    pec_out_2 = my_pec_advanced(yEq.', [alpha_snd, alpha_snd2], const);
    pec_out_3 = my_pec_advanced(yEq.', [alpha_snd, alpha_snd2, alpha_snd3], const);

    pnc_out = my_pnc(yEq, alpha_pnc, const, 2);
    pnc_out_2 = my_pnc_enhanced(yEq, alpha_pnc, const);
    % decode
    clockRecovery.PAM_ExecuteDecoding(pnc_out_2);
end

if 0
    % 噪声白化
    noise = yEq - refEq(1:length(yEq));
    lambda=0.9999;
    switch log2(M)
        case 1 % OOK
            Coe = zeros(1,8);
            step = 0.005;
        case 2 % PAM
            Coe = zeros(1,4);
            step = 0.005;
    end
    % 噪声白化滤波器
    % 失效,有待探究
    % [Shaping_Tap, error] = Modified_PR_PostFilter_RLS_Train( yEq, label(1:2e3), lambda, Coe )
    % [Shaping_Tap, error] = Modified_APR_PostFilter_RLS_Train( yEq, label(1:2e3),  lambda, Coe );


    % [Shaping_Tap, error] = PR_PostFilter_RLS_Train( yEq, label(1:2e3), lambda, Coe );
    % [Shaping_Tap, error] = Modified_APR_PostFilter_LMS_Train( yEq, label(1:2e3), step, Coe );
    [Shaping_Tap, error] = Modified_PR_PostFilter_LMS_Train(yEq, label(1:2e3), step, Coe );
    % [Shaping_Tap, error] = PR_PostFilter_LMS_Train( yEq, label(1:2e3), step, Coe );
    yEq_NoiseWhiting =filter(Shaping_Tap,1,yEq); % 对信号白化
    noise_NoiseWhiting  =filter(Shaping_Tap,1,noise); % 对噪声白化，查看效果
    RxDown =mlseeq(yEq_NoiseWhiting,Shaping_Tap,const,1000,'rst',1);%mlse

    % decode
    clockRecovery.PAM_ExecuteDecoding(real(RxDown));
end
%% 无Eq 解码

% downsample
outSignal=downsample(matchOut,sps);

if 0
    % THP
    outSignal = THP.modulo(outSignal, modN);  % 接收端模运算恢复信号
end

% decode
clockRecovery.Button.Train='off';
[decodedData,ber]=clockRecovery.PAM_ExecuteDecoding(outSignal);
% Tx.PAM_ExecuteDecoding(outSignal,const);

% 实部的分布
[value,percentReal] = amplitude_distribution(yEq);


% FiberLen=param.Ltotal*1e3;
% Tx.cdChannelResponse(FiberLen);