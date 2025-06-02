% DD system Train
clc;close all;clear;
addpath('Plot\')
addpath('Dsp\')
addpath('Sync\')
addpath('Phase_Sync\')

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
skip = 0.1 * Tx.TxPHY.NSym; % 跳过初始瞬态（20%符号）

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
%% 经过信道
hch = [0.207, 0.815, 0.207];
sigRx=Tx.channelApply(hch,signal);
matchOutO=Tx.matchFiltering(sigRx);
% downsample
outSignal=downsample(matchOutO,sps);
% decode
clockRecovery.PAM_ExecuteDecoding(outSignal);
%% LUT 计算
% 计算LUT误差error
err = (signal - sigRx).';
% 记忆深度(考虑前3符号)
mem_len=3;
% LUT每次遍历三个符号
N=3;
% 总共模式数
L = M^N;
%模式索引
Idx=1:L;
% 将索引生成矩阵形式
Idx = reshape(Idx, (repmat(M, 1, N)));
% 将模式数组进行转置
for i=1:(N+1)
    Idx(:,:,i)=Idx(:,:,i).';
end

% 处理信号
% refSig=sigTxo-mean(sigTxo);
% refSig=pnorm(refSig);
refSig=signal;
% 得到的矩阵是，发送端每个符号，与参考符号的差值，每个符号对应一列
k=abs(bsxfun(@minus, refSig, const*Ksym));
% 符号到最近参考电平的映射
[~, sig_idx]=min(k);
% 滑动窗口提取历史符号模式(N=3);以发射信号为中心构建表格
sig_rwin = rolling_window_central(sig_idx, N, 'true');
sig_rwin=sig_rwin.';
% sig_rwin = rolling_window(sig_idx, N, 'true');% 以发射信号为第一个，后续需要移位
% sig_rwin=sig_rwin.';
% shift=floor(mem_len/2);
% pattern_idx=circshift(pattern_idx,shift);

% 计算模式索引(1-64)
pattern_idx=zeros(length(sig_rwin),1);
for i=1:length(sig_rwin)
    pattern_idx(i) = Idx(sig_rwin(2,i),sig_rwin(3,i),sig_rwin(1,i));
end
% 计算各模式平均误差(LUT)
ea = cal_lut_avg(err, pattern_idx, pattern_idx, L);
Error=ea(pattern_idx).';