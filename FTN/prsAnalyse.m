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
[~,time] = freq_time_set(Tx.TxPHY.NSym *Tx.TxPHY.sps,fs);
% sps
sps=Tx.TxPHY.sps;

% 信号生成
[signal,pamsignal,pulse]=Tx.dataOutput();
% DSP参考信号
label=real(pamsignal);
% PAM Signal
pamsignal=pnorm(pamsignal);

% 参考星座图
[const,Ksym]=Tx.creatReferenceConstellation();

% 数据性能起始位置
skip = 0.1 * Tx.TxPHY.NSym; % 跳过初始瞬态（20%符号）
Tx.Nr.ncut_index=skip;
refEq=Tx.createReferenceSignal(label);

% 绘图时间轴
index=skip:skip+100;

% PRS code
% linear
D=1;
taps=3;
[prsPulse,prs_sig,lcoeff]=Tx.prsSignal(D,taps);

%% 匹配滤波
% 滤波器赋值
Tx.TxPHY.hsqrt=hsqrt;
matchOut=Tx.matchFiltering(prsPulse);
dataout=downsample(matchOut,sps);
%% Eq
EQ=struct();
EQ.u=0.005;
EQ.k1=21;
EQ.k2=15;
EQ.ref=11;
EQ.sps=2;
EQ.lamda=0.9999;
EQ.delta=0.01;

% FFE_LMS
[yEq,en,w] = FFE_LMS(EQ, matchOut.', refEq.');
plotEquParam(matchOut,yEq,refEq,w,en)
Tx.Button.Train='on';
Tx.PAM_QuantifyDecoding(yEq);
% 针对PRS的接收端(PR的解码)
yout = mlseeq(yEq, lcoeff, const, 1000, 'rst',1);
yout=real(yout);
Tx.PAM_QuantifyDecoding(yout);