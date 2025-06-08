% FTN Signal Test
clc;clear;close all;

% 数据符号
nSymbols=50;
rng(1001);
M=4;
symbols_rand=randi([0,1],log2(M),nSymbols);
data_symbols=qammod(symbols_rand,M,'InputType','bit','UnitAveragePower',1) ;


% FTN 参数
ftnParam = 0.6;         % FTN压缩因子τ (0 < τ < 1)
sps_nyquist = 50;       % 奈奎斯特速率下的每符号采样数
sps_ftn = ftnParam * sps_nyquist;   % FTN速率下的每符号采样数(τ×Nyquist)
rollOff = 0.5;          % SRRC滚降因子
NISI = 10;              % ISI干扰长度

% 上采样
data_symbols=upsample(data_symbols,sps_ftn);
nSymbols=length(data_symbols);


% 生成平方根升余弦(SRRC)滤波器
% 滤波器长度：2*NISI + 1个符号周期
g = rcosdesign(rollOff, 2*NISI, sps_nyquist);

% 计算滤波器自相关函数
g_xcorr = xcorr(g);

% FTN速率下采样自相关函数(脉冲的一半，进行下采样);
% 滤波器的正常长度为L
g_xcorr_samp = g_xcorr((length(g_xcorr)+1)/2 : sps_ftn: end);

% 添加边界符号 + 数据符号
% ISI影响范围由脉冲响应长度 L决定；每个数据符号会受前后 L-1个符号的影响
%为准确模拟连续传输，需确保：边界长度 > L-1

nSymbolsExtraHalf = 0; % 可以为比较小的数
nSymbolsExtraHalf=nSymbolsExtraHalf+length(g_xcorr_samp);
dummy_pre = zeros(1 ,nSymbolsExtraHalf);
dummy_post = zeros( 1,nSymbolsExtraHalf);

% 生成ISI-信道矩阵；保持数据长度一致
G = toeplitz([g_xcorr_samp, zeros(1, nSymbols + 2*nSymbolsExtraHalf - length(g_xcorr_samp))]);


symbols = [dummy_pre, data_symbols, dummy_post];
symbols=symbols.';
% 通过ISI信道矩阵生成FTN信号
txSymbols = G * symbols;
txSymbols1=filter(g_xcorr_samp,1,symbols);
% AWGN 加入
% 接收端裁剪;除去边界
rxSymbols = txSymbols(nSymbolsExtraHalf+1 : end-nSymbolsExtraHalf);
rxSymbols=downsample(rxSymbols,sps_ftn);
% MAP 接收
llr = qamdemod(rxSymbols, M, 'gray', 'OutputType', 'llr');              % Estimate LLRs
rxBits = double(llr < 0);
% 转换为列向量
label =symbols_rand (:);