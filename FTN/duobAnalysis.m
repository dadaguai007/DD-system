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
% pamsignal=pnorm(pamsignal);

% 参考星座图
[const,Ksym]=Tx.creatReferenceConstellation();

% 数据性能起始位置
skip = 0.1 * Tx.TxPHY.NSym; % 跳过初始瞬态（20%符号）
Tx.Nr.ncut_index=skip;
refEq=Tx.createReferenceSignal(label);

% 绘图时间轴
index=skip:skip+100;

% 1+D PR滤波器
reverse='False';
rolloff=0.01;
psfLength=16;
pulseDoub = Duob(sps, psfLength, reverse, rolloff);
symbolsUp = upsample(pamsignal, sps);
symbolsUp=pnorm(symbolsUp);
prsPulse=conv(symbolsUp,pulseDoub,'same');

% prs
taps=1;
pulsePrs=getPrsPulse(sps, psfLength, reverse, rolloff,taps);
prsPulse1=conv(symbolsUp,pulsePrs,'same');


%眼图
% eyediagram(prsPulse1(1:1e4),4*sps)

%% 匹配滤波
% 滤波器赋值
pulseDoub_match = Duob(sps, psfLength, 'true', rolloff);
Tx.TxPHY.hsqrt=pulseDoub_match;
matchOut=Tx.matchFiltering(prsPulse);

%% Eq
EQ=struct();
EQ.u=0.003;
EQ.k1=31;
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
yout = mlseeq(yEq, [1,1], const, 100, 'rst',1);
yout=real(yout);
Tx.PAM_QuantifyDecoding(yout);




%% 预编码 双二进制
% clc;clear;close all;
% SpS = 2;
% Rs  = 32e9;
% Ts  = 1/Rs;
% Fs  = 1/(Ts/SpS);
% Ta  = 1/Fs;
% rolloff = 0.01;
% 
% M=2;
% % rand
% bits = randi([0, 1], 1, 10000*log2(M));
% c_bits=zeros(size(bits,1));
% for i = 1:length(bits)
%     if i==1
%         c_bits(i)=bits(i);
%     else
%         %异或运算
%         c_bits(i) =xor(bits(i),c_bits(i-1));
%     end
% end
% c_bits=reshape(c_bits,log2(M),[]);
% 
% symbols = 2.^(0:log2(M)-1)*c_bits;
% symbTx = pammod(symbols,M,0,'gray');
% symbolsUp = upsample(symbTx, SpS);
% 
% reverse='False';
% pulseDoub = Duob(SpS, 16, reverse, rolloff);
% pulse = rcosdesign(rolloff,2048,SpS,'sqrt');
% % sigTx=filter(pulseDoub,1,symbolsUp);
% sigTx=conv(symbolsUp,pulseDoub,'same');
% t = (0:length(symbTx)-1) * (Ta / 1e-9);
% 
% %眼图
% eyediagram(sigTx,4*SpS)