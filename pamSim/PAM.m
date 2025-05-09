% PAM 通信系统
clc;clear;close all;
addpath('\001-处理中\相干算法\optical_communication\')
addpath('\001-处理中\相干算法\optical_communication\DFE-volterra\')
% 信号参数
SpS = 8;
Rs  = 40e9;
Ts  = 1/Rs ;
Fs  = SpS*Rs;
Ta  = 1/Fs;

% 调制器参数
Vpi = 10;
Vb = -Vpi/2;

% 输入光功率
Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W

% fiber
param=struct();
param.Ltotal = 40; %km
param.Lspan =20;
param.hz= 0.1;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='none';
param.Fs=Fs;

% PD接收
paramPD=struct();
paramPD.B =Rs;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=Fs;

rng(123);
%PAM
M=4;
data_2bit=randi([0,1],log2(M),80000);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;
% Mapeia bits para pulsos
symbTx = pammod(symbols,M,0,'gray');
% symbTx=real(symbTx);
symbTx = pnorm(symbTx);

% 参考向量
label=symbTx;
% Upsampling
symbolsUp = upsample(symbTx, SpS);

% Puls
hsqrt = rcosdesign(0.01,256,SpS,'sqrt');
% pulse shaping
sigTx=conv(symbolsUp,hsqrt,'same');

% 对信号进行重复
k=2;
sigTx=repmat(sigTx,1,k);
% 信号放大
PAM_amptype;

% mod index
m=Modulation_index(Sig_Tx,Vpi,'pam');
fprintf(' the module index =%.3f \n', m);

%mzm
Ai= sqrt(Pi);
sigTxo = mzm(Ai, Sig_Tx, Vb,Vpi);
power=signalpower(sigTxo);
fprintf(' after module signal power: %.2f dBm\n', 10 * log10(power / 1e-3));

plot_spectrum(sigTxo,Fs);
type='none';
if strcmp(type,'Tran')
% PAM_tran
sigRxo=ssfm(sigTxo,param);
power2=signalpower(sigRxo);
fprintf(' after ssfm signal power: %.2f dBm\n', 10 * log10(power2 / 1e-3));
else
    sigRxo=sigTxo;
end
% 
% 增加信道响应
PAM_Channel;

%pd
ipd = pd(sigRxo, paramPD);

% match
sigRx_E=ipd.';
sigRx_E=sigRx_E-mean(sigRx_E);
sigRx_E=pnorm(sigRx_E);

% 同步
% Label=repelem(label,SpS); % 别说，使用这样的方法进行同步，效果更为准确！！！
% [sigRx_E,ref_sync,ff_ssfm] = sync(sigRx_E,Label);
if 1
% 均衡
PAM_Equ;
else
% 下采样
sigRx_E=downsample(sigRx_E,SpS);
end
figure;hold on;
plot(sigRx_E,LineWidth=1)
eyediagram(sigRx_E(1000:10000),2*SpS)
PAM_Decode;