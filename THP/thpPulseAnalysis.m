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
[~,time]   =  freq_time_set(Tx.TxPHY.NSym *Tx.TxPHY.sps,fs);
% sps
sps=Tx.TxPHY.sps;
% 绘图时间轴
index=1:200;

% 信号生成
[signal,pamsignal,pulse]=Tx.dataOutput();

% PAM Signal 
pamsignal=real(pamsignal);
% DSP参考信号
label=real(pamsignal);
% 参考星座图
[const,Ksym]=Tx.creatReferenceConstellation();

% 数据性能起始位置
skip = 0.3 * Tx.TxPHY.NSym; % 跳过初始瞬态（20%符号）
Tx.Nr.ncut_index=skip;
refEq=Tx.createReferenceSignal(label);


% 应用脉冲成型
% signal=Tx.applyShapingFilter(pamsignal,pulse);
%% 信道创建
syms z;                     % 声明符号变量z（用于Z域分析）

c1=1;
c2=-0.4;
c3=0.7;
channel = c1+c2*z^{-1}+c3*z^{-2};       % 定义信道传递函数：H(z) = 1 - z⁻¹（典型ISI信道模型）
% channel = 1-z^{-1}; 
% channel=1-z^{-1}+3*z^{-2}+0.1*z^{-10};
% 创建类
THP=THPClase(1, 1);
N=3;  % 2倍的（M-1）

% 定义信道响应系数
b = [c1,c2,c3];  % 1,0.5,0.7
% b = [1,-1];% H(z) = 1 + z⁻¹ 
a = 1;  % 分母系数 (FIR系统)
%% 分析信道

if 0

% FFT 点数
N1=4096;
% 计算频率响应
[H, freq] = freqz(b, a, N1, 'whole', fs);
% 转换为单边频谱
H_mag = abs(H(1:N1/2+1));
H_phase = angle(H(1:N1/2+1));
freq = freq(1:N1/2+1);
figure('Position', [100, 100, 900, 700])

% 1. 幅度响应
subplot(2,1,1)
plot(freq, 20*log10(H_mag)), grid on
title('信道幅度响应 |H(f)| (dB)')
xlabel('频率 (Hz)')
ylabel('增益 (dB)')
xlim([0, fs/2])
ylim([-10, 20])

% 2. 相位响应
subplot(2,1,2)
plot(freq, unwrap(H_phase)*180/pi), grid on
title('信道相位响应 \angle H(f)')
xlabel('频率 (Hz)')
ylabel('相位 (度)')
xlim([0, fs/2])

end

%% without THP Code Data Train
% % Z Tran
% z_x=THP.Transz(signal);
% 
% % 信号通过信道传输
% z_output = z_x * channel;   % 信道输出Y(z) = X(z)·H(z)
% output_without_precoder=THP.zTrans(z_output,signal);% 逆Z变换得到时域表达式

% 时域冲击响应
output_without_precoder=filter(b,a,signal);


%%  THP Code Data Train
% 将时域信号转换为Z域表达式
z_x=THP.Transz(signal);

% 信道逆预均衡
channel_inverse = 1/channel; % 计算信道逆
z_pre_equalised = z_x * channel_inverse; % 预均衡：X'(z) = X(z)·H⁻¹(z)

% 预均衡后的信号
pre_equalised=THP.zTrans(z_pre_equalised,signal);% 逆Z变换得到时域表达式

% 模运算控制信号幅度
pre_equalised_with_mod = THP.modulo(pre_equalised, N); % 将信号限制在[-N/2, N/2)

% 经过信道
% % 将预编码信号重新转换为Z域
% zz_pre_equalised=THP.Transz(pre_equalised_with_mod);
% 
% % 预编码信号通过实际信道
% z_output = zz_pre_equalised * channel; % 
% output=THP.zTrans(z_output,signal); % 接收端时域信号


output=filter(b,a,(pre_equalised_with_mod));
output_with_mod = THP.modulo(output, N);  % 接收端模运算恢复信号



%% CCDM
sz_window=10;
[~,ccdfy1]=ccdf(signal,sz_window,0);
[ccdfx,ccdfy2]=ccdf(pre_equalised_with_mod,sz_window,0);

color=distinguishable_colors(20);
marker = 'so^d>v*phx';
CCDF=[ccdfy1;ccdfy2];

figure;
for i=1:2
semilogy(ccdfx,CCDF(i,:),LineWidth=1.25,Color=color(i,:),Marker=marker(i));
hold on;
end
legendArrary={'Without THP','With THP'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('','PAPR_0 (dB)','Pr(PAPR > PAPR_0)',[0 8],[10^-3 10^0],legendArrary,flag,FontSize)




%% Plot
n = 0:length(signal)-1; % 时间轴
figure
stem(n(index), [signal(index).', output_without_precoder(index).'],LineWidth=1.25);
legendArrary={'Tx sequence','received sequence'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('Transmitted Sequence and without Prcode','n','',[0 200],[-1.5 1.5],legendArrary,flag,FontSize)


figure
stem(n(index),[signal(index).', pre_equalised(index)],LineWidth=1.25);
legendArrary={'Tx sequence','Precode transmitted sequence without modulo'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('Transmitted Sequence and without Modulo','n','',[0 200],[-4 4],legendArrary,flag,FontSize)

figure
stem(n(index),[signal(index).', pre_equalised_with_mod(index)],LineWidth=1.25);
legendArrary={'Tx sequence','Precode transmitted sequence with modulo'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('Transmitted Sequence and with Modulo','n','',[0 200],[-1.5 1.5],legendArrary,flag,FontSize)

figure
stem(n(index),[signal(index).',output(index)],LineWidth=1.25);
legendArrary={'Tx sequence','thp code received sequence without modulo'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('Received Sequence without Modulo','n','',[0 200],[-3 3],legendArrary,flag,FontSize)


figure
stem(n(index),[signal(index).',output_with_mod(index)],LineWidth=1.25);
legendArrary={'Tx sequence','thp code received sequence with modulo'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('Received Sequence with Modulo','n','',[0 200],[-1.5 1.5],legendArrary,flag,FontSize)
