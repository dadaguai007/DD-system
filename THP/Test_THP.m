clc;clear;close all;
% addpath('D:\BIT_PhD\DD-system\Fncs\');
addpath('D:\PhD\DD-system\Fncs\');
addpath('D:\PhD\DD-system\Plot\')
% 系统参数设置
N = 6;                    % 模运算的模数
x = ones(8,1)*0.75;         % 生成8点原始信号序列，每个元素为0.75
syms z;                     % 声明符号变量z（用于Z域分析）
channel = 1 - z^(-1);       % 定义信道传递函数：H(z) = 1 - z⁻¹（典型ISI信道模型）
% channel=1-z^{-1}+3*z^{-2}+0.1*z^{-10};
% 创建类
THP=THPClase(1, 1);


% 信号参数
sps = 2;
Rs  = 40e9;
Ts  = 1/Rs ;
Fs  = sps*Rs;
Ta  = 1/Fs;

%PAM
M=4;
data_2bit=randi([0,1],log2(M),800);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;
% Mapeia bits para pulsos
symbTx = pammod(symbols,M,0,'gray');
% symbTx=real(symbTx);
symbTx = pnorm(symbTx);


% 将时域信号转换为Z域表达式
z_x=THP.Transz(symbTx);

% 信道逆预均衡
channel_inverse = 1/channel; % 计算信道逆：H⁻¹(z) = 1/(1 - z⁻¹)
z_pre_equalised = z_x * channel_inverse; % 预均衡：X'(z) = X(z)·H⁻¹(z)

% 预均衡后的信号
pre_equalised1=THP.zTrans(z_pre_equalised,symbTx);% 逆Z变换得到时域表达式

% 模运算控制信号幅度
pre_equalised_with_mod1 = THP.modulo(pre_equalised1, N); % 将信号限制在[-N/2, N/2)


% Upsampling
symbolsUp1 = upsample(pre_equalised1, sps);

% Puls
hsqrt = rcosdesign(0.5,20,sps,'sqrt');
% pulse shaping
sigTx1=conv(symbolsUp1,hsqrt,'same');

% Upsampling
symbolsUp2 = upsample(pre_equalised_with_mod1, sps);

% Puls
hsqrt = rcosdesign(0.5,20,sps,'sqrt');
% pulse shaping
sigTx2=conv(symbolsUp2,hsqrt,'same');



% Upsampling
symbolsUp = upsample(symbTx, sps);

% Puls
hsqrt = rcosdesign(0.5,20,sps,'sqrt');
% pulse shaping
sigTx=conv(symbolsUp,hsqrt,'same');

sz_window=8;
[~,ccdfy6]=ccdf(sigTx,sz_window,0);
[ccdfx,ccdfy7]=ccdf(sigTx2,sz_window,0);




color=distinguishable_colors(20);
marker = 'so^d>v*phx';
CCDF=[ccdfy6;ccdfy7];

for i=1:2
semilogy(ccdfx,CCDF(i,:),LineWidth=1.25,Color=color(i,:),Marker=marker(i));
hold on;
end
legendArrary={'Without THP','With THP'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('','PAPR_0 (dB)','Pr(PAPR > PAPR_0)',[0 8],[10^-3 10^0],legendArrary,flag,FontSize)
%% ========== 无预编码的传输仿真 ==========
% 将时域信号转换为Z域表达式
z_x=THP.Transz(sigTx);

% 信号通过信道传输
z_output = z_x * channel;   % 信道输出Y(z) = X(z)·H(z)
output_without_precoder=THP.zTrans(z_output,sigTx);% 逆Z变换得到时域表达式


n = [0:length(sigTx)-1]; % 时间轴

index=1:200;
% 图1：原始信号与无预编码接收信号对比
figure
stem(n(index), [sigTx(index).', output_without_precoder(index)],LineWidth=1.25);
% title('Transmitted Sequence and without Prcode');
% xlabel('n');
% legend({'original sequence','received sequence'}, 'Location','northwest');
legendArrary={'original sequence','received sequence'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('Transmitted Sequence and without Prcode','n','',[0 200],[-1.5 1.5],legendArrary,flag,FontSize)



sz_window=8;
[~,ccdfy1]=ccdf(output_without_precoder,sz_window,0);
%% ========== THP预编码处理 ==========
% 信道逆预均衡
channel_inverse = 1/channel; % 计算信道逆：H⁻¹(z) = 1/(1 - z⁻¹)
z_pre_equalised = z_x * channel_inverse; % 预均衡：X'(z) = X(z)·H⁻¹(z)

% 预均衡后的信号
pre_equalised=THP.zTrans(z_pre_equalised,sigTx);% 逆Z变换得到时域表达式

[~,ccdfy2]=ccdf(pre_equalised,sz_window,0);
% 模运算控制信号幅度
pre_equalised_with_mod = THP.modulo(pre_equalised, N); % 将信号限制在[-N/2, N/2)
[~,ccdfy3]=ccdf(pre_equalised_with_mod,sz_window,0);
% 将预编码信号重新转换为Z域
zz_pre_equalised=THP.Transz(pre_equalised_with_mod);

% 预编码信号通过实际信道
z_output = zz_pre_equalised * channel; % Y(z) = X'(z)·H(z)
output=THP.zTrans(z_output,sigTx); % 接收端时域信号
[~,ccdfy4]=ccdf(output,sz_window,0);
output_with_mod = THP.modulo(output, N);  % 接收端模运算恢复信号
[ccdfx,ccdfy5]=ccdf(output_with_mod,sz_window,0);



figure
stem(n(index),[sigTx(index).', pre_equalised(index)],LineWidth=1.25);
% title('Transmitted Sequence without Modulo'); 
% xlabel('n');
% legend({'original sequence','transmitted sequence without modulo'},'Location','northwest')
legendArrary={'original sequence','transmitted sequence without modulo'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('Transmitted Sequence and without Modulo','n','',[0 200],[-4 10],legendArrary,flag,FontSize)






figure
stem(n(index),[sigTx(index).', pre_equalised_with_mod(index)],LineWidth=1.25);
% title('Transmitted Sequence with Modulo'); 
% xlabel('n');
% legend({'original sequence','transmitted sequence with modulo'},'Location','southwest')
legendArrary={'original sequence','transmitted sequence with modulo'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('Transmitted Sequence and with Modulo','n','',[0 200],[-1.5 1.5],legendArrary,flag,FontSize)


figure
stem(n(index),[sigTx(index).',output(index)],LineWidth=1.25);
% title('Received Sequence without Modulo'); 
% xlabel('n');
% legend({'original sequence','received sequence without modulo'},'Location','southwest')
legendArrary={'original sequence','received sequence without modulo'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('Received Sequence without Modulo','n','',[0 200],[-3 3],legendArrary,flag,FontSize)


figure
stem(n(index),[sigTx(index).',output_with_mod(index)],LineWidth=1.25);
% title('Received Sequence with Modulo'); 
% xlabel('n');
% legend({'original sequence','received sequence with modulo'},'Location','southwest')
legendArrary={'original sequence','received sequence with modulo'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('Received Sequence with Modulo','n','',[0 200],[-1.5 1.5],legendArrary,flag,FontSize)



color=distinguishable_colors(20);
marker = 'so^d>v*phx';
% 无pre code
for i=1:1
semilogy(ccdfx,ccdfy1,LineWidth=1.25,Color=color(i,:),Marker=marker(i));
hold on;
end
legendArrary={'Without THP'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('','PAPR_0 (dB)','Pr(PAPR > PAPR_0)',[0 9],[10^-3 10^0],legendArrary,flag,FontSize)

%%
figure;
% 无mod
for i=1:1
semilogy(ccdfx,ccdfy2,LineWidth=1.25,Color=color(i,:),Marker=marker(i));
hold on;
end
legendArrary={'THP without Mod'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('','PAPR_0 (dB)','Pr(PAPR > PAPR_0)',[0 6],[10^-2 10^0],legendArrary,flag,FontSize)

figure;
% mod
for i=1:1
semilogy(ccdfx,ccdfy3,LineWidth=1.25,Color=color(i,:),Marker=marker(i));
hold on;
end
legendArrary={'THP with Mod'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('','PAPR_0 (dB)','Pr(PAPR > PAPR_0)',[0 8],[10^-3 10^0],legendArrary,flag,FontSize)


figure;
CCDM1=[ccdfy2;ccdfy3;];
% mod and with mod 
for i=1:2
semilogy(ccdfx,CCDM1(i,:),LineWidth=1.25,Color=color(i,:),Marker=marker(i));
hold on;
end
legendArrary={'Tx: THP without Mod','Tx: THP with Mod'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('','PAPR_0 (dB)','Pr(PAPR > PAPR_0)',[0 8],[10^-3 10^0],legendArrary,flag,FontSize)

%%

figure;
% 无mod
for i=1:1
semilogy(ccdfx,ccdfy4,LineWidth=1.25,Color=color(i,:),Marker=marker(i));
hold on;
end
legendArrary={'THP without Mod'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('','PAPR_0 (dB)','Pr(PAPR > PAPR_0)',[0 9],[10^-3 10^0],legendArrary,flag,FontSize)



figure;
% mod
for i=1:1
semilogy(ccdfx,ccdfy5,LineWidth=1.25,Color=color(i,:),Marker=marker(i));
hold on;
end
legendArrary={'THP with Mod'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('','PAPR_0 (dB)','Pr(PAPR > PAPR_0)',[0 7],[10^-3 10^0],legendArrary,flag,FontSize)



figure;
CCDM2=[ccdfy4;ccdfy5;];
% mod and with mod 
for i=1:2
semilogy(ccdfx,CCDM2(i,:),LineWidth=1.25,Color=color(i,:),Marker=marker(i));
hold on;
end
legendArrary={'Rx: THP without Mod','Rx: THP with Mod'};
FontSize=12;
flag.LegendON_OFF=1;
flag.Legendflage=0;
Plotter('','PAPR_0 (dB)','Pr(PAPR > PAPR_0)',[0 8],[10^-3 10^0],legendArrary,flag,FontSize)