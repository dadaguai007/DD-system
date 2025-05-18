clc;clear;close all;
addpath('Fncs\');
addpath('THP\');
% 系统参数设置
N = 3;                    % 模运算的模数
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

% Upsampling
symbolsUp = upsample(symbTx, sps);

% Puls
hsqrt = rcosdesign(0.5,20,sps,'sqrt');
% pulse shaping
sigTx=conv(symbolsUp,hsqrt,'same');


%% ========== 无预编码的传输仿真 ==========
% 将时域信号转换为Z域表达式
z_x=THP.Transz(sigTx);

% 信号通过信道传输
z_output = z_x * channel;   % 信道输出Y(z) = X(z)·H(z)
output_without_precoder=THP.zTrans(z_output,sigTx);% 逆Z变换得到时域表达式


n = [0:length(sigTx)-1]; % 时间轴

% 图1：原始信号与无预编码接收信号对比
figure
stem(n, [sigTx.', output_without_precoder]);
title('发射与接收序列对比（无预编码）');
xlabel('时间索引 n');
legend({'发射序列','接收序列'}, 'Location','northwest');


%% ========== THP预编码处理 ==========
% 信道逆预均衡
channel_inverse = 1/channel; % 计算信道逆：H⁻¹(z) = 1/(1 - z⁻¹)
z_pre_equalised = z_x * channel_inverse; % 预均衡：X'(z) = X(z)·H⁻¹(z)

% 预均衡后的信号
pre_equalised=THP.zTrans(z_pre_equalised,sigTx);% 逆Z变换得到时域表达式

% 模运算控制信号幅度
pre_equalised_with_mod = THP.modulo(pre_equalised, N); % 将信号限制在[-N/2, N/2)

% 将预编码信号重新转换为Z域
zz_pre_equalised=THP.Transz(pre_equalised_with_mod);

% 预编码信号通过实际信道
z_output = zz_pre_equalised * channel; % Y(z) = X'(z)·H(z)
output=THP.zTrans(z_output,sigTx); % 接收端时域信号

output_with_mod = THP.modulo(output, N);  % 接收端模运算恢复信号




figure
stem(n,[sigTx.', pre_equalised]);
title('Transmitted Sequence without Modulo'); 
xlabel('n');
legend({'original sequence','transmitted sequence without modulo'},'Location','northwest')

figure
stem(n,[sigTx.', pre_equalised_with_mod]);
title('Transmitted Sequence with Modulo'); 
xlabel('n');
legend({'original sequence','transmitted sequence with modulo'},'Location','southwest')


figure
stem(n,[sigTx.',output]);
title('Received Sequence without Modulo'); 
xlabel('n');
legend({'original sequence','received sequence without modulo'},'Location','southwest')


figure
stem(n,[sigTx.',output_with_mod]);
title('Received Sequence with Modulo'); 
xlabel('n');
legend({'original sequence','received sequence with modulo'},'Location','southwest')

