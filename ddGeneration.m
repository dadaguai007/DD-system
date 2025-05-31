addpath('Fncs\')
addpath('Plot\')
addpath('Dsp\')


% 发射机参数
Tx=DDTx(     ...
    64e9, ...                % 发射信号的采样率
    32e9, ...                % 发射信号的波特率
    2, ...                   % 随机信号的阶数
    15, ...                  % prbs码的阶数
    4, ...                   % 调制格式 M
    64e9/32e9, ...           % 每符号采样点
    1e5, ...                 % 码元数目
    'sqrt', ...              % 脉冲形式
    0.2, ...                 % 滚降系数
    10, ...                  % 影响长度
    'nrz', ...               % 用户的成型滤波器
    'rand', ...              % 选择模式
    'system');               % 成型滤波器的生成方式

fs=Tx.TxPHY.fs;
fb=Tx.TxPHY.fb;
Ts=1/fb;
Ta= 1/fs;
M=Tx.TxPHY.M;

% 成型滤波器
hsqrt = Tx.systemHsqrt();


Fc      = 193.1e12; % central optical frequency
% mzm 调制器参数
Vpi = 5;
Vb = -Vpi/2;
% PD 参数
paramPD=struct();
paramPD.B =fb;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=fs;
% ssfm
param=struct();
param.Ltotal = 40; %km
param.Lspan =10;
param.hz= 0.5;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = Fc;
param.NF = 4.5;
param.amp='edfa';
param.Fs=fs;

% 信号功率设置模式
channelPowerType='output';
% 信号功率
Pout_dBm=0;
