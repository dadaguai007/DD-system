clc;clear;close all;
current_date = date;
disp(current_date);

addpath('D:\PhD\Project\Base_Code\Base\')
addpath("Data\pam4 Data\")
addpath("Fncs\")
addpath("DSP\")
% 读取参考序列
load("AWG_PAM4_PRBS9_Reference_Sequence.mat")
addpath("Data\20220616_PRBS9_16G_10km2\")
fs = 20e9; % sampling rate
sps = 2; % sample per symbol in dsp
osr = 4; % oversampling rate (compared with baud rate)
fb = 16e9; % 信号的波特率 baud rate

% 参考序列
ref2 = repelem(ref,sps);
ref_seq = repmat(ref,1,1000);
ref_seq = ref_seq(:);
label=ref_seq;
PRE='Att-';
load('pd_inpower.mat')
Amp_NUM=pd_inpower;
load('att_vec.mat')

% 文件存储
datapath='Output\PAM4_16G_10km_ffe';

% 装载数据保存模块
ds=DataSaver([], datapath,[]);
ds.createFolder();

% 装载DSP-CDR 和解码模块
SignalRec  =  DspSyncDecoding(...
    fs,...    % 接收信号的采样率
    fb,...    % 接收信号的波特率
    2,...     % 接收信号的调制格式
    sps,...   % 上采样率
    osr,...   % 时钟信号的上采样率
    1e5,...   % 误码计算起始位置
    ref);      % 参考信号


BER=zeros(length(Amp_NUM),1);
power=zeros(length(Amp_NUM),1);
for index=1:length(Amp_NUM)
    close all;
    Title=strcat(PRE,num2str(att_vec(index)));
    dataname = strcat(Title,'.mat');
    fprintf('the signal power is %.2f dBm\n',Amp_NUM(index));
    % 读取数据
    load(dataname);

    % 数字时钟恢复
    xn=SignalRec.Time_recovery(data);
    % 均衡
    PAM_Equ;
    % 解码
    [decodedData,ber] = SignalRec.NRZ_ExecuteDecoding(sigRx_E);


    % 数据存储
    name=strcat('BER_',num2str(Amp_NUM(index)),'dBm');

    ds.name=name;
    ds.data=ber;
    ds.saveToMat();
    
    BER(index)=ber;
    power(index)=-Amp_NUM(index);

end


BER=flip(BER);
power=flip(power);


berplot = BERPlot_David();

berplot.interval=2;
berplot.flagThreshold=1;
berplot.flagRedraw=1;
berplot.flagAddLegend=0;
berplot.orderInterpolation=1;
berplot.plot(power,BER);
