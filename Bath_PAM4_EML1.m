clc;clear;close all;
current_date = date;
disp(current_date);

addpath('D:\PhD\Project\Base_Code\Base\')
addpath("Data\2020.02.18\EML")
addpath("Data\2020.02.18\EML\28G")
addpath("Fncs\")
addpath("DSP\")
% 读取参考序列
load("refPAM4.mat")

fs = 80e9; % sampling rate
sps = 2; % sample per symbol in dsp
osr = 4; % oversampling rate (compared with baud rate)
fb = 28e9; % 信号的波特率 baud rate

% 参考序列
ref2 = repelem(ref,sps);
ref_seq = repmat(ref,1,1000);
ref_seq = ref_seq(:);
label=ref_seq;

Amp_NUM=12:-1:1;
% Amp_NUM=2;
% 文件存储
datapath='Output\PAM4_28G_EML_BTB_ffe';

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

% 数据所在文件位置
foldername='0km';
BER=zeros(length(Amp_NUM),1);
power=zeros(length(Amp_NUM),1);
for index=1:length(Amp_NUM)
    BER_NUM=zeros(length(Amp_NUM),1);
    for fileID = 1:2
        filename = sprintf('%s\\12dB-att-283mA_-1.5-DAC_0km_-%1.1fdBm-PRBS215-ch4-0.50-28G-%d.mat',...
            foldername,Amp_NUM(index),fileID);
        load(filename);
        data = data';

        % 数字时钟恢复
        xn=SignalRec.Time_recovery(data);
        % 均衡
        PAM_Equ;
        % 解码
        [decodedData,ber] = SignalRec.NRZ_ExecuteDecoding(sigRx_E);

        BER_NUM(fileID)=ber;

        close all;
    end
    ber=mean(BER_NUM);
    % 数据存储
    name=strcat('BER_',num2str(Amp_NUM(index)),'dBm');

    ds.name=name;
    ds.data=ber;
    ds.saveToMat();

    % 绘图准备
    BER(index)=ber;
    power(index)=-Amp_NUM(index);

end

berplot = BERPlot_David();

berplot.interval=2;
berplot.flagThreshold=1;
berplot.flagRedraw=1;
berplot.flagAddLegend=0;
berplot.orderInterpolation=1;
berplot.plot(power,BER);
