% 以此为模版，更改所有的算法数据

clc;clear;close all;
current_date = date;
disp(current_date);

addpath('D:\PhD\Project\Base_Code\Base\')
addpath("D:\PhD\IM_DD_DSP\Data\DML-NRZ-50G\")
addpath("D:\PhD\IM_DD_DSP\Data\DML-NRZ-50G\252mA")
addpath("Fncs\")
addpath("DSP\")
addpath('Plot\')

% 读取参考序列
load ref32768 % can be replaced by prbs.m
Ref = ref32768;

fs = 80e9; % sampling rate
sps = 2; % sample per symbol in dsp
osr = 4; % oversampling rate (compared with baud rate)
fb = 50e9; % 信号的波特率 baud rate

% 参考序列
ref_seq = repmat(Ref,1000,1);
ref_seq=ref_seq(:);
label=ref_seq;

% 文件存储路径
datapath='Output\DML_NRZ_50G_BTB_ffe';

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
    Ref);      % 参考信号

% 文件存储
% Vpp = [4.0,5.1,6.1,7.1,8.1,9.1,10.1,11.1,12.0,13.0,14.0,15.0];
Vpp=4:1:15;
Vpp=flip(Vpp);
BER=zeros(length(Vpp),1);
V_PP=zeros(length(Vpp),1);

% 数据所在文件位置
foldername = '252mA';

for powID = 1:length(Vpp)
    BER_NUM=zeros(length(Vpp),1);
    for fileID = 1:1
        filename = sprintf('%s\\12dB-at252mA@NRZ_0km_-%1.1fdBm-PRBS215-ch4-50G-%d.mat',...
            foldername,Vpp(powID),fileID);

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
    name=strcat('BER_',num2str(Vpp(powID)),'Vpp');

    ds.name=name;
    ds.data=ber;
    ds.saveToMat();
    
    % 绘图准备
    BER(powID)=ber;
    V_PP(powID)=-Vpp(powID);
end

% BER 绘图
berplot = BERPlot_David();

berplot.interval=2;
berplot.flagThreshold=1;
berplot.flagRedraw=1;
berplot.flagAddLegend=0;
berplot.orderInterpolation=1;
berplot.plot(V_PP,BER);