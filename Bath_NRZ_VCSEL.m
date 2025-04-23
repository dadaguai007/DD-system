clc;clear;close all;
current_date = date;
disp(current_date);

addpath('D:\PhD\Project\Base_Code\Base\')
addpath("D:\PhD\IM_DD_DSP\Data\DataSet_28Gbps_VCSEL\")
addpath("D:\PhD\IM_DD_DSP\Data\DataSet_28Gbps_VCSEL\Data\20161223_VCSEL\")
addpath("Fncs\")
addpath("DSP\")

% 读取参考序列
load PRBS15 % can be replaced by prbs.m
ref = +bitSeq;

fs = 80e9; % sampling rate
sps = 2; % sample per symbol in dsp
osr = 4; % oversampling rate (compared with baud rate)
fb = 28e9; % 信号的波特率 baud rate

% 参考序列
ref2 = repelem(ref,sps);
ref_seq = repmat(ref,1,1000);
ref_seq = ref_seq(:);
label=2*ref_seq-1;


% 文件存储
datapath='Output\NRZ_28G_BTB_ffe';
if ~exist(datapath,'dir')
    mkdir(datapath);
end

% 文件存储
Vpp = 0.5:0.1:1.1;
BER=zeros(length(Vpp),1);
V_PP=zeros(length(Vpp),1);

foldername = '20161223_VCSEL';


for powID = 1:length(Vpp)
    BER_NUM=zeros(length(Vpp),1);
    for fileID = 1:5
        filename = sprintf('Data\\%s\\OBTB_PRBS15_28G_RP-9.3_Vpp_%1.1f-%d.mat',...
            foldername,Vpp(powID),fileID);
        load(filename);
        data = data';
        %数字时钟恢复
        Time_recovery;
        % 均衡
        PAM_Equ;
        % 解码
        NRZ_Decode;
        BER_NUM(fileID)=ber;
        close all;
    end
    
    ber=mean(BER_NUM);
    % 数据存储
    name=strcat('BER_',num2str(Vpp(powID)),'Vpp');
    %     save(sprintf('%s\\%s.mat',datapath,name),'ber');
    BER(powID)=ber;
    V_PP(powID)=Vpp(powID);
end


berplot = BERPlot_David();

berplot.interval=2;
berplot.flagThreshold=1;
berplot.flagRedraw=1;
berplot.flagAddLegend=0;
berplot.orderInterpolation=1;
berplot.plot(V_PP,BER);