clc;clear;close all;
current_date = date;
disp(current_date);

addpath('D:\PhD\Project\Base_Code\Base\')
addpath("Data\SOA1\")
addpath("Fncs\")
addpath("DSP\")
% 读取参考序列
load("ref_d2.mat")
fs = 80e9; % sampling rate
sps = 2; % sample per symbol in dsp
osr = 4; % oversampling rate (compared with baud rate)
fb = 28e9; % 信号的波特率 baud rate

% 参考序列
ref=xx;
ref2 = repelem(ref,sps);
ref_seq = repmat(ref,1,100);
ref_seq = ref_seq(:);
label=ref_seq;

PRE='20240526_pam4_28G_150mA_SOAin_';
Amp_NUM=-8;
pre='ROP-';


% 文件存储
Datapath='Output\PAM4_SOA_28G_BTB_ffe';
if ~exist(Datapath,'dir')
    mkdir(Datapath);
end



for index=1:length(Amp_NUM)
    close all;
    Title=strcat(PRE,num2str(Amp_NUM(index)));
    datapath = strcat(Title,'_btb');
    fprintf('the SOA_in is %d dBm\n',Amp_NUM(index));
    addpath(strcat('D:\PhD\IM_DD_DSP\Data\SOA1\',datapath))
    load('pd_inpower.mat')
    BER=zeros(length(pd_inpower),1);
    for Index=1:length(pd_inpower)
        close all;
%         pd_inpower=flip(pd_inpower);
        power = sprintf('%.1f.mat', pd_inpower(Index));
        title=strcat(pre,power);
        fprintf('the signal power is %.3f dBm\n',pd_inpower(Index));
        % 读取数据
        load(title);

        %数字时钟恢复
        Time_recovery;
        % 均衡
        PAM_Equ;
        % 解码
        PAM_Decode;

        % 数据存储
        name=strcat('BER_',num2str(floor(pd_inpower(Index))),'dBm');
%         save(sprintf('%s\\%s.mat',Datapath,name),'ber');
        BER(Index)=ber;


    end
end
pd_inpower=flip(pd_inpower);
BER=flip(BER);
berplot = BERPlot_David();

berplot.interval=1;
berplot.flagThreshold=1;
berplot.flagRedraw=1;
berplot.flagAddLegend=0;
berplot.orderInterpolation=1;
berplot.plot(pd_inpower,BER);
