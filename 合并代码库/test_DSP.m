clear;
close all;
% addpath(genpath('Fcns'));
addpath('algorithm存储\');
% datapath = 'Data\20230619_ook_PRBS15_7G_500mV_10km';
% prefix of file
% prefix = 'Att';
% load(sprintf('%s\\att_vec.mat',datapath));
% load(sprintf('%s\\pd_inpower.mat',datapath));

fs = 20e9; % sampling rate
sps = 2; % sample per symbol in dsp
osr = 4; % oversampling rate (compared with baud rate)
ncut = 2e3;
prbsOrder = 15;
M = 4;

fb = 5e9; % baud rate
baud = 4;
current = 300;
pathlength = '10km';
tap = 10;
delay = 2;

ber = [];
num = [];
%导入采集数据
rop = 10;
    % generate file name
    filename = sprintf('test_ook_5g_prbs15');%sprintf('%s\\%s-%1.1f.mat',datapath,prefix,att_vec(iAtt+8));
    % load
    load(filename);

    % load the reference
    load('awg_prbs15_reference.mat');
    ref2 = repelem(ref,sps);
    ref_seq = repmat(ref,1,1000);
    ref_seq = ref_seq(:);

    %数字时钟恢复
    data = resample(data,osr*fb,fs);  %四倍上采样，为了满足后续时钟恢复频率
    f_clk_rec = cr(data.',fb,osr,0); % 注意：f_clk_rec是对fb的估计，不是对fb*osr的估计！
    t_cur = 1/osr/fb*(0:length(data)-1); % 注意：插值前data1的采样率是osr*fb
    t_new = 0:1/f_clk_rec/osr:t_cur(end); % 注意：插值后data2的采样率是osr*f_clk_rec
    data1 = interp1(t_cur,data,t_new);
    data1 = resample(data1,sps,osr);
    data2 = std(ref)*normalization(data1.');  %归一化
    %     eyediagram(data2(end-1000+1:end),4)

    [data3,ref2] = sync(data2,ref2);

    %均衡
    eqMethod = 'QKLMS';
    switch eqMethod
        case 'LMS'
            [y0,e0,~,t] = ...
                LMS(data3(1:12000),ref_seq(1:end),0.001,tap,delay,sps);
            delay1 = delay;
        case 'LMK'
            [y0,e0,~,t] = ...
                LMK(data3(1:60000),ref_seq(1:end),0.0001,tap,delay,sps);
            delay1 = delay;
        case 'KLMS'
            [y0,e0,t] = ...
                KLMS(data3(1:12000),ref_seq(1:end),0.25,tap,delay,sps);
            delay1 = delay+1;
        case 'QKLMS_fixed_budget'
            db=1000;
            [y0,e0,t] = ...
                QKLMS_fixed_budget(data3(1:12000),ref_seq(1:end),0.25,tap,delay,sps,db);
            delay1 = delay+1;
        case 'DFE-KLMS'
            [y0,e0,t] = ...
                DFE_KLMS(data3(1:30000),ref_seq(1:end),[0.2 0.005],[tap,10],delay,sps);
            delay1 = delay+1;
        case 'KLMK'
            [y0,e0,t] = ...
                KLMK(data3(1:60000),ref_seq(1:end),0.07,tap,delay,sps);
            delay1 = delay+1;
        case 'QKLMS'
            [y0,e0,t] = ...
                QKLMS(data3(1:12000),ref_seq(1:end),0.2,tap,delay,sps);
            delay1 = delay+1;
        case 'QKLMK'
            [y0,e0,t] = ...
                QKLMK(data3(1:50000),ref_seq(1:end),0.01,tap,delay,sps);
            delay1 = delay+1;
        case 'VNLE3'
            [y0,e0,~,t] = ...
                VNLE3_RLS(data3(1:12000),ref_seq(1:end),tap,delay,sps);
            delay1 = delay;
        case 'RGKLMK'
            [y0,e0,~,t] = ...
                RGKLMK(data3(1:60000),ref_seq(1:end),0.002,tap,delay,sps);
            delay1 = delay;
        case ['NQK' ...
                'LMS']
            [y0,e0,t] = ...
                NQKLMS(data3(1:60000),ref_seq(1:end),0.5,tap,delay,sps);
            delay1 = delay+1;
        case 'MKLMS'
            [y0,e0,t] = ...
                MKLMS(data3(1:70000),ref_seq(1:end),0.5,tap,delay,sps);
            delay1 = delay+1;
        case 'MKLMS_NS'
            [y0,e0,t] = ...
                MKLMS_NS(data3(1:70000),ref_seq(1:end),0.5,tap,delay,sps);
            delay1 = delay+1;
    end
    e1 = e0(1:end).^2;
    e = smooth(e1,91);
    figure;semilogy(e(1:end));hold on;xlabel('Iteration');ylabel("MSE");
    ylim([0.005,10]);grid on;%xlim([1 end]);
    MSE = mean(e1(ncut:end));
    p1=find(e<=MSE);
    p=p1(1);
%     mean_e(1:length(e0))=MSE;
%     semilogy(mean_e)
    %判决
    [index,yyy] = quantiz(y0,[-2,0,2],[-3,-1,1,3]);
    %     [ref_seq,yyy] = sync(ref_seq,yyy);
    %计算误码率
    [ber,num] = CalcBER(ref_seq(ncut+delay1:end),yyy(ncut:end));

    %输出信息
    %     fprintf('%s: BER = %e, # of errors = %d\n',filename,ber(iAtt),num(iAtt));
    fprintf('%s; BER=%0.6f; NUM=%d; MES=%1.4f; P=%d; time=%2.4fs\n',eqMethod,ber,num,MSE,p,t);

%%

PowVec = pd_inpower(3:end-3);
ber1 = ber;
ber1(ber==0)=[];
ber(ber==0) = min(ber1);
% 
outputFilename = sprintf(['Outputs\\d20230308-CZ-PAM%d-PRBS%d-' ...
    '%dGbps-%dmV-%s-%s.mat'],M,prbsOrder,baud,current,pathlength,eqMethod);
if exist(outputFilename,'file')
    fprintf('同名文件存在，按任意键覆盖！\n');
    pause();
    save(outputFilename,'ber','PowVec');
else
    save(outputFilename,'ber','PowVec');
end

% sound(sin(2*pi*25*(1:4000)/100));
% eyediagram(data2(end-1000+1:end),4)
berplot = BERPlot();
h=berplot.plot(fliplr(-PowVec),fliplr(ber),1,1);
hold on