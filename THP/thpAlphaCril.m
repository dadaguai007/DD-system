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

% alpha=[0.01,0.2,0.5,0.8,0.99];
alpha=0.01:0.0245:0.99;
for i=1:length(alpha)
    Tx.Nr.psfRollOff=alpha(i);
    % 时间轴和频率轴创建
    [~,time] = freq_time_set(Tx.TxPHY.NSym *Tx.TxPHY.sps,fs);
    % sps
    sps=Tx.TxPHY.sps;

    % 信号生成
    [signal,pamsignal,pulse]=Tx.dataOutput();
    % DSP参考信号
    label=real(pamsignal);
    % PAM Signal
    pamsignal=pnorm(pamsignal);

    % 参考星座图
    [const,Ksym]=Tx.creatReferenceConstellation();

    % 数据性能起始位置
    skip = 0.1 * Tx.TxPHY.NSym; % 跳过初始瞬态（20%符号）
    Tx.Nr.ncut_index=skip;
    refEq=Tx.createReferenceSignal(label);

    % 绘图时间轴
    index=skip:skip+100;
    %% 信道创建
    syms z;                     % 声明符号变量z（用于Z域分析）
    % 第一种选择
    % c1=1;
    % c2=-0.4;
    % c3=0.5;
    % c4=0;

    c1=1;
    c2=0.2;
    c3=0.1;
    c4=0.05;

    channel = c1+c2*z^{-1}+c3*z^{-2};       % 定义信道传递函数：H(z) = 1 - z⁻¹（典型ISI信道模型）
    % channel = 1-z^{-1};
    % channel=1-z^{-1}+3*z^{-2}+0.1*z^{-10};
    % 创建类
    THP=THPClase(1, 1);

    % 定义信道响应系数
    b = [c1,c2,c3,c4];  % 1,0.5,0.7
    % b = [1,-2];% H(z) = 1 + z⁻¹
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

    % 时域冲击响应
    output_without_precoder=filter(b,a,signal);

    % 匹配滤波
    % 滤波器赋值
    Tx.TxPHY.hsqrt=pulse;
    matchOut=Tx.matchFiltering(output_without_precoder);
    % Eq
    EQ=struct();
    EQ.u=0.001;
    EQ.k1=21;
    EQ.k2=11;
    EQ.ref=11;
    EQ.sps=2;
    EQ.lamda=0.9999;
    EQ.delta=0.01;

    % FFE_LMS
    % [yEq,en,w] = FFE_LMS(EQ, matchOut.', refEq.');
    % DFE_LMS
    [yEq,en,w] = DFE_LMS(EQ,  matchOut.', refEq.');
    %plotEquParam(matchOut,yEq,refEq,w,en)
    % decode
    Tx.Button.Train='on';
    [decodedData,ber] =Tx. PAM_QuantifyDecoding(yEq);

    % THP
    tapsTHP=w(EQ.k1+1:end);

    % mod
    modN=4;
    % Pre THP
    [pre_equalised_with_mod,xInput]=THP.preTHP((pamsignal),tapsTHP,modN);
    % Pulse
    thpSignal=Tx.applyShapingFilter((pre_equalised_with_mod),pulse);

    % prs thp
    symbolsUp = upsample(pamsignal, sps);
    rolloff=alpha(i);
    psfLength=10;
    reverse='False';
    taps=3;
    pulsePrs=getPrsPulse(sps, psfLength, reverse, rolloff,taps);
    prsSignal=conv(symbolsUp,pulsePrs,'same');

    % pre thp
    thpsUp = upsample(pre_equalised_with_mod, sps);
    thpPrsSignal=conv(thpsUp,pulsePrs,'same');

    %%
    sz_window=1000;
    [~,ccdfy1]=ccdf(signal,sz_window,0);
    [ccdfx,ccdfy2]=ccdf(thpSignal,sz_window,0);
    %     [~,ccdfy3]=ccdf(thpPrsSignal,sz_window,0);
    %     [ccdfx,ccdfy4]=ccdf(prsSignal,sz_window,0);

    CCDF_woTHP(i,:)=ccdfy1;
    CCDF_THP(i,:)=ccdfy2;


    PAPR_dB_woTHP(i)= calc_papr(signal);
    PAPR_dB_wTHP(i)= calc_papr(thpSignal);


end

%%
Button_plot='plotpapr';
if strcmp(Button_plot,'plotpapr')
    color=distinguishable_colors(20);
    marker = 'so^d>';
    figure;hold on;
    plot(alpha,PAPR_dB_woTHP,LineWidth=1.25,Color=color(1,:),Marker=marker(1),MarkerFaceColor = 'w')
    plot(alpha,PAPR_dB_wTHP,LineWidth=1.25,Color=color(2,:),Marker=marker(2),MarkerFaceColor = 'w')
    legendArrary={'w/o THP','w THP'};
    FontSize=12;
    flag.LegendON_OFF=1;
    flag.Legendflage=0;
    Plotter('','alpha','PAPR(dB)',[0 1],[5 10],legendArrary,flag,FontSize)
end

%%
if strcmp(Button_plot,'plotccdm')
    color=distinguishable_colors(20);
    marker = 'so^d>';
    CCDFwoTHP=[CCDF_woTHP(1,:);CCDF_woTHP(2,:);CCDF_woTHP(3,:);CCDF_woTHP(4,:);CCDF_woTHP(5,:)];
    CCDFTHP=[CCDF_THP(1,:);CCDF_THP(2,:);CCDF_THP(3,:);CCDF_THP(4,:);CCDF_THP(5,:)];
    figure;
    for j=1:5
        semilogy(ccdfx,CCDFwoTHP(j,:),LineWidth=1.25,Color=color(j,:),Marker=marker(j),MarkerFaceColor = 'w');
        hold on;
    end
    for j=1:5
        semilogy(ccdfx,CCDFTHP(j,:),LineWidth=1.25,Color=color(j,:),Marker=marker(j),MarkerFaceColor=color(j,:));
    end
    legendArrary={'w/o THP alpha=0.01','w/o THP alpha=0.2','w/o THP alpha=0.5','w/o THP alpha=0.8',...
        'w/o THP alpha=0.99','w THP alpha=0.01','w THP alpha=0.2','w THP alpha=0.5','w THP alpha=0.8',...
        'w THP alpha=0.99'};
    FontSize=12;
    flag.LegendON_OFF=1;
    flag.Legendflage=0;
    Plotter('','PAPR_0 (dB)','Pr(PAPR > PAPR_0)',[0 10],[10^-3 10^0],legendArrary,flag,FontSize)
end