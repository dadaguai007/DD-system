% whitening Noise
% 1 sps 进行相应操作
% 在均衡之后使用
% 信号路径：补偿信道频率响应偏差
% 噪声路径：实现噪声白化（使噪声频谱平坦化）

% 白化滤波器参数
switch log2(M)
    case 1
        Coe = zeros(1,8);
        step = 0.005;
    case 2
        Coe = zeros(1,4);
        step = 0.005;
    case 3
        Coe = zeros(1,4);
        step = 0.004;
    case 4
        Coe = zeros(1,2);
        step = 0.004;
end


% 提取发射数据
txSym=label(1:length(sigRx_E));
txSym=pnorm(txSym);

% 噪声提取
noise = sigRx_E - txSym;

% 获取白化滤波参数
[Shaping_Tap, error] = APR_PostFilter_LMS_Train(sigRx_E, label, step, Coe );
% [Shaping_Tap, error] = Modified_APR_PostFilter_LMS_Train( RxDown, label, step, Coe );
% [Shaping_Tap, error,RxDown] = PR_FFE_RLS_Train( RxDown, TrainSeq, delay, lambda, Up, Coe );
% lambda = 0.9999;
% [Shaping_Tap, error] = PR_PostFilter_RLS_Train( RxDown, TrainSeq, lambda, Coe );
% [Shaping_Tap, error] = Modified_PR_PostFilter_RLS_Train( RxDown, TrainSeq, lambda, Coe );

% 接收信号白化
sigNoiseWhiten =filter(Shaping_Tap,1,sigRx_E);
% 噪声白化
noise =filter(Shaping_Tap,1,noise);

% 创建频率轴
[fre,~]=freq_time_set(noise,fs);

% 噪声谱
mon_ESA(noise,fs);
title('noise');

% figure;
% plot(fre,length(noise)),10*log10(abs(fftshift(fft(noise)))));
% title('noise');

% 选取相应的参考星座图
constellation = unique(label);
% 噪声白化后，使用MLSE进行解码
sigRx=mlseeq(sigNoiseWhiten,Shaping_Tap,constellation,1000,'rst',1);



%% PNC
alpha = 0:0.01:1; % 参数值 可改变

sndsym=sigRx_E; % 接收信号 为滤波器后的输出 未进行下采样 为：2sps
sndsym = sndsym/sqrt(mean(abs(sndsym).^2))*sqrt(5);
pncsym=sigRx_E;
pncsym = pncsym/sqrt(mean(abs(pncsym).^2))*sqrt(5);


pec_out = pec(sndsym, .36, constellation);

% 三阶段的判决循环
pnc_out_1 = pnc(pncsym, .29, constellation, 1);

pnc_out_2 = pnc(pncsym, .11, constellation, 2);

pnc_out_3 = pnc(pncsym, .05, constellation, 3);


% Test

if 1
    % Encoding 经过编码
    % Linear encoding
    alpha = 0.5;

    % 对信号进行编码，可有可无
    % sig信号为初始的信号，后续要经过脉冲成型
    lcoeff = [1, alpha];
    % Flitering
    prs_sig = filter(lcoeff, 1, sig);

    % ffe 过后的数据，经过PNC或者PEC

    pnc_out = pnc(sigRx_E, alpha, constellation);

    % 增强PNC
    my_pnc_out = pnc_enhanced(sigRx_E, alpha, constellation);
    my_pnc_out = pnc_enhanced(my_pnc_out, alpha, constellation);

else
    % PRS_PNC

    % Linear encoding
    D = 0.8;
    lcoeff = prsFilter(D, 2, 1);
    lcoeff = lcoeff / sum(lcoeff);

    % Flitering
    prs_sig = filter(lcoeff, 1, sig);

    pnc_out = pnc(sigRx_E, alpha, constellation);

    % 对比白化滤波之后，经过MLSE的性能
    % out_ffe = mlseeq(ffe_out, lcoeff, constellation, 1000, 'rst');

end
