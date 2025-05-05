%% Noise whitening
TxSym =pammod(TxSym_temp,M); % 为最普通的PAM序列
%         RxDown = RxDown/sqrt(mean(abs(RxDown).^2))*sqrt(5);
switch BitPerSym
    case 1
        Coe = zeros(1,8);
        step = 0.005;
    case 1.5
        Coe = zeros(1,4);
        step = 0.005;
    case 2
        Coe = zeros(1,4);
        step = 0.005;
    case 2.5
        Coe = zeros(1,4);
        step = 0.004;
    case 3
        Coe = zeros(1,4);
        step = 0.004;
    case 3.5
        Coe = zeros(1,3);
        step = 0.004;
    case 4
        Coe = zeros(1,2);
        step = 0.004;
end
Up = 1;
% [Shaping_Tap, error,RxDown] = PR_FFE_RLS_Train( RxDown, TrainSeq, delay, lambda, Up, Coe );
lambda = 0.9999;
% [Shaping_Tap, error] = PR_PostFilter_RLS_Train( RxDown, TrainSeq, lambda, Coe );
% [Shaping_Tap, error] = Modified_PR_PostFilter_RLS_Train( RxDown, TrainSeq, lambda, Coe );

%% 抽头越多，step越小 0.001-0.0001
noise = RxDown - TxSym;
[Shaping_Tap, error] = APR_PostFilter_LMS_Train(RxDown, TrainSeq, step, Coe );
% [Shaping_Tap, error] = Modified_APR_PostFilter_LMS_Train( RxDown, TrainSeq, step, Coe );
Shaping_Tap
%             Shaping_Tap =[1 0.0];
%              Shaping_Tap =[ 1.0000    0.1371   -0.3429    0.1477   -0.1285    0.0864];

% 信号路径：补偿信道频率响应偏差
% 噪声路径：实现噪声白化（使噪声频谱平坦化）

% 接收信号白化
RxDown =filter(Shaping_Tap,1,RxDown);
% 噪声白化
noise =filter(Shaping_Tap,1,noise);
figure;
plot(linspace(-BitRateDefault/2,BitRateDefault/2,length(noise)),10*log10(abs(fftshift(fft(noise)))));
title('noise');
constellation = unique(TrainSeq);
% 噪声白化后，使用MLSE进行解码
RxDown =mlseeq(RxDown,Shaping_Tap,constellation,1000,'rst',1);
