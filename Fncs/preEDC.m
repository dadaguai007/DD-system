%% pre EDC
clear all;
close all;
RanSeed = 16;
BitRateDefault = 50e9;
FrmLen = BitRateDefault/1e9*800;
SampleRateDefault = 2 * BitRateDefault;
UpRate = SampleRateDefault/BitRateDefault;
AWGSamplingRate = 224e9;
RTOSamplingRate = 224e9;
TrainLen = 1024;
N_T = 128;
RollOff =0.1;

% EDC 参数
PRE_EDC = 1;  %% 1：GS  2：FIR
CSPR = 15;
DER = 22;  %% 对于GS来说需要22左右   对于FIR来说在1左右
iteration = 50;


FiberLen = 50e3;
BitPerSym=2;
DREmode=2;

%% bit generation
rand('seed', RanSeed);
TxBit = randi([0 1],1, FrmLen*BitPerSym);
switch BitPerSym
    case 2 
        %% PAM4
        M = 2^BitPerSym;
        TxSym_temp = pam4(TxBit);
        TxSym =pammod(TxSym_temp,M);
    case 2.5
        M=32;
        TxSym_temp = qammod(TxBit.',M,'InputType','Bit');
        TxSym_temp = TxSym_temp.';
        TxSym = [real(TxSym_temp)  imag(TxSym_temp)];
    case 3
        M = 2^BitPerSym;
        TxSym_temp = pam8(TxBit);
        TxSym =pammod(TxSym_temp,M);
end
TrainSeq = TxSym(1:1024);

%% upsampling and RC filtering
TxWfm = upsample(TxSym,UpRate);
RRCFilt = rcosfir(RollOff, N_T, UpRate, 1, 'sqrt');
TxWfm = conv(TxWfm, RRCFilt, 'same');
disp("PAPR_beforeEDC")
PAPR = PAPR_cal(TxWfm)
figure;
plot(linspace(-SampleRateDefault/2,SampleRateDefault/2,length(TxWfm)),10*log10(abs(fftshift(fft(TxWfm)))));
title('Txwfm'); 
TxWfm = TxWfm/sqrt(mean(abs(TxWfm).^2));


%% PRE-EDC
 C_speed = 3e8;
switch PRE_EDC
    case 0
      TxWfm_PreEDC = TxWfm;
    case 1
        %% Digital extinction ratio Emulation
        TxWfm = TxWfm+ 10^(DER/20);
        disp(10*log(max(abs(TxWfm))/min(abs(TxWfm))))
        Pt=TxWfm;
        s0_t=sqrt(Pt);
        st_temp = s0_t;
        for i=1:iteration
            st = abs(st_temp);
            rt = FiberDispersion( st, C_speed, FiberLen, SampleRateDefault  );
            rt_temp = s0_t.*exp(1j*angle(rt));
             BW = BitRateDefault*(1+RollOff)/2+0.6e9;
            rt_temp = Brick_wall_filter(rt_temp,BW,SampleRateDefault);
            st_temp =FiberDispComp(rt_temp,C_speed, FiberLen, SampleRateDefault );
        end
        % 取功率信号
        TxWfm_PreEDC = st_temp.*conj(st_temp);
        TxWfm_PreEDC = TxWfm_PreEDC- mean(TxWfm_PreEDC);
        figure;
        plot(linspace(-SampleRateDefault/2,SampleRateDefault/2,length(TxWfm_PreEDC)),10*log10(abs(fftshift(fft(TxWfm_PreEDC)))));
        title('Txwfm');
    case 2
        % 使用FIR滤波器其进行计算
        Pre_filter=FIR_cal(FiberLen,SampleRateDefault,DER);
        TxWfm_PreEDC = ifft(fft(TxWfm).*fft(Pre_filter,length(TxWfm)));
        TxWfm_PreEDC = conv(TxWfm,Pre_filter,"same");
        TxWfm_PreEDC = TxWfm_PreEDC- mean(TxWfm_PreEDC);
        figure;
        plot(linspace(-SampleRateDefault/2,SampleRateDefault/2,length(TxWfm_PreEDC)),10*log10(abs(fftshift(fft(TxWfm_PreEDC)))));
        title('Txwfm');
end
disp("PAPR_afterEDC")
PAPR = PAPR_cal(TxWfm_PreEDC)
% C_speed = 3e8;
% FiberLen = 0*80e3;
% SampleRateDefault = BaudRate*Up;
lamda = C_speed/193.1e12;
CD_value = 17e-6 * FiberLen;
N = 30000;
TimeW = N/SampleRateDefault;
beta2 = -lamda^2*CD_value/2/pi/C_speed;
w = 2*pi*[(0:N/2-1),(-N/2:-1)]/TimeW;
% 理想的信道响应
Response_fft=cos(-beta2*(w.^2)/2);
% 转回时域
Response_Tdomain=ifft(Response_fft);
figure;
plot(linspace(-SampleRateDefault/2,SampleRateDefault/2,length(Response_Tdomain)),10*log10(abs(fftshift(fft(Response_Tdomain)))));
title('Response_Tdomain');