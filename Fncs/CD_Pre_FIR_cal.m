function Pre_filter=CD_Pre_FIR_cal(FiberLen,SampleRateDefault,DER)
% 实现了一个预均衡滤波器设计系统，通过迭代优化生成一个预加重滤波器（Pre_filter），
% 用于在发射端预先补偿光纤传输中的色散（CD）和非线性效应。
a=1;
b= 10^(DER/10)*a; % 根据DER计算中间带增益
N=900;  % 滤波器半长度
H_FIR_Orin=sqrt([a*ones(1,N) b a*ones(1,N+1)]); % 初始频响：中间带增益√b，两侧√a
iteration = 40;
C_speed = 3e8;
H_FIR= H_FIR_Orin;
for i=1:iteration
    % % 步骤1：发送端滤波（仅幅度影响）
    H_FIR_Tx = abs(H_FIR);
     % 步骤2：模拟光纤色散效应
    rt = FiberDispersion( H_FIR_Tx, C_speed, FiberLen, SampleRateDefault  );
%     BW = BitRateDefault*(1+RollOff)/2+0.6e9;
%     rt = Brick_wall_filter(rt,BW,SampleRateDefault);
    % 步骤3：接收端相位响应提取
    H_FIR_Rx = abs(H_FIR_Orin).*exp(1j*angle(rt));
     % 步骤4：色散补偿计算
    H_FIR =FiberDispComp(H_FIR_Rx,C_speed, FiberLen, SampleRateDefault );
end
% % 去直流处理（取幅度，适用于DD系统）
Pre_filter = abs(H_FIR).^2 - mean(abs(H_FIR).^2);
figure;
plot(linspace(-SampleRateDefault/2,SampleRateDefault/2,length(Pre_filter)),10*log10(abs(fftshift(fft(Pre_filter)))));
title('Pre_filter');

end