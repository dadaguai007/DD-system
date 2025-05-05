%% post-GS
Pt=Waveform;
% 幅度
s0_t=sqrt(Pt);
iteration = 100;
% 初始输入
rt_temp = s0_t;
C_speed = 3e8;
FiberLen = 58.45e3;
SampleRateDefault=RTOSamplingRate;
for i=1:iteration
    % 补偿色散
    rt = FiberDispComp( rt_temp, C_speed, FiberLen, SampleRateDefault  );
    % 取幅度
    rt_amp=  abs(rt);
    % 滤波
    BW = BitRateDefault*(1+RollOff)/2+0.1e9;
    rt_amp = Brick_wall_filter(rt_amp,BW,SampleRateDefault);
    % 色散仿真
    st_temp =FiberDispersion(rt_amp,C_speed, FiberLen, SampleRateDefault );
    % 取接收端相位
    rt_temp = s0_t.*exp(1j*angle(st_temp));
    rt_temp = Brick_wall_filter(rt_temp,BW,SampleRateDefault);
end
if iteration==0
    Waveform = Pt;
else 
    Waveform = rt_temp;
    % 再补充一次色散补偿
    Waveform = FiberDispComp(Waveform, C_speed, FiberLen, SampleRateDefault  );
    Waveform = abs(Waveform).^2;
end