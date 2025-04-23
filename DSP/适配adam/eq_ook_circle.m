clear;
close all;
addpath('MLSE\');
%信号和通道参
% 
% 数
%系统仿真参数
Fs = 1;                         %采样频率
Rb = 1;
N = 10000;                      %迭代次数
%调制信号参数
Rs = Fs;                        %符号速率
nSamp = Fs/Rs;                  %每符号采样

%信道参数
chnl = [0.2; 0.4; 0.7; 0.4;0.2]';       %信道冲激响应 可修改

%产生OOK信号
rng(1);  %随机种子
txSig = randi([0,1],1,N)*2-1;  %产生随机序列

%信号通过信道
filtSig = filter(chnl,1,txSig);

%给信号加awgn噪声
SNR = 30;
noisySig = awgn(filtSig,SNR,'measured');
% eyediagram(noisySig(end-1000:end),3)  %看眼图

%同步
[noisySig,txSig] = sync(noisySig,txSig);

%FFE
t_fe = 5;   %前向抽头系数 自己设定   
% t_fb = 3;
% ntrain = 10000;
m=0.84:0.001:0.999;
tic;f = waitbar(0,'Please wait...');
for i=1:length(m)
% m = 0.1;  %步长 一般为0.02 1e-3 2e-3 可自行修改
[y_ffe,e_ffe(i,:)] = adam(noisySig(2:end),txSig,t_fe,m(i));
% figure;plot(abs(e_ffe)); title('均方误差'); %MSE的值
% eyediagram(y_ffe(end-1000:end),3);
%判决

th = 0;
y(y_ffe>=th) = 1;  
y(y_ffe<=th) = -1;
%同步
[y1,txSig] = sync(y,txSig);
%计算误码率
[ber(i),k] = CalcBER(y1,txSig);
waitbar(i/length(m),f,'Please wait...');
end
waitbar(1,f,'Finishing');
pause(1)

close(f)
time = toc;
figure;plot(m,ber);xlabel ('\mu');ylabel ('BER');grid on;
e=abs(e_ffe);
% save('5taps 10SNR LMS','m','ber','time');

