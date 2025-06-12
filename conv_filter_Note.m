% conv 与 filter 函数如此使用，效果等价
channel_result = filter(channel, 1, trSymVec);

channel_result1 = conv(trSymVec, channel);
channel_result1 = channel_result1(1:end-channel_length+1);

%% 信道对信号的作用，及作用函数
clc;clear;close all;
% 定义输入信号和信道响应
x = [2, -1, 3, 4, 0, 7];        % 任意输入信号
h = [1, -1, 3, 0,0,0,0,0,0,0,0.1]; % 您的信道冲激响应

% 方法1: 直接卷积 (时域)
y_conv = conv(h, x);

% 方法2: 使用filter函数
y_filter = filter(h, 1, [x, zeros(1, length(h)-1)]); % 补零避免截断

% 方法3: 频域相乘 (Z域操作)
N = length(x) + length(h) - 1;   % 输出长度
X_z = fft(x, N);                 % 通过FFT计算Z变换采样
H_z = fft(h, N);
Y_z = H_z .* X_z;                % 频域相乘
y_freq = ifft(Y_z, 'symmetric'); % 逆变换回时域

% 比较结果
fprintf('卷积与filter误差: %.2e\n', max(abs(y_conv - y_filter)));
fprintf('卷积与频域误差: %.2e\n', max(abs(y_conv - y_freq)));
