clc;clear;close all;
%% 信道响应定义与频率响应分析
Fs = 1000;  % 采样频率 (Hz)
N = 1024;   % FFT点数

% 定义信道响应系数
b = [1, -1, 3, zeros(1,7), 0.1];  % H(z) = 1 - z⁻¹ + 3z⁻² + 0.1z⁻¹⁰
a = 1;  % 分母系数 (FIR系统)

% 计算频率响应
[H, freq] = freqz(b, a, N, 'whole', Fs);

% 转换为单边频谱
H_mag = abs(H(1:N/2+1));
H_phase = angle(H(1:N/2+1));
freq = freq(1:N/2+1);

% 计算冲激响应
n = 0:20;  % 时间序列
h = zeros(size(n));
h(1) = b(1);  % n=0
for k = 2:length(b)
    if n(k) <= max(n)
        h(k) = b(k);
    end
end

%% 绘制结果
figure('Position', [100, 100, 900, 700])

% 1. 幅度响应
subplot(3,1,1)
plot(freq, 20*log10(H_mag)), grid on
title('信道幅度响应 |H(f)| (dB)')
xlabel('频率 (Hz)')
ylabel('增益 (dB)')
xlim([0, Fs/2])
ylim([-10, 20])

% 2. 相位响应
subplot(3,1,2)
plot(freq, unwrap(H_phase)*180/pi), grid on
title('信道相位响应 \angle H(f)')
xlabel('频率 (Hz)')
ylabel('相位 (度)')
xlim([0, Fs/2])

% 3. 冲激响应
subplot(3,1,3)
stem(n, h, 'filled', 'LineWidth', 1.5)
title('信道冲激响应 h[n]')
xlabel('采样点 n')
ylabel('幅度')
xlim([-0.5, 20.5])
grid on

%% 关键频率点分析
fprintf('===== 关键频率点分析 =====\n');
fprintf('直流响应 (0 Hz): %.2f dB\n', 20*log10(abs(H(1))));
fprintf('Nyquist频率 (%d Hz): %.2f dB\n', Fs/2, 20*log10(abs(H(N/2+1))));

% 找到峰值频率
[~, idx] = max(H_mag);
fprintf('最大增益频率: %.1f Hz (增益: %.2f dB)\n', freq(idx), 20*log10(H_mag(idx)));

%% 系统特性分析
% 计算零极点
zeros = roots(b);
poles = roots(a);

figure('Position', [200, 200, 800, 600])
zplane(zeros, poles)
title('零极点分布图')
grid on

%% 系统类型判断
% 计算直流和Nyquist频率增益
dc_gain = abs(H(1));
nyq_gain = abs(H(N/2+1));

if dc_gain > nyq_gain
    fprintf('\n===== 系统特性: 低通滤波器 =====\n');
elseif dc_gain < nyq_gain
    fprintf('\n===== 系统特性: 高通滤波器 =====\n');
else
    fprintf('\n===== 系统特性: 带通/全通滤波器 =====\n');
end

% 计算带宽
bw = sum(H_mag > max(H_mag)/sqrt(2)) * (Fs/2)/(N/2);
fprintf('3-dB带宽: %.1f Hz\n', bw);
