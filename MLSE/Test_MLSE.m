%% ---------------------- 系统初始化与参数设置 ----------------------
clc; clear; close all;   % 清空工作区、命令窗口，关闭所有图形窗口
transfer_function = [0.6, -1, 0.8] / sqrt(2);  % 信道冲激响应（归一化）
No_symbols = 1e5;        % 仿真符号数量
mu = 2;                  % 回溯步长基数
symbols = [-3 -1 1 3];   % 4-PAM符号集 % 也可以写成OOK模式[-1 1]
initialstate = [-1 -1];  % 编码器初始/终止状态
M = length(symbols);     % 调制阶数（4）
Eb = mean(symbols.^2);   % 符号平均能量计算
SNR_db = 0:2:16;         % 信噪比范围（dB）
N_array = mu * [1, 2, 4, 5, 10]; % 回溯深度参数
%% ---------------------- 生成发射符号序列 ----------------------
% 生成随机符号序列，并添加首尾初始状态
tx_symbols = symbols(randi([1, M], 1, No_symbols)); % 随机生成1e5个符号
tx_symbols = [initialstate tx_symbols initialstate]; % 添加首尾状态约束

% MLSE初始化
MLSE=mlseEqu(length(transfer_function), N_array(2),M,symbols,mu);

% 状态转置图
MLSE.getStateTrellis(transfer_function);

MLSE.getRectStateTrellis(transfer_function);

%% ---------------------- 信道传输模拟 ----------------------
% 计算信道输出（含首尾状态）
channel_output = zeros(1, length(tx_symbols)-2);
for i = 1:length(tx_symbols)-2
    channel_output(i) = sum(tx_symbols(i+2:-1:i) .* transfer_function);
end
% 具体来说，在处理长度为 10004 的序列时，从第 1 个元素开始，要一直到第 3 个元素才能和长度为 3 的传递函数完全重叠，
% 而最后一个完全重叠的位置是第 10004 个元素。
% 这样算下来，总共就有 10002 个完全重叠的位置。

c1 = conv(tx_symbols, transfer_function, "valid"); % 验证卷积结果，完全重叠的卷积结果
c2 = conv(tx_symbols, transfer_function, "same"); % 验证卷积结果

No_symbols = No_symbols + 2; % 调整符号计数（包含首尾状态）
%% ---------------------- 维特比译码器初始化 ----------------------
[arr,input_from_states,outputs,states]=MLSE.getViterbi(transfer_function);

%% 维特比译码
% 所需变量初始化（展示）
matrix = inf(length(states), length(states));  % 路径度量矩阵
weight_of_survivor = zeros(1, length(states)); % 幸存路径度量
survivor_paths = zeros(length(states), No_symbols+1); % 幸存路径存储

SNR = 10.^(SNR_db / 10);   % 转换SNR为线性值
N0 = Eb ./ SNR;             % 计算噪声功率
SEP = zeros(1, length(N0)); % 符号错误概率存储
for a = 1:length(N0)
    % 生成AWGN噪声
    noise = sqrt(N0(a))/2 * randn(1, No_symbols);
    rx_signal = channel_output + noise; % 接收信号
    decoded_symbols=MLSE.runViterbi(states,rx_signal,No_symbols,arr,input_from_states,outputs,initialstate);
    SEP(a)=MLSE.getSEP(tx_symbols,decoded_symbols,No_symbols);
end
figure;
semilogy(SNR_db, SEP, 'LineWidth',2);
xlabel("信噪比(dB)", "FontWeight","bold");
ylabel("符号错误概率(SEP)", "FontWeight","bold");
title('不同回溯深度N下的SEP性能比较');
grid on;