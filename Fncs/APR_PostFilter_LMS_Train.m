function [Shaping_Tap, error] = APR_PostFilter_LMS_Train( input, TrainSeq, step, Coe )
% APR_PostFilter_LMS_Train - 基于LMS算法的自适应后置滤波器训练函数
% 输入参数：
%   input     : 接收端原始信号（含信道失真/噪声）
%   TrainSeq  : 已知的理想训练序列（参考信号）
%   step      : LMS算法步长（控制收敛速度）
%   Coe       : 滤波器系数初始值（不含首项1的系数数组）
% 输出参数：
%   Shaping_Tap : 训练完成的后置滤波器系数（结构为[1 Coe]的FIR滤波器）
%   error       : 训练过程中的误差序列（用于监控收敛）

%% 参数初始化
iter = 3;               % 设置全数据集遍历迭代次数（默认3轮）
TrainLen = length(TrainSeq); % 获取训练序列长度
error = [];             % 初始化误差记录数组

%% 主训练循环（迭代优化系数）
for kk = 1 : iter       % 多轮迭代提升收敛稳定性
    for ii = length(Coe)+1 : TrainLen % 遍历有效样本（考虑滤波器阶数）
        %% 信号预处理：对理想序列和接收信号进行滤波
        % 用当前系数生成参考信号和目标信号的滤波版本
        % [1 Coe]构成滤波器系数（首项固定为1的FIR结构）
        TrainSeq_PR = filter([1 Coe], 1, TrainSeq); % 理想序列滤波
        input_PR = filter([1 Coe], 1, input);       % 接收信号滤波

        %% 误差计算
        err = TrainSeq_PR(ii) - input_PR(ii); % 计算当前时刻的瞬时误差

        %% 梯度估计与系数更新
        % 获取历史差值：取前length(Coe)个时刻的理想序列与接收信号差值
        temp = TrainSeq(ii-1:-1:ii-length(Coe)) - input(ii-1:-1:ii-length(Coe));
        % LMS系数更新公式：Coe = Coe - μ * e(n) * x(n)
        % 其中x(n)为输入向量（此处使用历史差值作为梯度估计）
        Coe = Coe - step * (conj(err) .* temp); 

        %% 记录误差
        error = [error err]; % 累积误差序列用于可视化分析
    end
end

%% 后处理
Shaping_Tap = [1 Coe];  % 构造完整滤波器系数（首项1 + 训练得到的Coe）

%% 误差曲线可视化（调试用）
figure;
plot(abs(error));       % 绘制误差绝对值变化曲线
xlabel('迭代次数');
ylabel('误差幅度');
title('LMS训练误差收敛曲线');
grid on;
end