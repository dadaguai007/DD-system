function [Shaping_Tap, error,output] = PR_FFE_RLS_Train( input, TrainSeq, delay, lambda, Up2, Coe )
% Training sequence based LMS equalization
% 用于白化滤波器的生成
% 精确逆矩阵更新：通过Delta/Delta2维护逆相关矩阵，避免直接矩阵求逆

taps = 2*delay + 1; % 通过延时确定抽头数，以中心抽头为基准
iter = 3;
TrainLen = length(TrainSeq);
DataLen = length(input);
output = zeros(1, DataLen);
w = zeros(1, taps);
w(delay + 1) = 1;
error = [];
input = [zeros(1,delay) input zeros(1,delay)]; % 边界保护
Delta = 0.1*eye(taps,taps); % 逆相关矩阵初始化
Delta2 = 0.1*eye(length(Coe),length(Coe)); % 成形白化滤波器逆相关矩阵初始化
%% training
for kk =1 : iter
for ii = length(Coe)+1 : TrainLen
    % % 步骤1：生成成形后的训练序列
    TrainSeq_PR = filter([1 Coe],1,TrainSeq); 
    % 输入信号
    x = input( (ii-1)*Up2+1: (ii-1)*Up2+taps );
    % 计算
    err = TrainSeq_PR(ii) - conj(w) * x.';  
    % 用于白化滤波器的更新参数
    temp = TrainSeq_PR(ii-1:-1:ii-length(Coe));
    
    % RLS增益计算（双模块）
    G = Delta*x.' / (lambda + conj(x) *Delta*x.'); % FFE 增益
    Delta = 1/lambda*(Delta - G*conj(x)*Delta); % FFE矩阵更新

    G2 = Delta2*temp.' / (lambda + conj(temp) *Delta2*temp.'); % 白化滤波器增益
    Delta2 = 1/lambda*(Delta2 - G2*conj(temp)*Delta2); % 白化矩阵更新

    w = w + G.'.*conj(err);
    Coe = Coe - G2.'.*conj(err);
    error = [error err];
end
end
Shaping_Tap = [1 Coe]; % 优化发射端波形，减少码间干扰(ISI)，提升眼图质量
% figure;
% plot(abs(error));

%% equalization
for jj = 1 : DataLen
    in = input(jj : jj + taps - 1);
    output(jj) = conj(w) * in.';
end

end

