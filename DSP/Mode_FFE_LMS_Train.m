function [ w1,output ] =Mode_FFE_LMS_Train( input, TrainSeq, delay, step, Up2 )
% Training sequence based LMS equalization
% 支持两种工作模式，通过幅度阈值判断优化弱信号的处理
taps = 2*delay + 1;  % 滤波器延迟参数（决定抽头数）
iter = 10;
TrainLen = length(TrainSeq);
DataLen = length(input);
output = zeros(1, DataLen);
w1 = zeros(1, taps);
w2 = zeros(1, taps);
w1(delay + 1) = 1;
w2(delay + 1) = 1;
error = [];
input = [zeros(1,delay) input zeros(1,delay)];   % 对信号进行相应填充

ADMode = 1;
% % 计算输入信号幅度中位数作为阈值，用于区分强弱信号
threshold = median(abs(input));
%% training
for kk =1 : iter
    for ii = 1 : TrainLen
        x = input( (ii-1)*Up2+1: (ii-1)*Up2+taps );
        switch  ADMode
            case 0
                % 标准LMS算法
                err = TrainSeq(ii) - conj(w1) * x.';
                w1 = w1 + step * conj(err) .* x;
                error = [error err];
            case 1  % 模式1：带幅度门限的自适应LMS
                if abs(x(delay+1))>threshold
                    % 强信号区域：使用主滤波器w1
                    err = TrainSeq(ii) - conj(w1) * x.';
                    w1 = w1 + step * conj(err) .* x;
                    error = [error err];
                else
                    % 弱信号区域：使用备用滤波器w2
                    err = TrainSeq(ii) - conj(w2) * x.';
                    w2 = w2 + step * conj(err) .* x;
                end

        end

    end
end
% figure; plot(abs(error));
%% equalization
for jj = 1 : DataLen
    in = input(jj : jj + taps - 1);
    switch  ADMode
        case 0
            output(jj) = conj(w1) * in.';
        case 1
            if abs(x(delay+1))>threshold
                output(jj) = conj(w1) * in.';
            else
                output(jj) = conj(w2) * in.';
            end
    end
end

