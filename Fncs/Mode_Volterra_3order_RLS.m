function [ output ] = Mode_Volterra_3order_RLS( input, TrainSeq,lambda, L1,L2,L3,Up2)
% Training sequence based RLS equalization
% Ampilitude_directed 3-order Volttera
% 两种模式的工作滤波器

TrainLen = length(TrainSeq);
DataLen = length(input);
first_nonlinear_length  =  2*L1+1; %一阶长度
second_nonlinear_length = (2*L2+1)*(2*L2+1+1)/2; % 二阶长度
third_memory_length = 2*L3+1;
third_nonlinear_length = 3*(third_memory_length)*(third_memory_length+1)/12+...
    third_memory_length*(third_memory_length+1)*(2*third_memory_length+1)/12;  % 三阶长度

w1 = zeros(1,first_nonlinear_length+second_nonlinear_length+third_nonlinear_length);
w2 = zeros(1,first_nonlinear_length+second_nonlinear_length+third_nonlinear_length);
w1(L1+1) = 1;
w2(L1+1) = 1;
input = [zeros(1,L1) input zeros(1,L1)]; %  对信号进行相应填充， 防止滑动窗口越界
Delta = 0.1*eye(length(w1));   % FFE逆相关矩阵初始化
iter = 2;
error = [];

ADMode = 0;
threshold = median(abs(input));     % 阈值
for kk =1:iter
    for ii = 1:TrainLen
        x1 = input((ii-1)*Up2+1:(ii-1)*Up2+1+L1+L1);
        x2 = input((ii-1)*Up2+1+L1-L2:(ii-1)*Up2+1+L1+L2);
        x3 = input((ii-1)*Up2+1+L1-L3:(ii-1)*Up2+1+L1+L3);
        y1 = x1;
        y2 = two_order(x2); % 二阶信号非线性排列组合
        y3 = third_order(x3); % 三阶信号非线性排列组合
        y=[y1 y2 y3];
        switch  ADMode
            case 0
                err = TrainSeq(ii) - conj(w1)*y.';
                G = Delta*y.' / (lambda + conj(y) *Delta*y.');
                Delta = 1/lambda*(Delta - G*conj(y)*Delta);
                w1 = w1 + G.'.*conj(err);
            case 1
                if    y(L1+1)>0
                    % 强信号滤波器
%                     y(L1+1)>0   abs(y(L1+1))>threshold
                    err = TrainSeq(ii) - conj(w1)*y.';
                    G = Delta*y.' / (lambda + conj(y) *Delta*y.');
                    Delta = 1/lambda*(Delta - G*conj(y)*Delta);
                    w1 = w1 + G.'.*conj(err);
                else
                    % 备用滤波器
                    err = TrainSeq(ii) - conj(w2)*y.';
                    G = Delta*y.' / (lambda + conj(y) *Delta*y.');
                    Delta = 1/lambda*(Delta - G*conj(y)*Delta);
                    w2 = w2 + G.'.*conj(err);
                end
        end
        error = [error err];
    end
end
% figure; plot(abs(error));
% plot(w1);hold on
% plot(w2)
%% equalization
for jj = 1 : DataLen
    x1 = input(jj:jj+L1+L1);
    x2 = input(jj+L1-L2:jj+L1+L2);
    x3 = input(jj+L1-L3:jj+L1+L3);
    y1 = x1;
    y2 = two_order(x2);
    y3 = third_order(x3);
    y=[y1 y2 y3];
    switch  ADMode
        case 0
            output(jj) = conj(w1) * y.';
        case 1
            if  y(L1+1)>0
%                 y(L1+1)>0 or abs(y(L1+1))>threshold
                output(jj) = conj(w1) * y.';
            else
                output(jj) = conj(w2) * y.';
            end
    end
end
end


function Output=two_order(input)

x = input;
L = length(input);
t = 0;
for k=1:L
    for m=k:L
        t=t+1;
        Output(t,:)=x(k)*x(m);
    end
end
Output = Output.';
end

function kernel = third_order(input)
x = input;
L = length(input);
t=0;
% kernel=zeros(455,1);
kernel = [];
for k = 1:L
    for m = k:L
        for n=m:L
            t=t+1;
            kernel(t,:)=x(k)*x(m)*x(n);
        end
        %         kernel = [kernel;x(k)*x(m)*x(m:L)];
    end
end
kernel = kernel.';
end