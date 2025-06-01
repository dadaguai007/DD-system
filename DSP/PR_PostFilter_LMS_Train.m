function [Shaping_Tap, error] = PR_PostFilter_LMS_Train( input, TrainSeq, step, Coe )
% Training sequence based LMS equalization

iter = 3;
TrainLen = length(TrainSeq);
error = [];
%% training
for kk =1 : iter
for ii = length(Coe)+1 : TrainLen
    % 步骤1：生成理想成形信号
    TrainSeq_PR = filter([1 Coe],1,TrainSeq); %%
    % 步骤2：接收信号后置滤波
    input_PR = filter([1 Coe],1,input);
    % 步骤3：计算误差
    err = TrainSeq_PR(ii) - input_PR(ii);  
    % 梯度估计 引入符号差异信息加速收敛
    temp = TrainSeq(ii-1:-1:ii-length(Coe))-input(ii-1:-1:ii-length(Coe));
    Coe = Coe - step* (conj(err).*temp); %% 符号和W相反
    error = [error err];
end
end
Shaping_Tap = [1 Coe];
% figure;
% plot(abs(error));
% plot(abs(error));

end

