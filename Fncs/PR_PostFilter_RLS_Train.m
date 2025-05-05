function [Shaping_Tap, error] = PR_PostFilter_RLS_Train( input, TrainSeq, lambda, Coe )
% Training sequence based LMS equalization
% RLS算法应用
iter = 3;
TrainLen = length(TrainSeq);
error = [];
Delta2 = 0.1*eye(length(Coe),length(Coe));
%% training
for kk =1 : iter
for ii = length(Coe)+1 : TrainLen
    % 白化滤波器的默认值
     index = find ((Coe)>1);
    Coe(index) =1; 
    index = find ((Coe)<-1);
    Coe(index) =-1; 
    % 通过白化滤波的理想信号 和 接收信号
    TrainSeq_PR = filter([1 Coe],1,TrainSeq); %%
    input_PR = filter([1 Coe],1,input);
    err = TrainSeq_PR(ii) - input_PR(ii);  
    temp = TrainSeq(ii-1:-1:ii-length(Coe))-input(ii-1:-1:ii-length(Coe));
    % RLS增益
    G2 = Delta2*temp.' / (lambda + conj(temp) *Delta2*temp.');
    % 白化滤波器矩阵
    Delta2 = 1/lambda*(Delta2 - G2*conj(temp)*Delta2);
    Coe = Coe - G2.'.*conj(err);
    error = [error err];
end
end
Shaping_Tap = [1 Coe];
% figure;
% plot(abs(error));
% plot(abs(error));

end

