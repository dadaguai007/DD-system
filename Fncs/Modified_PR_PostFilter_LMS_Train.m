function [Shaping_Tap, error] = Modified_PR_PostFilter_LMS_Train( input, TrainSeq, step, Coe )
% Training sequence based LMS equalization
% 归一化白化滤波器参数
iter = 1;
TrainLen = length(TrainSeq);
error = [];
%% training
for kk =1 : iter
for ii = length(Coe)+1 : TrainLen
    TrainSeq_PR = filter([1 Coe],1,TrainSeq); %%
    input_PR = filter([1 Coe],1,input);
    temp = TrainSeq(ii-1:-1:ii-length(Coe))-input(ii-1:-1:ii-length(Coe));
    err = TrainSeq_PR(ii) - input_PR(ii);  
    % 白化滤波器的参数进行归一化
    Coe = Coe - step* (conj(err).*temp)./(TrainSeq(ii).^2); %% 符号和W相反
    error = [error err];
end
end
Shaping_Tap = [1 Coe];
% figure;
% plot(abs(error));
% plot(abs(error));

end

