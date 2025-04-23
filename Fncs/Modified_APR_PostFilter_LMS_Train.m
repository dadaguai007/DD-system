function [Shaping_Tap, error] = Modified_APR_PostFilter_LMS_Train( input, TrainSeq, step, Coe )
% Training sequence based LMS equalization
% ZhuYixiao @ 2019.8.16
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
    Coe = Coe - step* (conj(err).*temp)./(TrainSeq(ii).^2); %% 符号和W相反
    error = [error err];
end
end
Shaping_Tap = [1 Coe];
figure;
plot(abs(error));
% plot(abs(error));

end

