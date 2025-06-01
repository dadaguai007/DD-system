function [Shaping_Tap, error] = Modified_PR_PostFilter_RLS_Train( input, TrainSeq, lambda, Coe )
% Training sequence based LMS equalization
% ZhuYixiao @ 2019.8.16
iter = 3;
TrainLen = length(TrainSeq);
error = [];
Delta2 = 0.1*eye(length(Coe),length(Coe));
%% training
for kk =1 : iter
for ii = length(Coe)+1 : TrainLen
     index = find ((Coe)>1);
    Coe(index) =1; 
    index = find ((Coe)<-1);
    Coe(index) =-1; 
    TrainSeq_PR = filter([1 Coe],1,TrainSeq); %%
    input_PR = filter([1 Coe],1,input);
    err= TrainSeq_PR(ii)-input_PR(ii);
    % 进行归一化
    temp = (TrainSeq_PR(ii-1:-1:ii-length(Coe))-input_PR(ii-1:-1:ii-length(Coe)))./input_PR(ii-1:-1:ii-length(Coe));
    G2 = Delta2*temp.' / (lambda + conj(temp) *Delta2*temp.');
    Delta2 = 1/lambda*(Delta2 - G2*conj(temp)*Delta2);
    Coe = Coe - G2.'.*conj(err);
    error = [error err];
end
end
Shaping_Tap = [1 Coe];
figure;
plot(abs(error));
% plot(abs(error));

end

