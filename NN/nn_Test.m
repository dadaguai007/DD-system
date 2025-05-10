% 神经网络均衡训练
% 采样率为1sps


%% Downsample  通过均衡后的数据，下采样，绘制强度分布图
RxDown = RxEqu(1:Up2:end);
scatterplot(RxDown)
amplitude_distribution(RxDown);
Up=1;
%% NN
if (eqmodei==100)
    %             TrainLen=40000;Up2=1;Inputneuron=51;Hiddenneuron=40;
    %             maxepoch=30;
    %             [NNinput_train] = SlidingSample_MD_XZP(RxDown,Inputneuron,TrainLen,Up2);
    %             NNoutput_train = TxSym(1:TrainLen);
    %             [NNinput_test] = SlidingSample_MD_XZP(RxDown,Inputneuron,FrmLen,Up2);
    %             %     [net,tr,NNoutput_train_predict,NNoutput_train,NNoutput_test_predict,NNoutput_test,Train_BER_N_input,Train_MSE_N_input,Test_BER_N_input,Test_MSE_N_input]=...
    %             %         NNfunc_2layerCFNN_PAM4(r_downsample,t,60000,20000,0,30000,25,21,'trainlm','tansig','purelin',50,2e-6,2e6,6,0)
    %             net=feedforwardnet(Hiddenneuron,'trainlm');% try trainscg(Scaled Conjugate Gradient) and trainrp (Resilient Backpropagation)
    %             %             net=cascadeforwardnet(N_hidden(kk),'trainlm');
    %             %             if (flag_C_FNN==1)
    %             %                 net.inputConnect=[1;1];
    %             %             end
    %
    %             net=configure(net,NNinput_train,NNoutput_train);
    %             net.layers{1}.transferFcn='tansig';
    %             %             net.layers{2}.transferFcn='tansig';
    %             net.layers{2}.transferFcn='purelin';
    %             net=init(net);
    %             view(net);
    %             net.trainParam.epochs=ii;
    %             net.trainParam.goal=2.3e-4;
    %             net.trainParam.min_grad=2e-6;     %for trainlm (may not have sufficient data)
    %             net.trainParam.mu_max=2e6;
    %             net.trainParam.max_fail=4;
    %             %             net.trainParam.showCommandLine=1;
    %             %             net.trainParam.show=10;
    %             net.divideFcn='dividerand';
    %
    %             [net,tr]=train(net,NNinput_train,NNoutput_train,'useParallel','yes');
    %             %             figure(); plotperform(tr);
    %             %             disp('N_input:');disp(N_input(i));%disp(' '); disp(tr.num_epochs);
    %             %             Train_MSE_N_input(i,kk)=tr.best_perf;
    %             %             disp('Training MSE: '); disp(Train_MSE_N_input(i));
    %             NNoutput_test_predict=sim(net,NNinput_test,'useParallel','yes'); %1*N_train
    %             RxSym = NNoutput_test_predict;

    NNnodes = 51;
    TrainLen = 10000;
    TrainData =SlidingSample_MD(RxDown,NNnodes,TrainLen,Up);% 在1sps下均衡
    TestData = SlidingSample_MD_forTestData(RxDown,NNnodes,FrmLen);
    TrainSeq = [TxSym(1:TrainLen)];
    [TrainLabel] = LabelGeneration(TrainSeq);
    options = trainingOptions('adam', ... %优化器
        'InitialLearnRate',1e-3, ....% 学习率
        'LearnRateSchedule','piecewise', ... % 学习策略。每N个epoch，就将学习率×factor
        'LearnRateDropFactor',0.5, ... %M
        'LearnRateDropPeriod',10, ...%N
        'MaxEpochs',50, ... % 最大训练集遍历数
        'MiniBatchSize',100, ... % 每batch更新一次网络参数                    'Shuffle','every-epoch',...
        'ExecutionEnvironment','cpu',... % CPU or GPU
        'Plots','training-progress'); % 画图
    %                                 'ValidationFrequency',50,...
    %                                 'ValidationData',ValidationData,...

    [trainnet,info]= trainNetwork(TrainData,TrainLabel,FNN(NNnodes),options);
    save ('parameters.mat','trainnet');
    RxEqu = predict(trainnet,TestData);
    temp0 = cell2mat(RxEqu);
    RxDown = temp0.';
end