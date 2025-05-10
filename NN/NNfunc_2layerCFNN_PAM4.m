function [net,tr,NNoutput_train_predict,NNoutput_train,NNoutput_test_predict,NNoutput_test,Train_BER_N_input,Train_MSE_N_input,Test_BER_N_input,Test_MSE_N_input]=NNfunc_2layerCFNN_PAM4(r_downsample,t,I_N_total,I_N_train,I_N_validation,I_N_test,Nin,Nhid,TrainMethod,L1TF,L2TF,NE,NGR,NMU,NMAXF,IFcascade)

% PAM4 2-layer (Cascade) FNN /can be generated to OOK, expect the demoduation part
% Input: r_downsample  received downsampled symbols
%        t  transmitted symbols (correspond to r_downsample)
%        I_N_total,I_N_train,I_N_validation,I_N_test  Number of total/train/validation/test symbols   
%        Nin,Nhid  Number of input/hidden neurons
%        TrainMethod  Traing method ,e.g., 'trainlm'
%        L1TF,L2TF  Transfer Function of layer 1 and
%        NE,NGR,NMU,NMAXF  epochs/min_grad/mu_max/max_fail (NGR for trainlm) 
%        IFcascade  set as 1 if NN has cascade structure 

% Output: net,tr  trained net and the corresponding parameter 
%         NNoutput_train_predict,NNoutput_train  predicted/real symbols for training
%         NNoutput_test_predict,NNoutput_test  predicted/real symbols for testing
%         Train_BER_N_input,Train_MSE_N_input  BER/MSE for training symbols
%         Test_BER_N_input,Test_MSE_N_input  BER/MSE for testing symbols

flag_FNN=1;flag_C_FNN=IFcascade;flag_RNN_self=0;flag_C_RNN=0;
% flag_FNN4=0;flagRNN=1; flagfeedback=1;
N_total=I_N_total; % #PAM4 symbols
N_train=I_N_train;
N_validation=I_N_validation;
N_test=I_N_test;


%% F-NN
%% Prepare received symbol for NN
r_pd_real=r_downsample;
[r_pd,ps]=mapminmax(r_pd_real);
 
%% NN performance: MSE versus N_input&N_hidden (2-layer)
if (flag_FNN==1)
    pool = gcp;pool.NumWorkers;
    N_input=Nin:2:Nin;
    N_hidden=Nhid:1:Nhid;
    % create MSE BER and epoch matrices for recording
    Train_MSE_N_input=zeros(length(N_input),length(N_hidden));Validation_MSE_N_input=zeros(length(N_input),length(N_hidden));Test_MSE_N_input=zeros(length(N_input),length(N_hidden));
    Train_BER_N_input=zeros(length(N_input),length(N_hidden));Validation_BER_N_input=zeros(length(N_input),length(N_hidden));Test_BER_N_input=zeros(length(N_input),length(N_hidden));
    Num_epochs=zeros(length(N_input),length(N_hidden));
    for i=1:1:length(N_input);
%% Creat data using different N_input
        NNinput_train=zeros(N_input(i),N_train);% Network input (Received data) N_input*N_train
        for ii=1:1:N_train
            NNinput_train(:,ii)=r_pd(ii:ii+N_input(i)-1).';
        end
        NNoutput_train=t((1+N_input(i))/2:N_train+(1+N_input(i))/2-1);% Network output/Target output (Transmitted data) predict the middle
%         NNoutput_train=t((N_input(i):N_train+N_input(i)-1));% predict the last
%         NNoutput_train=t((1:N_train+1-1));% predict the first
        NNinput_validation=zeros(N_input(i),N_validation);
        for ii=1:1:N_validation
            NNinput_validation(:,ii)=r_pd(N_train+ii:N_train+ii+N_input(i)-1).';
        end
        NNoutput_validation=t(N_train+(1+N_input(i))/2:N_train+(1+N_input(i))/2+N_validation-1);
%         NNoutput_validation=t(N_train+N_input(i):N_validation+N_train+N_input(i)-1);%last
%         NNoutput_validation=t(N_train+1:N_validation+N_train+1-1);%first
        NNinput_test=zeros(N_input(i),N_test);
        for ii=1:1:N_test
            NNinput_test(:,ii)=r_pd(N_train+N_validation+ii:N_train+N_validation+ii+N_input(i)-1).';
        end
        NNoutput_test=t(N_train+N_validation+(1+N_input(i))/2:N_train+N_validation+(1+N_input(i))/2+N_test-1); 
%         NNoutput_test=t(N_train+N_validation+N_input(i):N_test+N_train+N_validation+N_input(i)-1); %last
%         NNoutput_test=t(N_train+N_validation+1:N_test+N_train+N_validation+1-1); %first
        for kk=1:1:length(N_hidden)
%% NN setting
            net=feedforwardnet(N_hidden(kk),TrainMethod);% try trainscg(Scaled Conjugate Gradient) and trainrp (Resilient Backpropagation)          
%             net=cascadeforwardnet(N_hidden(kk),'trainlm');
            if (flag_C_FNN==1)
                net.inputConnect=[1;1];
            end
            
            net=configure(net,NNinput_train,NNoutput_train);
            net.layers{1}.transferFcn=L1TF;
%             net.layers{2}.transferFcn='tansig';
            net.layers{2}.transferFcn=L2TF;
            net=init(net);
%             view(net);
            net.trainParam.epochs=NE;
%             net.trainParam.goal=2.3e-4; 
            net.trainParam.min_grad=NGR;     %for trainlm (may not have sufficient data)
            net.trainParam.mu_max=NMU; 
            net.trainParam.max_fail=NMAXF; 
        %     net.trainParam.showCommandLine=1;
        %     net.trainParam.show=10;
            net.divideFcn='dividerand';

            [net,tr]=train(net,NNinput_train,NNoutput_train,'useParallel','yes');
%             figure(); plotperform(tr);
        %     disp('N_input:');disp(N_input(i));%disp(' '); disp(tr.num_epochs);
%             Train_MSE_N_input(i,kk)=tr.best_perf;
%             disp('Training MSE: '); disp(Train_MSE_N_input(i));
%%  Cal train BER and MSE         
            NNoutput_train_predict=sim(net,NNinput_train,'useParallel','yes'); %1*N_train

            NNoutput_train_predict_pam4_withoutdc=NNoutput_train_predict-0;%1*N_train -3,-1,1,3
            h=modem.pamdemod('M',4, 'SymbolOrder', 'gray', 'OutputType', 'Bit');
            NNoutput_train_predict_binarysequence=demodulate(h,NNoutput_train_predict_pam4_withoutdc.').';
            
            NNoutput_train_pam4_withoutdc=NNoutput_train-0; %1*N_train, -3,-1,1,3 Real (Transmitted data)
            h=modem.pamdemod('M',4, 'SymbolOrder', 'gray', 'OutputType', 'Bit');
            NNoutput_train_binarysequence=demodulate(h,NNoutput_train_pam4_withoutdc.').';    
            Train_BER_N_input(i,kk)=sum(abs(NNoutput_train_predict_binarysequence-NNoutput_train_binarysequence))/N_train/2; 
            Train_MSE_N_input(i,kk)=sum(sum((NNoutput_train_predict-NNoutput_train).^2))/N_train;
%% Cal test BER and MSE
            NNoutput_test_predict=sim(net,NNinput_test,'useParallel','yes');

            NNoutput_test_predict_pam4_withoutdc=NNoutput_test_predict-0;
            h=modem.pamdemod('M',4, 'SymbolOrder', 'gray', 'OutputType', 'Bit');
            NNoutput_test_predict_binarysequence=demodulate(h,NNoutput_test_predict_pam4_withoutdc.').';

            NNoutput_test_pam4_withoutdc=NNoutput_test-0;
            h=modem.pamdemod('M',4, 'SymbolOrder', 'gray', 'OutputType', 'Bit');
            NNoutput_test_binarysequence=demodulate(h,NNoutput_test_pam4_withoutdc.').';    
            Test_BER_N_input(i,kk)=sum(abs(NNoutput_test_predict_binarysequence-NNoutput_test_binarysequence))/N_test/2;   
%             disp('Test BER: '); disp(Test_BER_N_input(i,kk));
            Test_MSE_N_input(i,kk)=sum(sum((NNoutput_test_predict-NNoutput_test).^2))/N_test;
        %     disp('Test MSE: '); disp(Test_MSE_N_input(i));
            Num_epochs(i,kk)=tr.num_epochs;
%             Test_BER_dBm(i_dBm)=Test_BER_N_input(i,kk);
        end
    end
end

