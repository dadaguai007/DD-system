function [output1,output2] = DD_RLS(input1,input2,delay,lambda,Up2,BitPerSym)
% 交叉干扰建模；双偏振（DP-QPSK/16QAM）传输

Len = length(input1);
taps = 2*delay + 1; % 抽头数
w1 = zeros(1,taps*2);
w2 = zeros(1,taps*2);
w1(delay+1) = 1;
w2(delay+1) = 1;
input1 = [zeros(1,delay) input1 zeros(1,delay)]; % 防止边界溢出
input2 = [zeros(1,delay) input2 zeros(1,delay)]; % 防止边界溢出
% lambda = 0.99;
Delta1 = 0.1*eye(taps*2,taps*2);
Delta2 = 0.1*eye(taps*2,taps*2);

output1 = zeros(1,Len);
output2 = zeros(1,Len);

for ii = 1:Len
    x1 = [input1( (ii-1)*Up2+1:(ii-1)*Up2+taps ) input2( (ii-1)*Up2+1:(ii-1)*Up2+taps )];
    x2 = [input2( (ii-1)*Up2+1:(ii-1)*Up2+taps ) input1( (ii-1)*Up2+1:(ii-1)*Up2+taps )];
    output1(ii) = conj(w1)*x1.';
    output2(ii) = conj(w2)*x2.';
    TrainSeq1(ii) = QAMdecision(output1(ii),BitPerSym);
    TrainSeq2(ii) = QAMdecision(output2(ii),BitPerSym);
    err1 = TrainSeq1(ii) -  output1(ii);
    err2 = TrainSeq2(ii) -  output2(ii);
    
    G1 = Delta1*x1.' / (lambda + conj(x1) *Delta1*x1.');
    Delta1 = 1/lambda*(Delta1 - G1*conj(x1)*Delta1);
    G2 = Delta2*x2.' / (lambda + conj(x2) *Delta2*x2.');
    Delta2 = 1/lambda*(Delta2 - G2*conj(x2)*Delta2);
    w1 = w1 + G1.'.*conj(err1);
    w2 = w2 + G2.'.*conj(err2);
    error1(ii) = err1;
    error2(ii) = err2;
end

% figure;plot(abs(error1)); hold on;plot(abs(error2),'r');hold off;


end

