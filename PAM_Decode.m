% 解码，计算误码
% 重新量化
A1=[-2 0 2];
% 参考序列
[~,label1] = quantiz(label,A1,[-3,-1,1,3]);
label1=decoding(label1,'gray');
% 接收序列
[~,I] = quantiz(sigRx_E,A1,[-3,-1,1,3]);
I=decoding(I,'gray');
ncut=1e5;
[ber,num,error_location] = CalcBER(I(ncut:end),label1(ncut:end)); %计算误码率
fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);