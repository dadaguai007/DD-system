offset=0;
% normalize
sigRx_E = sigRx_E - mean(sigRx_E);
sigRx_E = sigRx_E./max(sigRx_E);
% make decision and convert into 0,1 sequence
out = sign(sigRx_E-offset);
out = 0.5*(out+1);

ncut=1e5;
[ber,num,error_location] = CalcBER(out(ncut:end),ref_seq(ncut:end)); %计算误码率
fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);