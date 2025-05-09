% PAM Equlation
EQ=struct();
EQ.u=0.015;
EQ.k1=31;
EQ.k2=15;
EQ.ref=8;
EQ.sps=2;
EQ.lamda=0.9999;
ref_seq = repmat(label,1,100);
ref_seq = ref_seq(:);
xn=resample(sigRx_E,2,SpS);
% FFE_LMS
[sigRx_E,en,w] = FFE_LMS(EQ, xn, ref_seq);
% FFE_RLS
% [sigRx_E,en,w] = FFE_RLS(EQ, xn, ref_seq);
% DFE_LMS
% [sigRx_E,en,w] = DFE_LMS(EQ, xn, ref_seq);
% DFE_RLS
% [sigRx_E,en,w] = DFE_RLS(EQ, xn, ref_seq);

% % Volterra
sps=2;
ref=8;
taps_list = [31 5 0 15 3 0];
taps_list_ABS = [31 0 0 15 5 3];
step_len = 0.01;
lamda=0.9999;
% volterra_lms
% [sigRx_E,en,w]=volterra_dfe_lms(xn,ref_seq,sps,ref,taps_list,step_len);
% volterra_rls
% [sigRx_E,en,w]=volterra_dfe_rls(xn,ref_seq,sps,ref,taps_list,lamda);
% abf_rls
% [sigRx_E,en,w]=abf_lms(xn,ref_seq,sps,ref,taps_list_ABS,step_len);

% pwl
% [sigRx_E,en,w] = pwl(xn,ref_seq,sps,ref,15,0.005,[-0.7661,0,0.7661]);

figure;hold on
plot(xn(1:length(label)),'.')
plot(sigRx_E,'k.')
plot(label,'.')
legend('接收信号','均衡后信号','发送信号')

figure
stem(w(:))
title("均衡器抽头")

figure;
semilogy(abs(en(1:1e4)).^2)
xlabel("迭代次数")
ylabel("误差")

[sigRx_E,ref_sync,ff_index] = sync(sigRx_E,label);