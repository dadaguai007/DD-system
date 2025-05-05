ML = 101;

constellation = unique(TxSym);
alpha = 0:0.01:1;

sndsym=RxDown; % 接收信号 为滤波器后的输出 未进行下采样 为：2sps
sndsym = sndsym/sqrt(mean(abs(sndsym).^2))*sqrt(5);
pncsym=RxDown;
pncsym = pncsym/sqrt(mean(abs(pncsym).^2))*sqrt(5);

% parfor i = 1:length(alpha)
%
%     pf_out = mlseeq(filter([1, alpha(i)], 1, sndsym), [1, alpha(i)], constellation, 1000, 'rst');
%     % [SNR_snd, noise_snd] =SNR_Calculation_v2(pec_out,TxSym) ;
%     % disp(SNR_snd);
%     M=2^BitPerSym;
%     pf_out = pamdemod(pf_out,M);
%     pfBit= pam4demod(pf_out);
%     pfBit = reshape(pfBit, 1,[]);
%     [~, BER_pf(i), ~] = biterr(pfBit,TxBit(ML*2-1:length(RxBit)+ML*2-2));
%     % disp(BER_snd);
% end
pec_out = my_pec(sndsym, .36, constellation);
% [SNR_snd, noise_snd] =SNR_Calculation_v2(pec_out,TxSym) ;
% disp(SNR_snd);
M=2^BitPerSym;
pec_out = pamdemod(pec_out ,M);
sndBit= pam4demod(pec_out);
sndBit = reshape(sndBit, 1,[]);
[~, BER_pec, ~] = biterr(sndBit,TxBit(ML*2-1:length(RxBit)+ML*2-2));
% disp(BER_snd);

pnc_out_1 = my_pnc(pncsym, .29, constellation, 1);
% [SNR_pnc, noise_pnc] =SNR_Calculation_v2(pnc_out_1,TxSym) ;
% disp(SNR_pnc);
M=2^BitPerSym;
pnc_out_1 = pamdemod(pnc_out_1,M);
pncBit_1 = pam4demod(pnc_out_1);
pncBit_1 = reshape(pncBit_1, 1,[]);
[~, BER_pnc_1, ~] = biterr(pncBit_1,TxBit(ML*2-1:length(RxBit)+ML*2-2));
% disp(BER_pnc);

pnc_out_2 = my_pnc(pncsym, .11, constellation, 2);
% [SNR_pnc, noise_pnc] =SNR_Calculation_v2(pnc_out_1,TxSym) ;
% disp(SNR_pnc);
M=2^BitPerSym;
pnc_out_2 = pamdemod(pnc_out_2,M);
pncBit_2= pam4demod(pnc_out_2);
pncBit_2 = reshape(pncBit_2, 1,[]);
[~, BER_pnc_2, ~] = biterr(pncBit_2,TxBit(ML*2-1:length(RxBit)+ML*2-2));
% disp(BER_pnc);

pnc_out_3 = my_pnc(pncsym, .05, constellation, 3);
% [SNR_pnc, noise_pnc] =SNR_Calculation_v2(pnc_out_1,TxSym) ;
% disp(SNR_pnc);
M=2^BitPerSym;
pnc_out_3 = pamdemod(pnc_out_3,M);
pncBit_3= pam4demod(pnc_out_3);
pncBit_3 = reshape(pncBit_3, 1,[]);
[~, BER_pnc_3, ~] = biterr(pncBit_3,TxBit(ML*2-1:length(RxBit)+ML*2-2));
% disp(BER_pnc);