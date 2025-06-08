function [w,out]=MMSE_Filter(tx,rx,lags)
% MMSE滤波器
% w_mmse = (Ryy)^-1 * Ryd ;
%%%输入输出为列向量
%%%lags控制w的长度，w数量为lags+1；
Z=xcorr(tx,'biased');%输入信号的自相关
Rxx = toeplitz(Z(end/2:end/2+lags));%输入信号自相关矩阵
Rxy=xcorr(tx,rx,"biased");%输入输出互相关
Rxy=Rxy(end/2:end/2+lags);%输入输出的互相关向量
% MMSER准则下的滤波器抽头
w=inv(Rxx)*Rxy;
equalizer_length=length(w);

% 输出
out=conv(tx,w);
out=out(1:end-equalizer_length+1);

end