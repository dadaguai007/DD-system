function [out,w]=MMSE_Filter(tx,rx,lags)

%%%输入输出为列向量
%%%lags控制w的长度
Z=xcorr(tx,'biased');%输入自相关
Rxx = toeplitz(Z(end/2:end/2+lags));%输入的自相关矩阵
Rxy=xcorr(tx,rx,"biased");%输入输出互相关
Rxy=Rxy(end/2:end/2+lags);%输出的互相关向量
% MMSER准则下的滤波器抽头
w=inv(Rxx)*Rxy;

% 输出
out=conv(tx,w,"valid");

end