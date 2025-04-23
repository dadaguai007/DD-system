function [X] = signal2mat(x,N)
%SIGNAL2MAT 此处显示有关此函数的摘要
%   此处显示详细说明
% 目的是将一个信号x转换为一个矩阵，
% 其中矩阵的每一列都是信号x中的一个长度为N的连续子序列
X=zeros(N,length(x)-N+1);
for ii=1:(length(x)-N+1)
    X(:,ii)=x(ii:ii+N-1);
end
end

