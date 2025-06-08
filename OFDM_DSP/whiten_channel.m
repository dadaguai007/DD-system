
function [y, h] = whiten_channel(y, h, s)
% 信道白化
%     输入:
%       y : 接收信号
%       h : 信道矩阵
%       s : 噪声协方差矩阵
%     输出:
%       y : 白化后的接收信号
%       h : 白化后的信道矩阵
% 计算白化矩阵 (Cholesky分解的逆)

% L_inv = inv_cholesky(s);
L = chol(s, 'lower');  % Cholesky 分解
% L_inv = inv(L);
% 计算下三角矩阵的逆
L_inv = L \ eye(size(s));
y = L_inv * y;
h = L_inv * h;
end
