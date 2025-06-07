function out = matrix_pinv(tensor)
% 计算矩阵的Moore-Penrose伪逆（列满秩情况）
% 输入：
%   tensor : [..., M, K] 实浮点数或复数数组
%       输入数组，最后两维为矩阵维度（M行K列）
% 输出：
%   out : [..., K, M] 实浮点数或复数数组
%       伪逆矩阵 A⁺，满足 A⁺A = I_K

% 获取输入数组尺寸
sz = size(tensor);
M = sz(end-1); % 行数
K = sz(end);   % 列数

% 计算共轭转置 Aᴴ
A_H = pagectranspose(tensor); % 分页共轭转置

% 计算 AᴴA (K×K 矩阵)
AHA = pagemtimes(A_H, tensor); % 分页矩阵乘法

% Cholesky分解（返回上三角矩阵R：AᴴA = R'*R）
R = chol(AHA);

% 解Cholesky方程：AᴴA·X = Aᴴ → R'R·X = Aᴴ
% 等效于TensorFlow的cholesky_solve(l, Aᴴ)
% 1. 解下三角方程：R'·Y = Aᴴ
Y = R' \ A_H; 
% 2. 解上三角方程：R·X = Y
X = R \ Y;

% 恢复原始批量维度（但最后两维从[M,K]变为[K,M]）
out_size = [sz(1:end-2), K, M];
out = reshape(X, out_size);
end
