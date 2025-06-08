function out = inv_cholesky(tensor)
% 计算矩阵的Cholesky分解的逆矩阵
% 输入：
%   tensor : [..., M, M] 实浮点数或复数数组
%       输入数组，最后两个维度必须是方阵（M×M）
% 输出：
%   out : [..., M, M] 实浮点数或复数数组
%       Cholesky分解的逆矩阵（L^{-1}），其中 A = LL^H

% 获取输入数组尺寸
sz = size(tensor);
M = sz(end); % 矩阵维度

% 计算Cholesky分解（返回上三角矩阵R：tensor = R' * R）
R = chol(tensor); % MATLAB的chol返回上三角矩阵，而TensorFlow返回下三角矩阵
% 注意：R' 等价于TensorFlow中的L（下三角矩阵）

% 构造单位矩阵（与输入数据类型一致）
I = eye(M, 'like', tensor);

% 解下三角方程：L*X = I → R'*X = I
% 这里使用左除运算符高效求解
L_inv = R' \ I; % 等价于TensorFlow的triangular_solve(l, rhs, lower=True)

% 恢复原始维度
out = reshape(L_inv, sz);
end

