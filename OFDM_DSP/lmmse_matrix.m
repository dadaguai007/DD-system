function g = lmmse_matrix(h, s)
    % MIMO LMMSE 均衡矩阵计算
    % 输入:
    %   h : M×K 复数矩阵，信道矩阵
    %   s : M×M 复数矩阵，噪声协方差矩阵 (若为 [] 表示单位矩阵)
    % 输出:
    %   g : K×M 复数矩阵，LMMSE 均衡矩阵
    
    if isempty(s)
        % 白噪声情况 (S=I)
        % 计算 G = (H^H *H + I)^-1 H^H
        hhs = h' * h + eye(size(h, 2));
        g = hhs \ h';  % 等价于 inv(hhs)*h'
    else
        % 通用噪声协方差情况
        % 计算 G = H^H (H*H^H + S)^-1
        hhs = h * h' + s;
        g_t = hhs \ h;  % 求解 (H H^H + S) G_t = H
        g = g_t';       % G = G_t^H
    end
end

% LMMSE 均衡：
% 最小化均方误差，平衡干扰消除和噪声增强
% 通过 Cholesky 分解高效求解矩阵逆
% 支持噪声白化预处理提升数值稳定性
% ZF 均衡：
% 完全消除干扰，但放大噪声
% 使用 Moore-Penrose 伪逆求解
% 计算复杂度低于 LMMSE
% MF 均衡：
% 最大化信噪比，但忽略干扰
% 通过信道共轭转置实现简单匹配
% 适合低干扰场景