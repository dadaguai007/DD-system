function [x_hat, no_eff] = zf_equalizer(y, h, s)
    % MIMO 迫零(ZF)均衡器
    % 输入参数同上
    % 输出参数同上
    
    % 计算伪逆矩阵
    g = pinv(h);  % Moore-Penrose 伪逆
    
    % 符号估计
    x_hat = g * y;
    
    % 计算噪声协方差 G *S* G^H
    gsg = g * s * g';
    no_eff = real(diag(gsg));  % 提取对角线元素
end