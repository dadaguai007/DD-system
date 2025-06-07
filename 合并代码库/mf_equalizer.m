function [x_hat, no_eff] = mf_equalizer(y, h, s)
    % MIMO 匹配滤波器(MF)均衡器
    % 输入参数同上
    % 输出参数同上
    
    % 计算 H^H H
    hth = h' * h;
    
    % 构建对角缩放矩阵 diag(H^H * H)^-1
    d_inv = diag(1 ./ diag(hth));
    g = d_inv * h';  % G = diag(H^H H)^-1 * H^H
    
    % 符号估计
    x_hat = g * y;
    
    % 计算噪声相关项
    gsg = g * s * g';  % G S G^H
    gh = g * h;        % G H
    I = eye(size(gh)); % 单位矩阵
    
    % 合并噪声项并提取方差
    noise_term = (I - gh) * (I - gh)' + gsg;
    no_eff = abs(diag(noise_term));
end