function [x_hat, no_eff] = lmmse_equalizer(y, h, s, whiten_interference)
    % MIMO LMMSE 均衡器
    % 输入:
    %   y : M×1 复数向量，接收信号
    %   h : M×K 复数矩阵，信道矩阵
    %   s : M×M 复数矩阵，噪声协方差矩阵
    %   whiten_interference : 是否进行白化处理
    % 输出:
    %   x_hat : K×1 复数向量，估计符号
    %   no_eff : K×1 实数向量，有效噪声方差
    
    if whiten_interference
        % 噪声白化处理
        [y, h] = whiten_channel(y, h, s);
        s = [];  % 白化后噪声协方差变为单位矩阵
    end
    
    % 计算均衡矩阵
    g = lmmse_matrix(h, s);
    
    % 计算 Gy (均衡后信号)
    gy = g * y;
    
    % 计算 GH (信道-均衡器组合)
    gh = g * h;
    
    % 提取对角线元素并求逆
    d = diag(gh);
    x_hat = gy ./ d;  % 符号估计
    
    % 计算有效噪声方差
    no_eff = real(1 ./ d - 1);
end