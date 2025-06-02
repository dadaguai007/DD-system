function err_avg = cal_lut_avg(err, idx_I, idx_Q, N)
    % Calculate average pattern lookup tables
    %计算平均模式查找表。着眼于I和Q路，计算所有模式的平均误差。索引需要对齐，以便误差和相应的模式具有相同的索引
    % 计算平均模式查找表（LUT）
    % 输入参数:
    %   err   : 误差向量（复数形式，实部为I路误差，虚部为Q路误差）
    %   idx_I : I路模式索引向量（每个元素表示对应误差点的模式索引）
    %   idx_Q : Q路模式索引向量（每个元素表示对应误差点的模式索引）
    %   N     : 查找表的大小（总模式数）
    % 输出:
    %   err_avg : 长度为N的复数向量，包含每个模式的平均误差值
    
    % Get the size of the error array
    L = length(err);
    
    % Initialize arrays
    err_avg_I = zeros(N, 1);
    err_avg_Q = zeros(N, 1);
    nI = zeros(N, 1);
    nQ = zeros(N, 1);
    
    % Loop through each error and corresponding indices
    % err_avg为每个模式的累计误差量，所有误差都经过后，得到累计量
    for i = 1:L
        err_avg_I(idx_I(i)) = err_avg_I(idx_I(i)) + real(err(i));
        err_avg_Q(idx_Q(i)) = err_avg_Q(idx_Q(i)) + imag(err(i));
        nI(idx_I(i)) = nI(idx_I(i)) + 1;
        nQ(idx_Q(i)) = nQ(idx_Q(i)) + 1;
    end
    
    % 某个模式没有出现，要避免除以0的情况发生
    % Ensure that counts are not zero to avoid division by zero
    nI(nI == 0) = 1;
    nQ(nQ == 0) = 1;
    
    % Calculate average errors
    % 计算平均的误差
    err_avg = (err_avg_I./nI) + 1i * (err_avg_Q./nQ);
end
