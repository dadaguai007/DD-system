function wout = pec(win, alpha, constellation)
% 通过分析相邻符号误差方向特征，抑制脉冲噪声引起的系统性符号偏移

% 输入参数:
%   win          - 输入信号向量（含噪声的符号估计值，实数/复数）
%   alpha        - 步长因子（修正强度系数，建议值0.1~0.5）
%   constellation - 星座图（用于硬判决，如QPSK/16QAM的复数向量）
%
% 输出参数:
%   wout         - 修正后的信号向量

% ================== 1. 初始化与第一阶段处理 ==================
Q = decision(win, constellation);     % 初始硬判决：映射到最近星座点
error = win - Q;                      % 计算原始误差向量（估计值-判决值）
error_judge = sign(error);            % 提取误差方向符号(+1/-1)
wout = win;                           % 初始化输出信号

% 遍历中间符号（跳过首尾无双侧邻居的符号）
for i = 2:length(Q)-1
    % 条件1：检查前向(i-1)和后向(i+1)符号的误差方向是否一致
    if error_judge(i-1) == error_judge(i+1)
        % 条件2：当前符号误差方向是否与相邻符号一致？
        if error_judge(i-1) == error_judge(i)
            % 连续三个符号同向误差 → 判定为脉冲噪声干扰
            temp = error(i-1) + error(i+1);  % 叠加相邻误差
        else
            % 当前符号方向与相邻不同 → 判定为正确判决
            temp = 0;               % 不进行修正（避免过校正）
        end
    else
        % 相邻符号误差方向不一致 → 判定为随机噪声
        temp = 0;                   % 不进行修正
    end
    % 加权修正：alpha控制修正强度，/2防止过冲
    wout(i) = wout(i) + alpha*temp/2; 
end

% ================== 2. 第二阶段处理 ==================
Q = decision(wout, constellation);    % 重新硬判决修正后的信号
error = wout - Q;                    % 计算新误差向量
error_judge = sign(error);           % 更新误差方向符号

for i = 2:length(Q) - 1
    % 计算相邻误差的综合方向（向量和符号）
    judge = sign(error(i-1) + error(i+1));
    
    % 判断当前符号误差方向是否与综合方向一致
    if error_judge(i) ~= judge
        temp = 0;                    % 方向不一致 → 不修正
    else
        temp = error(i-1) + error(i+1); % 方向一致 → 叠加相邻误差
    end
    % 二次加权修正
    wout(i) = wout(i) + alpha*temp/2;
end

end
