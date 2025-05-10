function wout = pnc_enhanced(win, alpha, constellation)

% 通过双阶段混合条件修正策略，抑制脉冲噪声和系统性符号偏移

% 输入参数:
%   win          - 输入信号向量（含噪声的符号估计值，实数/复数）
%   alpha        - 步长因子（修正强度系数，建议值0.1~0.3）
%   constellation - 星座图（用于硬判决，如QPSK/16QAM的复数向量）
%
% 输出参数:
%   wout         - 修正后的信号向量

% ================== 1. 初始化与第一阶段处理 ==================
Q = decision(win, constellation);     % 初始硬判决：将信号映射到最近星座点
error = win - Q;                      % 计算原始误差向量（估计值 - 判决值）
error_judge = sign(error);            % 提取误差方向符号（+1表示高估，-1表示低估）
wout = win;                           % 初始化输出信号为输入副本

% 遍历中间符号（跳过首尾无两侧邻居的符号）
for i = 2:length(Q)-1
    % 条件1：检查前向(i-1)和后向(i+1)符号的误差方向是否一致
    if error_judge(i-1) == error_judge(i+1)
        % 条件2：当前符号误差方向是否与相邻符号一致？
        if error_judge(i-1) == error_judge(i)
            % 场景：连续三个符号同向误差（如+++或---）
            % 判定为强脉冲噪声，执行叠加修正
            temp = error(i-1) + error(i+1); 
        else
            % 场景：+-+或-+-
            % 当前符号方向与相邻不同，视为正确判决，不修正
            temp = 0;               
        end
    else
        % 场景：相邻符号误差异向（如++-或--+）
        % 计算相邻误差综合方向（向量和符号）
        judge = sign(error(i-1) + error(i+1));
        % 判断当前符号误差方向是否与综合方向一致
        if error_judge(i) ~= judge
            temp = 0;               % 方向不一致 → 不修正（避免干扰随机噪声）
        else
            temp = error(i-1) + error(i+1); % 方向一致 → 叠加修正
        end
    end
    % 应用修正：alpha/2防止过冲，平衡收敛速度与稳定性
    wout(i) = wout(i) + alpha*temp/2; 
end

% ================== 2. 第二阶段处理（迭代优化）==================
Q = decision(wout, constellation);    % 重新硬判决修正后的信号
error = wout - Q;                    % 计算新误差向量
error_judge = sign(error);           % 更新误差方向符号

% 再次遍历中间符号（复用相同逻辑）
for i = 2:length(Q) - 1
    % 此处逻辑与第一阶段完全一致（实现迭代细化修正）
    if error_judge(i-1) == error_judge(i+1)
        if error_judge(i-1) == error_judge(i)
            temp = error(i-1) + error(i+1);
        else
            temp = 0;
        end
    else
        judge = sign(error(i-1) + error(i+1));
        if error_judge(i) ~= judge
            temp = 0;
        else
            temp = error(i-1) + error(i+1);
        end
    end
    % 二次加权修正（进一步逼近理想值）
    wout(i) = wout(i) + alpha*temp/2;
end

end
