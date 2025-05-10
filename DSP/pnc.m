function wout = pnc(win, alpha, constellation, stage)

% 通过对接收信号进行迭代修正，抑制脉冲噪声引起的符号失真

% 输入参数:
%   win          - 输入信号向量（含噪声的符号估计值）
%   alpha        - 学习率/步长因子（控制修正幅度，建议范围0.1~0.5）
%   constellation - 星座图（用于硬判决，如QPSK/16QAM的复数向量）
%   stage        - 处理阶段数（1~3，默认值=2）

% 输出参数:
%   wout         - 修正后的输出信号向量

% ================== 1. 参数初始化 ==================
if nargin < 4               % 若未指定stage参数
    stage = 2;              % 设置默认处理阶段为2
end
wout = win;                 % 初始化输出信号为输入信号的副本

% ================== 2. 第一阶段处理 ==================
% 目标：基础噪声平滑，处理明显突发噪声
Q = decision(wout, constellation);  % 对当前信号进行硬判决
error = wout - Q;                   % 计算误差向量（估计值 - 判决值）
error_judge = sign(error);          % 提取误差方向符号（+1/-1）

% 遍历信号中间符号（跳过首尾无双侧邻居的符号）
for i = 2:length(win)-1
    % 策略：若相邻符号误差同向，用平均误差修正当前符号
    if error_judge(i-1) == error_judge(i+1)
        temp = (error(i-1)+error(i+1))/2; % 计算相邻误差均值
    else
        temp = 0;               % 误差方向不一致时不修正
    end
    wout(i) = wout(i) + alpha*temp; % 更新当前符号值
end

% ================== 3. 第二阶段处理（stage>=2时执行）==================
if stage >= 2
    Q = decision(wout, constellation);  % 重新硬判决修正后的信号
    error = wout - Q;                   % 更新误差向量
    error_judge = sign(error);          % 更新误差符号
    
    for i = 2:length(win)-1
        if error_judge(i-1) == error_judge(i+1)
            % 相邻误差同向：增强修正力度（误差直接叠加）
            temp = error(i-1)+error(i+1); 
        else
            % 相邻误差异向：选择较小误差的50%进行保守修正
            if abs(error(i-1)) >= abs(error(i+1))
                temp = error(i+1) / 2;  % 取后侧误差的一半
            else
                temp = error(i-1) / 2;  % 取前侧误差的一半
            end
        end
        wout(i) = wout(i) + alpha*temp; % 二次修正
    end 
end

% ================== 4. 第三阶段处理（stage==3时执行）==================
if stage == 3
    Q = decision(wout, constellation);  % 第三次硬判决
    error = wout - Q;                   % 最终误差计算
    error_judge = sign(error);          % 最终误差符号
    
    for i = 2:length(win)-1
        if error_judge(i-1) == error_judge(i+1)
            % 实验性策略：同向误差直接叠加（可启用注释代码清零）
            temp = error(i-1)+error(i+1); 
            % temp = 0;                 % 保留的备选策略（当前未启用）
        else
            % 维持第二阶段的保守修正策略
            if abs(error(i-1)) >= abs(error(i+1))
                temp = error(i+1) / 2;
            else
                temp = error(i-1) / 2;
            end
        end
        wout(i) = wout(i) + alpha*temp; % 最终修正
    end
end

end
