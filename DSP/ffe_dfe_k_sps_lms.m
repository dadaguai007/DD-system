function [wout, err, w] = ffe_dfe_k_sps_lms(win, ML, step, ref, sps, iter)
% 分数间隔采样FFE-DFE联合均衡器 - LMS自适应实现
% 输入:
%   win   : 输入信号向量(接收信号)
%   ML    : 二维向量 [FFE抽头数, DFE抽头数] (符号长度)
%   step  : 二维向量 [FFE步长, DFE步长]
%   ref   : 训练序列(已知参考信号)
%   sps   : 每符号采样数(支持1.0,1.2,1.4,1.6,1.8,2.0)
%   iter  : 训练迭代次数
% 输出:
%   wout  : 均衡后的输出信号
%   err   : 训练误差序列
%   w     : 最终FFE权重矩阵(每列对应一个相位)

    % === 1. 采样模式定义 ===
    tb = [1, 1, 1, 1, 1;   % sps=1.0
          2, 1, 1, 1, 1;   % sps=1.2
          2, 1, 2, 1, 1;   % sps=1.4
          2, 1, 2, 1, 2;   % sps=1.6
          2, 2, 2, 2, 1;   % sps=1.8
          2, 2, 2, 2, 2];  % sps=2.0
    
    % === 2. 选择采样模式 ===
    switch sps
        case 1
            tbref = tb(1, :);
        case 1.2
            tbref = tb(2, :);
        case 1.4
            tbref = tb(3, :);
        case 1.6
            tbref = tb(4, :);
        case 1.8
            tbref = tb(5, :);
        case 2
            tbref = tb(6, :);
        otherwise
            error('仅支持 1, 1.2, 1.4, 1.6, 1.8, 2sps');
    end
    
    % === 3. 初始化参数 ===
    numtaps = round(ML(1) * sps);    % FFE实际抽头数
    w = zeros(numtaps, 5);           % FFE权重矩阵(5相位)
    w(round(numtaps/2), :) = 1;      % FFE中心抽头初始化为1
    
    trainlen = length(ref);          % 训练序列长度
    err = zeros(iter*trainlen, 1);   % 误差向量初始化
    
    % 计算输出长度
    enum = floor(length(win) / sps); % 输入信号包含的符号数
    L = enum - round(max(ML)*sps/2) + 1; % 有效输出长度
    wout = zeros(L, 1);              % 输出向量初始化
    
    % DFE初始化
    wb = zeros(ML(2), 1);            % DFE权重向量
    wb(round(ML(2)/2)) = 1;          % DFE中心抽头初始化
    
    constellation = unique(ref);     % 获取调制星座点(用于判决)
    
    % === 4. FFE-DFE联合训练 ===
    block = zeros(numtaps, 1);       % FFE输入缓冲区
    temp = zeros(ML(2), 1);          % DFE输入缓冲区(存储判决符号)
    for i = 1:iter                   % 多轮迭代训练
        number = 1;                  % 输入信号索引
        for j = 1:trainlen           % 遍历训练序列
            % -- 相位计算(1-5循环) --
            index = mod(j, 5);
            if index == 0
                index = 5;
            end
            
            % -- 动态采样更新FFE缓冲区 --
            if tbref(index) == 1     % 单采样
                block = [win(number); block(1:end-1)];
                number = number + 1;
            else                     % 双采样
                block = [win(number+1); win(number); block(1:end-2)];
                number = number + 2;
            end
            
            % -- 联合均衡计算 --
            yf = w(:, index).' * block;  % FFE输出
            yb = wb.' * temp;            % DFE输出
            y = yf + yb;                 % 联合输出
            
            % -- 误差计算 --
            err((i-1)*trainlen+j) = ref(j) - y;
            
            % -- LMS权重更新 --
            w(:, index) = w(:, index) + step(1)*err((i-1)*trainlen+j)*conj(block); % FFE更新
            wb = wb + step(2)*err((i-1)*trainlen+j)*conj(temp);                   % DFE更新
            
            % -- 更新DFE缓冲区(使用已知参考信号) --
            temp = [ref(j); temp(1:end-1)]; % 滑动存储参考符号
        end
    end
    
    % === 5. 均衡输出阶段 ===
    block = zeros(numtaps, 1);       % 重置FFE缓冲区
    temp = zeros(ML(2), 1);          % 重置DFE缓冲区
    number = 1;                      % 重置输入索引
    for i = 1:L                      % 生成最终输出
        % -- 相位同步 --
        index = mod(i, 5);
        if index == 0
            index = 5;
        end
        
        % -- 动态采样更新FFE缓冲区 --
        if tbref(index) == 1
            block = [win(number); block(1:end-1)];
            number = number + 1;
        else
            block = [win(number+1); win(number); block(1:end-2)];
            number = number + 2;
        end
        
        % -- 联合均衡计算 --
        yf = w(:, index).' * block;  % FFE输出
        yb = wb.' * temp;            % DFE输出
        wout(i) = yf + yb;           % 最终输出
        
        % -- 判决反馈更新DFE --
        dec_symbol = decision(wout(i), constellation); % 星座判决
        temp = [dec_symbol; temp(1:end-1)]; % 用判决符号更新DFE缓冲区
    end
end