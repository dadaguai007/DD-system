function [wout, err, w] = ffe_k_sps_lms(win, ML, step, ref, sps, iter)
% 分数间隔采样FFE均衡器 - LMS自适应实现
% 输入:
%   win   : 输入信号向量
%   ML    : 均衡器抽头数(符号长度)
%   step  : LMS步长
%   ref   : 训练序列(参考信号)
%   sps   : 每符号采样数(支持1.0,1.2,1.4,1.6,1.8,2.0)
%   iter  : 训练迭代次数
% 输出:
%   wout  : 均衡后的输出信号
%   err   : 训练误差序列
%   w     : 最终权重矩阵(每列对应一个相位)

    % === 1. 采样模式定义 ===
    tb = [1, 1, 1, 1, 1;   % sps=1.0: 5符号取5样本
          2, 1, 1, 1, 1;   % sps=1.2: 5符号取6样本
          2, 1, 2, 1, 1;   % sps=1.4: 5符号取7样本
          2, 1, 2, 1, 2;   % sps=1.6: 5符号取8样本
          2, 2, 2, 2, 1;   % sps=1.8: 5符号取9样本
          2, 2, 2, 2, 2];  % sps=2.0: 5符号取10样本
    
    % === 2. 选择采样模式 ===
    switch sps
        case 1
            tbref = tb(1, :);  % 选择sps=1.0模式
        case 1.2
            tbref = tb(2, :);  % 选择sps=1.2模式
        case 1.4
            tbref = tb(3, :);  % 选择sps=1.4模式
        case 1.6
            tbref = tb(4, :);  % 选择sps=1.6模式
        case 1.8
            tbref = tb(5, :);  % 选择sps=1.8模式
        case 2
            tbref = tb(6, :);  % 选择sps=2.0模式
        otherwise
            error('仅支持 1, 1.2, 1.4, 1.6, 1.8, 2sps');
    end
    
    % === 3. 初始化参数 ===
    numtaps = round(ML * sps);       % 计算实际抽头数
    w = zeros(numtaps, 5);           % 初始化5相位权重矩阵
    w(round(numtaps/2), :) = 1;      % 中心抽头初始化为1(单位脉冲)
    trainlen = length(ref);          % 训练序列长度
    err = zeros(iter*trainlen, 1);   % 误差向量初始化
    
    % === 4. 计算输出长度 ===
    enum = floor(length(win) / sps); % 输入信号包含的符号数
    L = enum - round(max(ML)*sps/2) + 1; % 有效输出长度(考虑延迟)
    wout = zeros(L, 1);              % 输出向量初始化
    
    % === 5. LMS训练阶段 ===
    block = zeros(numtaps, 1);       % 初始化FFE输入缓冲区
    for i = 1:iter                   % 多轮迭代训练
        number = 1;                  % 输入信号索引复位
        for j = 1:trainlen           % 遍历训练序列
            % -- 相位计算(1-5循环) --
            index = mod(j, 5);       % 计算当前相位索引
            if index == 0
                index = 5;           % 相位5处理
            end
            
            % -- 动态采样控制 --
            if tbref(index) == 1     % 单采样模式
                % 插入新样本并移位: [新样本 | 旧样本(1:end-1)]
                block = [win(number); block(1:end-1)]; 
                number = number + 1; % 移动单个样本指针
            else                     % 双采样模式
                % 插入两个新样本: [新样本2 | 新样本1 | 旧样本(1:end-2)]
                block = [win(number+1); win(number); block(1:end-2)]; 
                number = number + 2; % 移动双样本指针
            end
            
            % -- LMS计算与更新 --
            y = w(:, index).' * block; % 均衡输出(点积)
            err((i-1)*trainlen+j) = ref(j) - y; % 误差计算
            % 权重更新: w_new = w_old + μ * e * x*
            w(:, index) = w(:, index) + step*err((i-1)*trainlen+j)*conj(block);
        end
    end
    
    % === 6. 均衡输出阶段 ===
    block = zeros(numtaps, 1);       % 重置缓冲区
    number = 1;                      % 重置输入索引
    for i = 1:L                      % 生成最终输出
        % -- 相位同步(与训练相同) --
        index = mod(i, 5);
        if index == 0
            index = 5;
        end
        
        % -- 动态采样(与训练相同逻辑) --
        if tbref(index) == 1
            block = [win(number); block(1:end-1)];
            number = number + 1;
        else
            block = [win(number+1); win(number); block(1:end-2)];
            number = number + 2;
        end
        
        % -- 计算均衡输出 --
        wout(i) = w(:, index).' * block; % 使用训练好的权重
    end
end
