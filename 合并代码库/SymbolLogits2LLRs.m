classdef SymbolLogits2LLRs
    % 将星座符号的对数似然值(logits)转换为比特级的对数似然比(LLRs)或硬判决比特
    % 支持两种计算方法: 'app'(精确计算) 和 'maxlog'(最大对数近似)
    
    properties (SetAccess = private)
        method          % LLR计算方法: 'app' 或 'maxlog'
        hard_out        % 输出硬判决标志
        num_bits_per_symbol % 每个符号的比特数
        c0              % 比特为0的星座点索引
        c1              % 比特为1的星座点索引
        a               % 星座标签(-1,1格式)
    end
    
    methods
        function obj = SymbolLogits2LLRs(method, num_bits_per_symbol, hard_out)
            % 构造函数初始化
            % 参数:
            %   method: 'app' 或 'maxlog'
            %   num_bits_per_symbol: 每个符号的比特数(如QPSK=2, 16QAM=4)
            %   hard_out: true输出硬判决, false输出LLR
            % 验证正确
            obj.method = method;
            obj.hard_out = hard_out;
            obj.num_bits_per_symbol = num_bits_per_symbol;
            num_points = 2^num_bits_per_symbol;
            
            % 生成所有星座点的二进制标签
            a = zeros(num_points, num_bits_per_symbol);
            for i = 0:num_points-1
                bin_str = dec2bin(i, num_bits_per_symbol);
                a(i+1, :) = bin_str - '0';  % 将二进制字符串转为数值数组
            end
            
            % 预计算每个比特位置的星座点索引
            obj.c0 = zeros(num_points/2, num_bits_per_symbol);
            obj.c1 = zeros(num_points/2, num_bits_per_symbol);
            for j = 1:num_bits_per_symbol
                % 找到第j比特为0和1的星座点索引
                obj.c0(:, j) = find(a(:, j) == 0);
                obj.c1(:, j) = find(a(:, j) == 1);
            end
            
            % 生成{-1,1}格式的星座标签
            obj.a = 2*a - 1;  % 0->-1, 1->1
        end
        
        function output = call(obj, logits, prior)
            % 核心计算函数
            % 输入:
            %   logits: [batch_size, num_points] 星座点对数似然值
            %   prior: 可选, [batch_size, num_bits] 或 [1, num_bits] 比特先验信息
            % 输出:
            %   output: LLR或硬判决比特 [batch_size, num_bits_per_symbol]
            
            % 需验证 先验信息 处理部分 和 判决处问题
            [batch_size, num_points] = size(logits);
            num_bits = obj.num_bits_per_symbol;
            
            % ===== 处理先验信息 =====
            if exist('prior', 'var') && ~isempty(prior)
                % 扩展先验维度以匹配logits
                if size(prior, 1) == 1
                    prior = repmat(prior, batch_size, 1);
                end
                
                % 计算星座点先验概率(对数域)
                log_prior_batch = zeros(batch_size, num_points);
                for i = 1:batch_size
                    % 扩展先验到星座点维度
                    p_expanded = repmat(prior(i, :), num_points, 1);
                    
                    % 计算对数先验: sum(log(sigmoid(a * prior))
                    term = obj.a .* p_expanded;
                    log_sigmoid_term = obj.log_sigmoid(term);
                    log_prior_i = sum(log_sigmoid_term, 2)';
                    log_prior_batch(i, :) = log_prior_i;
                end
                
                % 将先验信息添加到logits
                logits = logits + log_prior_batch;
            end
            
            % ===== 计算LLR =====
            llr = zeros(batch_size, num_bits);
            
            % 遍历每个比特位置
            for j = 1:num_bits
                % 获取当前比特的星座点索引
                idx0 = obj.c0(:, j);
                idx1 = obj.c1(:, j);
                
                % 收集比特0/1对应的logits
                exp0 = logits(:, idx0);
                exp1 = logits(:, idx1);
                
                % 根据方法选择约简操作
                if strcmp(obj.method, 'app')
                    % APP方法: log(sum(exp(x)))
                    reduce0 = obj.logsumexp(exp0, 2);
                    reduce1 = obj.logsumexp(exp1, 2);
                else
                    % Max-log方法: max(x)
                    reduce0 = max(exp0, [], 2);
                    reduce1 = max(exp1, [], 2);
                end
                
                % LLR = log(Pr(b=1)) - log(Pr(b=0))
                llr(:, j) = reduce0 - reduce1;
            end
            
            % ===== 输出选择 =====
            if obj.hard_out
                % 硬判决: LLR<0 -> 1, 否则0
                output = llr < 0;
            else
                output = llr;
            end
        end
    end
    
    methods (Static, Access = private)
        function y = log_sigmoid(x)
            % 数值稳定的log(sigmoid)计算
            y = zeros(size(x));
            idx_positive = x > 0;
            
            % x > 0: log(sigmoid(x)) = -log(1 + exp(-x))
            y(idx_positive) = -log(1 + exp(-x(idx_positive)));
            
            % x <= 0: log(sigmoid(x)) = x - log(1 + exp(x))
            y(~idx_positive) = x(~idx_positive) - log(1 + exp(x(~idx_positive)));
        end
        
        function s = logsumexp(x, dim)
            % 数值稳定的log(sum(exp(x)))计算
            % 沿指定维度dim操作
            xmax = max(x, [], dim);
            x = x - xmax;  % 数值稳定处理
            s = log(sum(exp(x), dim)) + xmax;
        end
    end
end
