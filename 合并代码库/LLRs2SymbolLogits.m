classdef LLRs2SymbolLogits
    % 将比特LLRs转换为星座符号的logits(未归一化对数概率)或硬判决符号索引
    %
    % 参数:
    %   num_bits_per_symbol: 每个符号的比特数 (如QPSK=2, 16QAM=4)
    %   hard_out: true输出硬判决符号索引, false输出符号logits
    %
    % 输入:
    %   llrs: [..., n, num_bits_per_symbol] 每个比特的对数似然比
    %
    % 输出:
    %   logits: [..., n, num_points] 星座符号的logits (当hard_out=false)
    %   或
    %   indices: [..., n] 硬判决符号索引 (当hard_out=true)
    
    properties (SetAccess = private)
        num_bits_per_symbol  % 每个符号的比特数
        hard_out             % 输出模式标志
        a                    % 星座标签矩阵 (-1/1格式)
    end
    
    methods
        function obj = LLRs2SymbolLogits(num_bits_per_symbol, hard_out)
            % 构造函数初始化
            obj.num_bits_per_symbol = num_bits_per_symbol;
            obj.hard_out = hard_out;
            
            % 计算星座点数
            num_points = 2^num_bits_per_symbol;
            
            % 生成所有星座点的二进制标签矩阵
            a = zeros(num_points, num_bits_per_symbol);
            for i = 0:num_points-1
                % 将符号索引转为固定长度的二进制表示
                bin_str = dec2bin(i, num_bits_per_symbol);
                a(i+1, :) = double(bin_str) - '0'; % 将字符串转为数值数组
            end
            
            % 将二进制标签转换为{-1,1}格式 (0->-1, 1->1)
            obj.a = 2*a - 1;
        end
        
        function output = call(obj, llrs)
            % 核心计算函数: 将LLRs转换为符号logits或硬判决索引
            %
            % 输入:
            %   llrs: 比特对数似然比矩阵 [..., n, num_bits_per_symbol]
            %
            % 输出:
            %   output: 符号logits或硬判决索引
            
            % 获取输入维度信息
            orig_size = size(llrs);
            K = obj.num_bits_per_symbol;  % 每个符号比特数
            M = size(obj.a, 1);           % 星座点数
            
            % 将输入重塑为二维矩阵 [N, K]
            % 其中N = 前面所有维度的乘积
            N = prod(orig_size(1:end-1)); 
            llrs_2d = reshape(llrs, [N, K]);
            
            % 扩展维度用于广播计算
            llrs_3d = reshape(llrs_2d, [N, 1, K]);  % [N, 1, K]
            a_3d = reshape(obj.a, [1, M, K]);       % [1, M, K]
            
            % 核心计算: 星座符号的logit = Σ log(sigmoid(LLR(k) * 星座标签(k))
            product = a_3d .* llrs_3d;        % [N, M, K]
            log_sig = obj.log_sigmoid(product); % 计算log(sigmoid)
            logits_2d = sum(log_sig, 3);       % 沿比特维度求和 [N, M]
            
            % 重塑为原始维度结构 (前面维度 + 星座点数)
            new_size = [orig_size(1:end-1), M];
            logits = reshape(logits_2d, new_size);
            
            % 输出选择
            if obj.hard_out
                % 硬判决: 沿最后一个维度取最大值索引
                [~, indices] = max(logits, [], ndims(logits));
                output = indices - 1; % 转换为0-based索引
            else
                output = logits; % 软输出: 所有符号的logits
            end
        end
    end
    
    methods (Access = private)
        function y = log_sigmoid(~, x)
            % 数值稳定的log(sigmoid)计算
            %
            % 对于x>0: log(sigmoid(x)) = -log(1 + exp(-x))
            % 对于x<=0: log(sigmoid(x)) = x - log(1 + exp(x))
            
            y = zeros(size(x)); % 初始化输出矩阵
            
            % 找出正数和负数部分
            idx_positive = x > 0;
            idx_non_positive = ~idx_positive;
            
            % 正数部分计算
            if any(idx_positive(:))
                y(idx_positive) = -log(1 + exp(-x(idx_positive)));
            end
            
            % 非正数部分计算
            if any(idx_non_positive(:))
                % 使用log1p提高数值稳定性 (log(1+exp(x)))
                y(idx_non_positive) = x(idx_non_positive) - log1p(exp(x(idx_non_positive)));
            end
        end
    end
end
