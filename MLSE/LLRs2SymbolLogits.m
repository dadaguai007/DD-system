classdef LLRs2SymbolLogits < handle
    properties (SetAccess = private)
        hardOut         % 硬判决标志
        numBitsPerSym   % 每个符号的比特数
        a               % 符号标签 [-1,1] 矩阵 [num_points, num_bits]
    end
    
    methods
        function obj = LLRs2SymbolLogits(numBitsPerSym, hardOut)
            obj.numBitsPerSym = numBitsPerSym;
            obj.hardOut = hardOut;
            numPoints = 2^numBitsPerSym;
            
            % 生成符号标签矩阵
            a = zeros(numPoints, numBitsPerSym);
            for i = 0:numPoints-1
                a(i+1,:) = dec2bin(i, numBitsPerSym) - '0';
            end
            obj.a = 2*a - 1; % 转换为[-1,1]
        end
        
        function logits = call(obj, llrs)
            numBits = obj.numBitsPerSym;
            numPoints = 2^numBits;
            n = size(llrs, 1); % 符号数
            
            % 计算符号对数概率 [n, numPoints]
            logits = sum(log(sigmoid(llrs * obj.a')), 2); % 矩阵乘法直接求和
            
            % 硬判决
            if obj.hardOut
                [~, logits] = max(logits, [], 2); % 返回1-based索引
            end
        end
    end
end