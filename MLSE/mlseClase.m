classdef mlseClase
    % MLSE均衡器类
    properties
        channelLength   % 信道冲激响应长度（抽头数）
        blockLength     % 数据块长度（比特数）
        states_viterbi  % 维特比算法状态矩阵 [2^channelLength × channelLength]
        states_mlse     % 全搜索算法状态矩阵 [2^blockLength × blockLength]
    end
    methods
        function this = mlseClase(channelLength, blockLength)
            % 初始化信道参数并生成状态矩阵
            this.channelLength = channelLength;  % 信道记忆长度（如3抽头）
            this.blockLength = blockLength;      % 传输块长度（如100比特）
            this.states_mlse = this.states_calc(blockLength);    % 全搜索状态（所有可能的发送序列）
            this.states_viterbi = this.states_calc(channelLength);% 维特比状态（信道状态组合）
        end

        function states = states_calc(varargin)
            % 生成所有可能的二进制状态序列并调制
            % BPSK调制后的状态矩阵;可用于全状态MLSE算法，也可用于维特比算法！
            len = varargin{2};  % 状态序列长度
            states_char = dec2bin(0:2^len - 1);         % 生成二进制序列（如'000','001'...）
            states_char_reshaped = reshape(states_char, 1, [])';
            states = reshape(str2num(states_char_reshaped), [], len); % 转换为数值矩阵
            states = 2*states - 1;  % BPSK调制（0→-1，1→+1）
            % len=3 时生成 8×3 矩阵，包含所有可能的 [-1, +1] 组合。
        end

        function equalized = run(this, receivedSignal, estimatedCir)
            assert(iscolumn(receivedSignal), "Input signal must be column vector");
            assert(iscolumn(estimatedCir), "Estimated CIR must be column vector");
            assert(length(receivedSignal) == this.channelLength + this.blockLength - 1, "Input signal has wrong size");
            assert(length(estimatedCir) == this.channelLength, "Estimated CIR has wrong size");
            % 全搜索计算：
            % 穷举所有 2^blockLength 种可能的发送序列（BPSK调制后为±1组合）。
            % 对每个候选序列与估计CIR卷积，生成理论接收信号。
            % 计算候选信号与实际接收信号的L1误差（绝对值之和），选择最小误差对应的序列。

            equalized = mlse_full(this, receivedSignal, estimatedCir);
            equalized_mex = MLSEmexFunc('run', this.prepare_data_for_mex(receivedSignal.'), this.prepare_data_for_mex(estimatedCir.'));

            if ~isequal(equalized, double(equalized_mex))
                error("MATLAB and C++ realisations show different results");
            else
                disp("MATLAB and C++ realisations are the same");
            end

            %             equalized = viterbi(this, receivedSignal, estimatedCir);
        end


        function data_out = prepare_data_for_mex(this, data_in)
            % 数据格式转换：复数→实虚分离+单精度
            data_out = [real(data_in); imag(data_in)];
            data_out = reshape(data_out, 1, []);
            data_out = single(data_out);% C++接口要求
        end

        function equalized = mlse_full(this, receivedSignal, estimatedCir)
           % 核心思想：暴力遍历所有可能的发送序列，选择与接收信号最匹配的序列。
           % 步骤：
           % 候选信号生成：对每个可能的发送序列（states_mlse 的行）与CIR进行卷积，模拟信道失真。
           % 误差计算：计算每个候选信号与接收信号的L1范数误差。
           % 最优选择：选择误差最小的序列作为均衡结果。

            candidates = zeros(size(this.states_mlse, 1), this.channelLength + this.blockLength - 1);
            for i = 1:size(this.states_mlse, 1)
                % 卷积模拟信道效应
                candidates(i, :) = conv(this.states_mlse(i, :), estimatedCir);
            end
            err = sum(abs(candidates - receivedSignal.'), 2);% L1范数误差
            [~, idx] = min(err); % 选择最小误差序列
            equalized = this.states_mlse(idx, :);% 返回最优解
        end
    end
end