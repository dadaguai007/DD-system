classdef MaximumLikelihoodDetector
    % MIMO最大似然检测器
    properties
        outputType          % 输出类型 ('bit'/'symbol')
        ConstTransType      % 直接使用星座图 或者 星座图与信道响应乘积
        hardOut             % 硬判决标志
        constellation       % 星座图对象
        vecs                % 所有符号组合 [num_vecs, K] (复数)
        vecsInd             % 符号索引 [num_vecs, K] (1-based)
        cTable              % 符号查找表 [num_indices, K, num_points]
        numStreams          % 发送流数，发射天线数
        numPoints           % 星座点数
        bitsPerSymbol       % 每个符号的比特数
    end
    methods

        % 如何获取正确索引，需要解决
        function obj=MaximumLikelihoodDetector(varargin)
            %初始化类的属性
            if numel(varargin) == 9
                obj.outputType                        = varargin{1} ;% 发射信号的采样率
                obj.hardOut                           = varargin{2} ;% 随机信号的阶数
                obj.constellation                     = varargin{3} ;% prbs码的阶数
                obj.vecs                              = varargin{4} ;% 调制格式
                obj.vecsInd                           = varargin{5} ;% 每符号采样点
                obj.cTable                            = varargin{6} ;% 码元数目
                obj.numStreams                        = varargin{7} ;% 脉冲形式
                obj.bitsPerSymbol                     = varargin{8} ;
                obj.ConstTransType                    = varargin{9} ;
            end
        end

        function [states,statesIndex,vecs_ind_mat,mm]=build_vecs(~,symbols,indexArrary,M,K)
            % 两天线情况,n>3往里填数据即可
            % symbols=[-3,-1,1,3];
            if K==2
                [s1, s2] = ndgrid(symbols,symbols);
                states = [s1(:), s2(:)]; % 生成所有状态组合

                % 索引值
                % indexArrary=[1,2,3,4];
                [s1Index, s2Index] = ndgrid( indexArrary,  indexArrary);
                statesIndex=[s1Index(:),s2Index(:)];
            elseif K==3
                [s1, s2,s3] = ndgrid(symbols,symbols);
                states = [s1(:), s2(:),s3(:)]; % 生成所有状态组合

                % 索引值
                % indexArrary=[1,2,3,4];
                [s1Index, s2Index,s3Index] = ndgrid( indexArrary,  indexArrary);
                statesIndex=[s1Index(:),s2Index(:),s3Index(:)];
            end
            % 天线标称
            %M=4;
            %K=2;
            tx_ind = 1:K;
            tx_ind=repmat(tx_ind,M^K,1);

            % 带有天线目标的索引
            for i=1:K
                vecs_ind=[tx_ind(:,i),statesIndex(:,i)];
                vecs_ind_mat{i}=vecs_ind;
            end

            % 每个天线发送第 s 个星座点的所有符号组合索引
            for j=1:K
                for i=1:length(symbols)
                    ind= find(states(:,j)==symbols(i));
                    mat(:,i)=ind;
                end
                mm{j}={mat};
            end

        end

        function  result= detectMAP(obj, yEq, H, refEq,prior)
            % SISO 系统
            % 都使用最佳算法app计算LLR，不使用maxlog。
            if nargin<5
                M=log2(obj.constellation);
                % 先验信息复制
                prior = 1 / M * ones(size(const));
            end
            % 星座图
            constSymb=obj.constellation;

            % 噪声
            noise = yEq - refEq(1:length(yEq));
            sigma2= var(noise);

            % 星族图对应的比特分布
            num_bits_per_symbol=obj.bitsPerSymbol;
            bitMap= obj.SymbolLogits2LLRs(num_bits_per_symbol);

            % 计算距离度量(接收减去信道与发射信号相作用)
            switch lower(obj.ConstTransType)
                case 'constellation'
                    % 获取unnormalized log-probability（软决策信息【symbol幅度，不涉及比特】）
                    logits=zeros(length(constSymb),length(yEq));
                    for i = 1:length(yEq)
                        logits(:,i) = exp(-abs(yEq(i) - constSymb.').^2 / sigma2) * prior;
                    end

                    % 输出处理
                    if strcmp(obj.outputType, 'bit')
                        % 使用LLR计算输出,MAP方法
                        % 最直接LLR计算
                        LLRs = calcLLR(yEq, sigma2, constSymb, bitMap, prior);
                        LLRs(isinf(LLRs) & LLRs > 0) = 500;    % 将正无穷大替换为500
                        LLRs(isinf(LLRs) & LLRs < 0) = -500;   % 将负无穷大替换为-500
                        % ===== 输出选择 =====
                        if obj.hard_out
                            % 硬判决: LLR<0 -> 1, 否则0
                            result = LLRs < 0;
                        else
                            % 软信息
                            result = LLRs;
                        end
                    else
                        % 直接对符号计算输出
                        if obj.hardOut
                            % 选取每列的最大符号输出
                            [~,constInd] = max(logits);
                            % 输出符号
                            result=constSymb(constInd);
                        else
                            result = logits;
                        end
                    end
                case 'channeltransfer'
                    % 获取unnormalized log-probability（软决策信息【symbol幅度，不涉及比特】）
                    % 此时H大小为 （length(constSymb),length(constSymb)）
                    logits=zeros(length(constSymb),length(yEq));
                    for i = 1:length(yEq)
                        logits(:,i) = exp(-abs(yEq(i) - H*constSymb.').^2 / sigma2) * prior;
                    end
                    % 输出处理
                    if strcmp(obj.outputType, 'bit')
                        % 使用LLR计算输出,MAP方法
                        % 最直接LLR计算
                        LLRs = calcLLR(yEq, sigma2, H*constSymb.', bitMap, prior);
                        LLRs(isinf(LLRs) & LLRs > 0) = 500;    % 将正无穷大替换为500
                        LLRs(isinf(LLRs) & LLRs < 0) = -500;   % 将负无穷大替换为-500
                        % ===== 输出选择 =====
                        if obj.hard_out
                            % 硬判决: LLR<0 -> 1, 否则0
                            result = LLRs < 0;
                        else
                            % 软信息
                            result = LLRs;
                        end
                    else
                        % 直接对符号计算输出
                        if obj.hardOut
                            % 选取每列的最大符号输出
                            [~,constInd] = max(logits);
                            % 输出符号
                            result=constSymb(constInd);
                        else
                            result = logits;
                        end
                    end
            end
        end


        function  result= detectML(obj, yEq, H, refEq)
            % SISO系统
            % 星座图
            constSymb=obj.constellation;

            % 噪声
            noise = yEq - refEq(1:length(yEq));
            sigma2= var(noise);

            % 计算距离度量(接收减去信道与发射信号相作用)
            switch lower(obj.ConstTransType)
                case 'constellation'

                    logits=zeros(length(constSymb),length(yEq));
                    for i = 1:length(yEq)
                        logits(:,i) = (abs(yEq(i) - constSymb.')).^2/sigma2 ;
                    end
                    % 直接对符号计算输出
                    if obj.hardOut
                        % 选取每列的最大符号输出
                        [~,constInd] = min(logits);
                        % 输出符号
                        result=constSymb(constInd);
                    else
                        result = logits;
                    end

                case 'channeltransfer'

                    % 此时H大小为 （length(constSymb),length(constSymb)）
                    logits=zeros(length(constSymb),length(yEq));
                    for i = 1:length(yEq)
                        logits(:,i) = (abs(yEq(i) - H*constSymb.')).^2 / sigma2;
                    end
                    % 直接对符号计算输出
                    if obj.hardOut
                        % 选取每列的最大符号输出
                        [~,constInd] = min(logits);
                        % 输出符号
                        result=constSymb(constInd);
                    else
                        result = logits;
                    end
            end
        end

        % MIMO
        function  result= detectMLMIMO(obj, yEq, H, refEq)
            % 系统
            % 星座图
            constSymb=obj.constellation;
            indexArrary=1:length(constSymb);
            M=log2(constSymb);
            K=2;% 两个天线
            % 噪声
            noise = yEq - refEq(1:length(yEq));
            sigma2= var(noise);
            % 创建Tx的发射组合及索引
            [states,~]=ovj.build_vecs(constSymb,indexArrary,M,K);

            % 计算距离度量(接收减去信道与发射信号相作用)
            switch lower(obj.ConstTransType)
                case 'constellation'
                    % 状态数
                    num_candidates = size(states, 1);
                    % 状态向量
                    Gamma = zeros(1, num_candidates);
                    % 输出向量
                    logits=zeros(length(Gamma),length(yEq));
                    for i = 1:length(yEq)
                        for j = 1:num_candidates
                            states_candidate = states(j,:);
                            Gamma(j)=norm((yEq(i) - states_candidate.')).^2/sigma2 ;
                        end
                        logits(:,i) = Gamma;
                    end
                    % 直接对符号计算输出
                    if obj.hardOut
                        % 选取每列的最小符号组合
                        [~,constInd] = min(logits);
                        % 输出符号
                        result=states(constInd,:);
                    else
                        result = logits;
                    end

                case 'channeltransfer'
                    % 状态数
                    num_candidates = size(states, 1);
                    % 状态向量
                    Gamma = zeros(1, num_candidates);
                    % 输出向量
                    logits=zeros(length(Gamma),length(yEq));
                    for i = 1:length(yEq)
                        for j = 1:num_candidates
                            states_candidate = states(j,:);
                            Gamma(j)=norm((yEq(i) - H*states_candidate.')).^2/sigma2 ;
                        end
                        logits(:,i) = Gamma;
                    end
                    % 直接对符号计算输出
                    if obj.hardOut
                        % 选取每列的最小符号组合输出
                        [~,constInd] = min(logits);
                        % 输出符号
                        result=states(constInd,:);
                    else
                        result = logits;
                    end
            end
        end

        % MIMO
        function  result= detectMAPMIMO(obj, yEq, H, refEq,prior)
            %  MIMO系统
            % 都使用最佳算法app计算LLR，不使用maxlog。
            if nargin<5
                M=log2(obj.constellation);
                % 先验信息复制
                prior = 1 / M * ones(size(const));
            end
            % 星座图
            constSymb=obj.constellation;
            % 索引
            indexArrary=1:length(constSymb);
            M=log2(constSymb);
            K=2;% 两个天线
            % 创建Tx的发射组合及索引
            [states,statesIndex]=ovj.build_vecs(constSymb,indexArrary,M,K);
            % 噪声
            noise = yEq - refEq(1:length(yEq));
            sigma2= var(noise);

            % 星族图对应的比特分布 * K
            num_bits_per_symbol=obj.bitsPerSymbol;
            bitMap= obj.SymbolLogits2LLRsMIMo(num_bits_per_symbol,K);
            % 计算距离度量(接收减去信道与发射信号相作用)
            switch lower(obj.ConstTransType)
                case 'constellation'
                    % 状态数
                    num_candidates = size(states, 1);
                    % 状态向量
                    Gamma = zeros(1,num_candidates);
                    % 输出向量
                    logits=zeros(length(Gamma),length(yEq));
                    for i = 1:length(yEq)
                        for j = 1:num_candidates
                            states_candidate = states(j,:);
                            Gamma(j)=sum(exp(-abs(yEq(i) - states_candidate.').^2 / sigma2) * prior);
                        end
                        logits(:,i) = Gamma;
                    end
                    % 输出处理
                    if strcmp(obj.outputType, 'bit')
                        % 使用LLR计算输出,MAP方法
                        % 最直接LLR计算
                        LLRs = calcLLR_MIMO(yEq, sigma2, K,states,statesIndex, bitMap, px);
                        LLRs(isinf(LLRs) & LLRs > 0) = 500;    % 将正无穷大替换为500
                        LLRs(isinf(LLRs) & LLRs < 0) = -500;   % 将负无穷大替换为-500
                        % ===== 输出选择 =====
                        if obj.hard_out
                            % 硬判决: LLR<0 -> 1, 否则0
                            result = LLRs < 0;
                        else
                            % 软信息
                            result = LLRs;
                        end
                    else
                        % 直接对符号计算输出
                        if obj.hardOut
                            % 选取每列的最大符号输出
                            [~,constInd] = max(logits);
                            % 输出符号
                            result=states(constInd,:);
                        else
                            result = logits;
                        end
                    end
                case 'channeltransfer'
                    % 状态数
                    num_candidates = size(states, 1);

                    states_H=zeros(size(states));
                    % 状态向量
                    Gamma = zeros(1,num_candidates);
                    % 输出向量
                    logits=zeros(length(Gamma),length(yEq));
                    for i = 1:length(yEq)
                        for j = 1:num_candidates
                            states_candidate = states(j,:);
                            Gamma(j)=sum(exp(-abs(yEq(i) -H*states_candidate.').^2 / sigma2) * prior);
                            % 信道乘上发射信号
                            states_H(:,j)=H*states_candidate.';
                        end
                        logits(:,i) = Gamma;
                    end
                    % 输出处理
                    if strcmp(obj.outputType, 'bit')
                        % 使用LLR计算输出,MAP方法
                        % 最直接LLR计算
                        LLRs = calcLLR_MIMO(yEq, sigma2, K,states_H,statesIndex, bitMap, px);
                        LLRs(isinf(LLRs) & LLRs > 0) = 500;    % 将正无穷大替换为500
                        LLRs(isinf(LLRs) & LLRs < 0) = -500;   % 将负无穷大替换为-500
                        % ===== 输出选择 =====
                        if obj.hard_out
                            % 硬判决: LLR<0 -> 1, 否则0
                            result = LLRs < 0;
                        else
                            % 软信息
                            result = LLRs;
                        end
                    else
                        % 直接对符号计算输出
                        if obj.hardOut
                            % 选取每列的最大符号输出
                            [~,constInd] = max(logits);
                            % 输出符号
                            result=states(constInd,:);
                        else
                            result = logits;
                        end
                    end
            end
        end


        function  bitMap= SymbolLogits2LLRs(~,num_bits_per_symbol)
            % 构造函数初始化
            % 参数:
            %   method: 'app'
            %   num_bits_per_symbol: 每个符号的比特数(如QPSK=2, 16QAM=4)
            %   hard_out: true输出硬判决, false输出LLR
            M = 2^num_bits_per_symbol;

            % 生成所有星座点的二进制标签
            bitMap = zeros(M, num_bits_per_symbol);
            for i = 0:M-1
                bin_str = dec2bin(i, num_bits_per_symbol);
                bitMap(i+1, :) = bin_str - '0';  % 将二进制字符串转为数值数组
            end
        end
        %MIMO的索引构造
        function  bitMap= SymbolLogits2LLRsMIMo(~,num_bits_per_symbol,K)
            % 构造函数初始化
            % 参数:
            %   method: 'app'
            %   num_bits_per_symbol: 每个符号的比特数(如QPSK=2, 16QAM=4)
            %   hard_out: true输出硬判决, false输出LLR
            M = 2^num_bits_per_symbol;
            indexArrary=1:M;
            symbols=real(pammod(indexArrary-1,M));
            % 生成状态矩阵
            [states,statesIndex]=obj.build_vecs(symbols,indexArrary,M,K);

            % 生成所有星座点的二进制标签
            bitMap = zeros(length(states), num_bits_per_symbol*K);
            for index =1:K
                for j = 1:length(statesIndex)
                    i=statesIndex(j,index)-1;
                    % bin解码
                    bin_str = dec2bin(i, num_bits_per_symbol);
                    bitMap(j, num_bits_per_symbol*(index-1)+1:num_bits_per_symbol*index) = bin_str - '0';  % 将二进制字符串转为数值数组
                end
            end
        end
    end

end
