classdef MaximumLikelihoodDetector
    % MIMO最大似然检测器
    properties
        outputType          % 输出类型 ('bit'/'symbol')
        demappingMethod     % 检测方法 ('app'/'maxlog')
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
            if numel(varargin) == 8
                obj.outputType                           = varargin{1} ;% 发射信号的采样率
                obj.demappingMethod                          = varargin{2} ;% 发射信号的波特率
                obj.hardOut                       = varargin{3} ;% 随机信号的阶数
                obj.constellation                    = varargin{4} ;% prbs码的阶数
                obj.vecs                            = varargin{5} ;% 调制格式
                obj.vecsInd                          = varargin{6} ;% 每符号采样点
                obj.cTable                         = varargin{7} ;% 码元数目
                obj.numStreams                        = varargin{8} ;% 脉冲形式

            end
            function [states,statesIndex,vecs_ind_mat,mm]=build_vecs(~,symbols,indexArrary,M,K)
                % 两天线情况,n>3往里填数据即可
                % symbols=[-3,-1,1,3];
                [s1, s2] = ndgrid( symbols,  symbols);
                states = [s2(:), s1(:)]; % 生成所有状态组合

                % 索引值
                % indexArrary=[1,2,3,4];
                [s1Index, s2Index] = ndgrid( indexArrary,  indexArrary);
                statesIndex=[s2Index(:),s1Index(:)];

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

            function distances = detect(obj, y, H, prior)

                % 将比特级先验 LLR 转换为符号级对数概率后，可与检测器计算的符号似然（logits）直接相加，实现 MAP
                % 计算距离度量(接收减去信道与发射信号相作用)
                % x1* H + x2*H = y
                for i=1:length(i)
                    diff(i,:) = y(i,:)- sum(H * states);
                    exponents(i) = -sum(conj(diff(i,:)) .* diff(i,:));
                end

                % 融合先验信息（一般先验概率为0）
                if ~isempty(prior)
                    if strcmp(obj.outputType, 'bit')
                        prior = obj.llrToSymbolLogits(prior);
                    end
                    % 扩展先验维度并查找
                    % 需要完善
                    priorVec = prior(vecs_ind_mat,:);
                    exponents = exponents + sum(priorVec);
                end

                % 第k个天线，发送第 s 个星座点的所有符号组合索引。
                for j=1:K
                    k=cell2mat(mat(j));
                end
                for i=1:length(symbols)
                    exponents(k(:,i))
                end
                % 两种方法计算度量值
                if strcmp(demappingMethod,'app')
                    logits = logsumexp(exponents);
                elseif strcmp(demappingMethod,'maxlog')
                    logits = max(exponents);
                end
                % 输出处理
                if strcmp(obj.outputType, 'bit')
                    % 使用LLR计算输出
                    result = logits2llrFunc(logits);
                    %
                    %LLRs = calcLLR(y, sigma2, constSymb, bitMap, px);
                else
                    % 直接对符号计算输出
                    if obj.hardOut
                        result = max(logits);
                    else
                        result = logits;
                    end
                end


            end


        end

    end
end