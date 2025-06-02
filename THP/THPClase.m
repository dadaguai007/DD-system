classdef THPClase
    % THP均衡器类
    properties
        channelLength   % 信道冲激响应长度（抽头数）
        blockLength     % 数据块长度（比特数）
    end

    methods
        % 此函数只给出频域用Z变换表示的例子
        % 实践中，应该测量信道响应，在发射端做预均衡

        % H_inv_f = 1 ./ H_f; % 信道逆频响
        % H_inv_f(abs(H_f) < eps) = 0; % 避免除以零
        % X_pre_f = X_f .* H_inv_f; % 预均衡信号频域

        function obj=THPClase(channelLength, blockLength)
            % 初始化信道参数并生成状态矩阵
            obj.channelLength = channelLength;  % 信道记忆长度（如3抽头）
            obj.blockLength = blockLength;      % 传输块长度（如100比特）
        end

        function z_x=Transz(~,x)
            % 时域转换为z域信号
            syms z;                     % 声明符号变量z（用于Z域分析）
            z_x = 0;                    % 初始化Z变换结果
            for i = 1:length(x)
                z_x = z_x + x(i)*z^(1-i); % 逐个元素累加计算Z变换：X(z) = Σx(n)z⁻ⁿ
            end
        end

        function output=zTrans(obj,z_output,x)
            % 从z变换转换为时域信号
            % z_output为表达式，x为初始时域信号
            o = iztrans(z_output);      % 逆Z变换得到时域表达式
            output = obj.delta2sequence(o, length(x)); % 转换为离散序列
        end

        function mod = modulo(~,sequence, N)
            % 模运算函数：将信号限制在[-N/2, N/2)范围内
            for i = 1:length(sequence)
                % 上界处理（大于N/2时循环减N）
                while sequence(i) > N/2
                    sequence(i) = sequence(i) - N;
                end
                % 下界处理（小于-N/2时循环加N）
                while sequence(i) < -N/2
                    sequence(i) = sequence(i) + N;
                end
            end
            mod = sequence;
        end

        function s = delta2sequence(~,expression, len)
            % 将逆Z变换结果转换为离散序列
            % 输入：expression-符号表达式，len-序列长度
            % 输出：s-长度为len的离散序列
            s = zeros(len, 1);
            for i = 0:len-1
                syms n; % 声明符号变量n
                % 代入n=i计算表达式值
                s(i+1) = double(subs(expression, n, i));
            end
        end

        % 时域抽头反馈THP
        function [pre_equalised_with_mod,xInput]=preTHP(obj,xn,w,N,sps)
            % THP
            if nargin < 5
                sps=1;
            end
            n = length(xn);
            % 抽头数
            taps=length(w);
            % 反馈向量
            x1=zeros(taps, 1);
            re_xn=repmat(xn,1,2);
            % THP预均衡
            pre_equalised_with_mod = zeros(n, 1);
            xInput = zeros(n, 1);
            for i=1:n
                yn=x1.'*w;
                % 反馈
                xInput(i)=xn(i)+yn;
                % 模运算控制信号幅度
                pre_equalised_with_mod (i)= obj.modulo(xInput(i), N); % 将信号限制在[-N/2, N/2)
                % 取taps数的
                x1 = cat(1,x1(sps+1:end),xn(sps*(i-1)+1:1:sps*i));
            end

        end
    end
end
