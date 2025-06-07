% ======================= NLMS类 (归一化最小均方算法) =======================
classdef NLMS
    properties
        step_size  % 学习步长 (μ)
        order     % 滤波器阶数 (抽头数)
    end
    
    methods
        function obj = NLMS(step_size, order)
            % 构造函数: 初始化步长和滤波器阶数
            obj.step_size = step_size;
            obj.order = order;
        end
        
        function [f, w_trk, e] = learn(obj, u, d)
            % NLMS自适应滤波算法
            % 输入:
            %   u - 输入信号向量
            %   d - 期望输出信号向量
            % 输出:
            %   f - 滤波输出信号
            %   w_trk - 权重演变轨迹(每行代表一个时间点的权重)
            %   e - 误差信号
            
            N = length(d);          % 信号长度
            f = zeros(1, N);        % 初始化滤波输出
            e = zeros(1, N);        % 初始化误差信号
            w = rand(1, obj.order); % 随机初始化权重向量
            w_trk = zeros(N - obj.order, obj.order); % 权重轨迹矩阵
            
            for i = 1:(N - obj.order - 1)
                % 提取当前输入向量
                u_segment = u(i:i + obj.order - 1);
                
                % 计算滤波器输出 (内积)
                y = dot(w, u_segment);
                
                % 计算误差信号
                e(i) = d(i + obj.order) - y;
                
                % 计算输入向量的能量 (防止除零)
                norm_val = norm(u_segment)^2 + 0.1;
                
                % NLMS权重更新 (归一化步长)
                w = w + (obj.step_size * e(i) * u_segment) / norm_val;
                
                % 存储结果
                f(i) = y;
                w_trk(i, :) = w;
            end
        end
    end
end

% ======================= 高斯核函数基类 =======================
classdef GaussianKernel
    methods
        function k = kernel(obj, a, b, sigma)
            % 高斯核函数计算
            % 输入:
            %   a, b - 输入向量
            %   sigma - 核宽度参数
            % 输出:
            %   k - 高斯核函数值
            numer = norm(a - b)^2;
            denom = 2 * sigma^2;
            k = exp(-numer / denom);
        end
    end
end

% ======================= KLMS类 (核最小均方算法) =======================
classdef KLMS < GaussianKernel
    properties
        step_size  % 学习步长
        sigma      % 高斯核宽度参数
        order      % 输入向量长度
    end
    
    methods
        function obj = KLMS(step_size, sigma, order)
            % 构造函数: 初始化算法参数
            obj.step_size = step_size;
            obj.sigma = sigma;
            obj.order = order;
        end
        
        function [estimates, coefficients, errors] = learn(obj, input_signal, desired_signal)
            % KLMS非线性自适应滤波算法
            % 输入:
            %   input_signal - 输入信号
            %   desired_signal - 期望信号
            % 输出:
            %   estimates - 估计输出
            %   coefficients - 系数向量
            %   errors - 误差信号
            
            N = length(input_signal);
            estimates = zeros(1, N);       % 初始化估计输出
            coefficients = zeros(1, N);    % 初始化系数向量
            errors = zeros(1, N);          % 初始化误差信号
            
            % 初始化第一个系数
            coefficients(1) = obj.step_size * desired_signal(1);
            
            for i = 2:(N - 2*obj.order)
                % 提取当前输入向量
                current_vec = input_signal(i:i + obj.order - 1);
                
                % 计算核展开输出
                for j = 1:(i-1)
                    % 提取历史输入向量
                    hist_vec = input_signal(j:j + obj.order - 1);
                    
                    % 计算核函数值并加权累加
                    k_val = obj.kernel(current_vec, hist_vec, obj.sigma);
                    estimates(i) = estimates(i) + coefficients(j) * k_val;
                end
                
                % 计算误差信号
                errors(i) = desired_signal(i + obj.order) - estimates(i);
                
                % 更新系数
                coefficients(i) = obj.step_size * errors(i);
            end
        end
    end
end

% ======================= QKLMS类 (量化核最小均方算法) =======================
classdef QKLMS < GaussianKernel
    properties
        step_size  % 学习步长
        sigma      % 高斯核宽度参数
        order      % 输入向量长度
        epsilon    % 量化阈值 (控制中心点数量)
    end
    
    methods
        function obj = QKLMS(step_size, sigma, order, epsilon)
            % 构造函数: 初始化算法参数
            obj.step_size = step_size;
            obj.sigma = sigma;
            obj.order = order;
            obj.epsilon = epsilon;
        end
        
        function [estimates, coefficients, errors, centers_size] = learn(obj, input_signal, desired_signal)
            % QKLMS算法 (带矢量量化的KLMS)
            % 输入:
            %   input_signal - 输入信号
            %   desired_signal - 期望信号
            % 输出:
            %   estimates - 估计输出
            %   coefficients - 系数向量
            %   errors - 误差信号
            %   centers_size - 量化中心数量演变
            
            N = length(input_signal);
            estimates = zeros(1, N);        % 初始化估计输出
            coefficients = zeros(1, N);     % 初始化系数向量
            errors = zeros(1, N);           % 初始化误差信号
            centers = zeros(1, N);          % 初始化量化中心存储
            centers_size = zeros(1, N);     % 记录中心数量变化
            centers_size(1) = 1;            % 初始中心数量
            
            % 初始化第一个中心和系数
            coefficients(1) = obj.step_size * desired_signal(1);
            centers(1) = input_signal(1);
            
            for i = 2:(N - 2*obj.order)
                % 提取当前输入向量
                current_vec = input_signal(i:i + obj.order - 1);
                distances = [];  % 存储与所有中心的距离
                
                % 计算与所有现有中心的距离
                for j = 1:(centers_size(i-1))
                    % 提取量化中心向量
                    center_vec = centers(j:j + obj.order - 1);
                    
                    % 计算核函数值并加权累加
                    k_val = obj.kernel(current_vec, center_vec, obj.sigma);
                    estimates(i) = estimates(i) + coefficients(j) * k_val;
                    
                    % 计算欧氏距离
                    dist_val = norm(current_vec - center_vec);
                    distances = [distances, dist_val];
                end
                
                % 计算误差信号
                errors(i) = desired_signal(i) - estimates(i);
                
                if ~isempty(distances)
                    % 找到最小距离及其位置
                    [min_dist, min_idx] = min(distances);
                    
                    if min_dist <= obj.epsilon
                        % 更新最近中心的系数 (量化操作)
                        coefficients(min_idx) = coefficients(min_idx) + obj.step_size * errors(i);
                        centers_size(i) = centers_size(i-1);  % 中心数量不变
                    else
                        % 添加新量化中心
                        new_center_pos = centers_size(i-1) + 1;
                        centers(new_center_pos:new_center_pos + obj.order - 1) = current_vec;
                        coefficients(new_center_pos) = obj.step_size * errors(i);
                        centers_size(i) = centers_size(i-1) + 1;  % 中心数量增加
                    end
                end
            end
        end
    end
end
