classdef mlseEqu
    % MLSE equalizer
    % 此MLSE 只适用于PAM4 ，信道长度为3，寄存器为2.
    properties
        channelLength % 信道长度
        blockLength % 回溯长度
        M % 调制格式
        symbols % 星座图
        mu; % 回溯步长基数
        % states_viterbi
        % states_mlse
    end

    methods
        function obj = mlseEqu(channelLength, blockLength,M,symbols,mu)
            obj.channelLength = channelLength;
            obj.blockLength = blockLength;
            obj.M=M;
            obj.symbols=symbols;
            obj.mu=mu;
            %             obj.states_mlse = obj.states_calc(blockLength);
            %             obj.states_viterbi = obj.states_calc(channelLength);
        end

        function getStateTrellis(obj,transfer_function)
            %% ---------------------- 构建状态转移矩阵 ----------------------
            % 生成所有可能状态组合（二维状态：[previous_symbol, current_symbol]）
            state1 = repmat( obj.symbols, 1, obj.M);   % 重复符号集M次：[ -3 -1 1 3 -3 -1 1 3 ... ]
            state2 = repelem( obj.symbols, obj.M);     % 元素重复M次：[ -3 -3 -3 -3 -1 -1 -1 -1 ... ]
            statesdiag = [state2; state1];    % 组合成状态矩阵
            states = flip(statesdiag, 2);     % 反转列顺序以匹配时间顺序

            %% ---------------------- 绘制状态图(State Diagram) ----------------------
            numStates = obj.M^2;  % 总状态数=16
            theta = linspace(0, 2*pi, numStates+1); % 将圆周长均分17个点
            statePositions = [cos(theta(1:end-1)); sin(theta(1:end-1))]'; % 状态节点极坐标转笛卡尔坐标

            figure; hold on;
            title('4-PAM双状态网格编码状态图');
            axis equal; % 保持坐标轴比例一致

            % 绘制所有状态节点及标签
            for i = 1:numStates
                plot(statePositions(i,1), statePositions(i,2), 'o', 'MarkerSize',12,...
                    'MarkerFaceColor','cyan'); % 画圆形节点
                text(statePositions(i,1), statePositions(i,2), ...
                    num2str([statesdiag(1,i), statesdiag(2,i)]), ... % 显示状态值
                    'HorizontalAlignment','center', 'FontSize',8);
            end

            % 绘制状态转移路径
            colors = lines(obj.M); % 获取4种颜色
            for i = 1:numStates
                for j = 1:obj.M
                    currentInput =  obj.symbols(j); % 当前输入符号
                    nextState = [currentInput; statesdiag(1,i)]; % 计算下一状态
                    nextIndex = find(ismember(statesdiag', nextState', 'rows')); % 查找状态索引

                    if ~isempty(nextIndex)
                        % 计算箭头坐标差
                        dx = statePositions(nextIndex,1) - statePositions(i,1);
                        dy = statePositions(nextIndex,2) - statePositions(i,2);

                        % 绘制带箭头的转移线
                        quiver(statePositions(i,1), statePositions(i,2), dx, dy, 0,...
                            'Color', colors(j,:), 'MaxHeadSize',0.3);

                        % 计算转移输出值（编码器响应）
                        f_value = transfer_function(1)* obj.symbols(j) + ...
                            transfer_function(2)*statesdiag(1,i) + ...
                            transfer_function(3)*states(2,i);

                        % 标注输入/输出标签
                        label = sprintf('%d/%.2f', currentInput, f_value);
                        text(statePositions(i,1)+0.5*dx, statePositions(i,2)+0.5*dy,...
                            label, 'FontSize',8, 'Color',colors(j,:));
                    end
                end
            end
            hold off;
        end

        function getRectStateTrellis(obj,transfer_function)
            %% ---------------------- 构建状态转移矩阵 ----------------------
            % 生成所有可能状态组合（二维状态：[previous_symbol, current_symbol]）
            state1 = repmat( obj.symbols, 1, obj.M);   % 重复符号集M次：[ -3 -1 1 3 -3 -1 1 3 ... ]
            state2 = repelem( obj.symbols, obj.M);     % 元素重复M次：[ -3 -3 -3 -3 -1 -1 -1 -1 ... ]
            statesdiag = [state2; state1];    % 组合成状态矩阵
            states = flip(statesdiag, 2);     % 反转列顺序以匹配时间顺序
            %% ---------------------- 构建网格图(Trellis Diagram) ----------------------
            numTimeSteps = 2;    % 时间步数
            ySpacing = 1;        % 状态间垂直间距
            xSpacing = 1;        % 时间步间水平间距
            numStates = obj.M^2;  % 总状态数=16
            statePositionsTrellis = zeros(numStates, numTimeSteps, 2); % 状态坐标矩阵

            % 计算网格节点坐标
            for t = 1:numTimeSteps
                for i = 1:numStates
                    statePositionsTrellis(i,t,:) = [(t-1)*xSpacing, (i-1)*ySpacing];
                end
            end

            figure; hold on;
            title('4-PAM双状态网格编码时序图');
            for t = 1:numTimeSteps
                for i = 1:numStates
                    % 绘制状态节点
                    plot(squeeze(statePositionsTrellis(i,t,1)),...
                        squeeze(statePositionsTrellis(i,t,2)), 'o',...
                        'MarkerSize',10, 'MarkerFaceColor','cyan');
                    % 标注状态值
                    text(squeeze(statePositionsTrellis(i,t,1)),...
                        squeeze(statePositionsTrellis(i,t,2)),...
                        num2str([statesdiag(1,i), statesdiag(2,i)]),...
                        'HorizontalAlignment','center', 'FontSize',6);
                end
            end
            colors = lines(obj.M); % 获取4种颜色
            % 绘制时间步间转移路径
            for t = 1:numTimeSteps-1
                for i = 1:numStates
                    for j = 1:obj.M
                        currentState = statesdiag(:,i);  % 当前状态
                        currentInput =  obj.symbols(j);       % 当前输入
                        nextState = [currentInput; currentState(1)]; % 计算下一状态
                        nextIndex = find(ismember(statesdiag', nextState', 'rows'));

                        if ~isempty(nextIndex)
                            % 绘制转移线
                            plot([statePositionsTrellis(i,t,1),...
                                statePositionsTrellis(nextIndex,t+1,1)],...
                                [statePositionsTrellis(i,t,2),...
                                statePositionsTrellis(nextIndex,t+1,2)],...
                                'Color',colors(j,:), 'LineWidth',1.2);

                            % 标注转移信息
                            mid_x = mean([statePositionsTrellis(i,t,1),...
                                statePositionsTrellis(nextIndex,t+1,1)]);
                            mid_y = mean([statePositionsTrellis(i,t,2),...
                                statePositionsTrellis(nextIndex,t+1,2)]);
                            f_value = transfer_function(1)* obj.symbols(j) + ...
                                transfer_function(2)*statesdiag(1,i) + ...
                                transfer_function(3)*states(2,i);
                            label = sprintf('%d/%.2f', currentInput, f_value);
                            text(mid_x, mid_y, label, 'FontSize',8, 'Color',colors(j,:));
                        end
                    end
                end
            end
            hold off;
        end

        function [arr,input_from_states,outputs,states]=getViterbi(obj,transfer_function)
            %% ---------------------- 维特比译码器初始化 ----------------------
            if obj.channelLength ==3
                % states = combinations(symbols, symbols).Variables; % 生成所有状态组合
                [s1, s2] = ndgrid( obj.symbols,  obj.symbols);
                states = [s2(:), s1(:)]; % 生成所有状态组合
            end
            % % 行数=最大符号-最小符号+1，列数同理
            arr = zeros( obj.symbols(end)- obj.symbols(1)+1,  obj.symbols(end)- obj.symbols(1)+1); % 状态转移索引矩阵
            counter = 1;
            % 构建状态转移查找表
            for stage =  obj.symbols  % 前一状态
                for j =  obj.symbols  % 当前输入
                    % 建立状态转移索引映射
                    % 计算矩阵索引（符号值→矩阵坐标）
                    % 前一状态映射行号  % 当前输入映射列号
                    arr(stage+1- obj.symbols(1), j+1- obj.symbols(1)) = counter;
                    counter = counter + 1;
                end
            end

            % 构建状态转移输入符号矩阵
            input_from_states = zeros(length(states), length(states));
            for stage = 1:length(states)
                for j = 1:length(states)
                    if (mod(ceil(stage/obj.M), obj.M) == mod(j, obj.M))
                        input_from_states(stage,j) =  obj.symbols(ceil(j/obj.M));
                    end
                end
            end

            % 预计算所有状态转移的输出值
            outputs = zeros(length(states), length(states));
            for stage = 1:length(outputs)
                for j = 1:length(outputs)
                    l1 = input_from_states(stage,j);
                    if l1 ~= 0
                        % 三个符号（当前输入，之前状态量【一个状态中有两个寄存】），得到输出大小
                        outputs(stage,j) = sum([l1, states(stage,:)] .* transfer_function);
                    end
                end
            end

        end


        function decoded_symbols=runViterbi(obj,states,rx_signal,No_symbols,arr,input_from_states,outputs,initialstate)
            % No_symbols % 信号长度
            % states     % 状态向量
            % rx_signal  % 接收信号
            % arr,input_from_states,outputs 维特比算法初始量
            % initialstate 用于处理信号边界

            % 初始化译码变量
            possible_states = [];
            weight_of_survivor = zeros(1, length(states)); % 幸存路径度量
            survivor_paths = zeros(length(states), No_symbols+1); % 幸存路径存储
            matrix = inf(length(states), length(states)); % 路径度量矩阵
            decoded_symbols = zeros(1, No_symbols); %解码信号
            c = 1; % 解码符号索引

            % 维特比译码主循环
            for stage = 1:No_symbols+1
                % 确定当前阶段可能状态
                if stage == 1
                    possible_states = 1; % 初始状态
                elseif stage == 2                   
                    possible_states = zeros(size(obj.symbols));
                    for p = 1:length(obj.symbols)
                        possible_states = obj.next_state_index(1,  obj.symbols(p), states, arr);
                    end
                elseif stage == No_symbols
                    possible_states = find(input_from_states(:,1) ~= 0);
                elseif stage == No_symbols+1
                    possible_states = 1; % 终止状态
                else
                    possible_states = 1:length(states);
                end

                % 更新幸存路径度量
                weight_of_survivor = inf(1, length(states));
                for j = possible_states
                    [val, previous_index] = min(matrix(:,j));
                    weight_of_survivor(j) = val;
                    survivor_paths(j, stage) = previous_index;
                end

                % 回溯机制（延迟决策）
                [~, index] = min(weight_of_survivor);
                a1 = [];
                if stage >= obj.blockLength+obj.mu && mod(stage,  obj.mu)==0
                    % 逆向回溯N+mu步
                    for iteration = stage:-1:stage-(obj.blockLength+obj.mu)+2
                        index = survivor_paths(index, iteration);
                        a1 = [index, a1];
                    end
                    % 提取解码符号
                    for n = 2:obj.mu+1
                        decoded_symbols(c) = input_from_states(a1(n-1), a1(n));
                        c = c + 1;
                    end
                end

                % 处理终止阶段
                if stage == No_symbols+1
                    index = 1;
                    a1 = [index];
                    for iteration = stage:-1:c+1
                        index = survivor_paths(index, iteration);
                        a1 = [index, a1];
                    end
                    for n = 2:length(a1)
                        decoded_symbols(c) = input_from_states(a1(n-1), a1(n));
                        c = c + 1;
                    end
                end

                % 确定当前输入符号范围
                if stage == No_symbols-1
                    input_symbol = initialstate(1);
                elseif stage == No_symbols
                    input_symbol = initialstate(2);
                else
                    input_symbol =  obj.symbols;
                end

                % 初始化首阶段度量
                if stage == 1
                    weight_of_survivor(1) = 0;
                end

                % 更新路径度量矩阵
                if stage ~= No_symbols+1
                    matrix = inf(length(states), length(states));
                    for j = possible_states
                        for input = 1:length(input_symbol)
                            k1 = obj.next_state_index(j, input_symbol(input),...
                                states, arr);
                            k2 = outputs(j, k1);
                            matrix(j, k1) = (rx_signal(stage) - k2)^2; % 欧氏距离平方
                        end
                    end
                    matrix = matrix + weight_of_survivor'; % 累积路径度量
                end
            end

        end

        function SEP=getSEP(~,tx_symbols,decoded_symbols,No_symbols)
            % 计算符号错误率
            SEP= sum(tx_symbols(3:end-2) ~= decoded_symbols(1:end-2)) / No_symbols;

        end

        %% ---------------------- 状态转移索引查询 ----------------------
        function cal = next_state_index(obj,current_state_index, input, states, arr)
            % 输入参数：
            %   current_state_index : 当前状态索引
            %   input               : 输入符号
            %   states              : 状态矩阵
            %   arr                 : 状态转移索引表
            %   symbols             : 符号集
            % 输出：
            %   下一状态的索引值

            cal = arr(input + 1 -  obj.symbols(1),...        % 输入符号对应行
                states(current_state_index,1) + 1 - obj.symbols(1)); % 当前状态前一符号对应列
        end


    end

end
