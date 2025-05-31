function x_demod=Equalizer_channel(signal_origin,signal_current,channel_h_coefficient,r_phase,coar_syn_point,num_of_windows,maxlocation,times,SNR,type)
x_demod = 0;
% signal_current为接收信号，signal_ori为发送信号，signal_received是同步之后信号
if strcmpi(type,'MRC')
    % 如果最大相关位置 maxlocation 大于或等于 12，则设置位置窗口
    if maxlocation >= 12
        location_in_window = (maxlocation-10:maxlocation+19);
    else
        %使用默认窗口
        location_in_window = (5:14);
    end
    % 根据粗同步点和窗口数量计算信号中的位置
    if coar_syn_point*times-(num_of_windows*times-1) > 0
        location_in_signal = r_phase+(location_in_window-1)*times+coar_syn_point*times-(num_of_windows*times-1)-1;
    else
        location_in_signal = r_phase+(location_in_window-1)*times;      % Find the sampling point of h1 ... h30 in received signal
    end
    % 计算相位 y_phase，它可能是信号采样点的偏移量
    y_phase = mod(location_in_signal,times);  % y_phase = r_phase
    if y_phase == 0
        y_phase = times;
    end
    y = signal_current(y_phase(1):times:end); % Received signal represented by symbols
    % 计算 h1 到 h30 在接收信号 y 中的符号点位置
    location_in_y = fix(location_in_signal/times)+1; % Find the symbol point of h1 ... h30 in received signal

    x_hat = 0;
    %每个信道冲激响应系数进行迭代，执行 MRC 等化
    for mrc_i = 1:length(channel_h_coefficient)
        y_delay = y(location_in_y(mrc_i):location_in_y(mrc_i)+length(y)-max(location_in_y));
        x_hat = x_hat+channel_h_coefficient(mrc_i)*[y_delay,zeros(1,length(signal_origin)-length(y_delay))]; % MRC
    end



    % 调整信号能量
    x_hat = x_hat./norm(x_hat,2)*sqrt(length(x_hat))*sqrt(bandpower(signal_origin)); % Make the symbolic energy of x_hat = the symbolic energy of transmitted signal
    % 提取等化信号中的导频部分，并进行归一化
    pilot = x_hat(1:pilot_length);
    pilot = pilot./norm(pilot,2)*sqrt(length(pilot));
    % 判断导频信号
    pilot(pilot <= 0) = 0;
    pilot(pilot > 0) = 1;
    zero = zeros(1,zero_length);
    % 选取信号区间
    signal_data = x_hat(pilot_length+zero_length+1:pilot_length+zero_length+1+9999);
    signal_data = signal_data./norm(signal_data,2)*sqrt(length(signal_data))*sqrt(5); % Make the symbolic energy of x_hat = the symbolic energy of transmitted signal
    % 判断信号
    signal_data = (signal_data+3)/2;
    signal_data(signal_data <= 0.5) = 0;
    signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
    signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
    signal_data(2.5 < signal_data) = 3;
    x_demod = [ pilot, zero, signal_data ];

elseif strcmpi(type,'ZF')
    % % % % % % % % % % % ZF begin % % % % % % % % % % %
     %     一个以向量 cloumn 为第一列，向量 channel_h 为第一行的 Toeplitz 矩阵
    % 信道参数进行翻转
    channel_h = [fliplr(channel_h_coefficient) zeros(1,length(signal_ori)-length(channel_h_coefficient))];
    cloumn = [channel_h(1) zeros(1,length(signal_received)-1)];
    % 构建toeplitz矩阵
    toep_h = toeplitz(cloumn,channel_h);

    % 伪逆 ，ZF原则下的均衡器
    toep_h_inv = pinv(toep_h);                                           %  H is m*n , m = n
    zf_filter = fliplr(toep_h_inv(1,:));

    x_hat = conv(signal_received,zf_filter);
    x_hat = x_hat(min(length(signal_received),length(zf_filter)):end);

    % % % % % % % % % % % ZF end % % % % % % % % % % %

    % 导频判断
    pilot = x_hat(1:pilot_length);
    pilot = pilot./norm(pilot,2)*sqrt(length(pilot));
    pilot(pilot <= 0) = 0;
    pilot(pilot > 0) = 1;
    zero = zeros(1,zero_length);
    % 信号判断
    signal_data = x_hat(pilot_length+zero_length+1:pilot_length+zero_length+1+9999);
    signal_data = signal_data./norm(signal_data,2)*sqrt(length(signal_data))*sqrt(5); % Make the symbolic energy of x_hat = the symbolic energy of transmitted signal
    signal_data = (signal_data+3)/2;
    signal_data(signal_data <= 0.5) = 0;
    signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
    signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
    signal_data(2.5 < signal_data) = 3;
    x_demod = [ pilot, zero, signal_data ];

elseif strcmpi(type,'MMSE')
    % % % % % % % % % % % MMSE begin % % % % % % % % % % %
    channel_h = [fliplr(channel_h_coefficient) zeros(1,length(signal_received)-length(channel_h_coefficient))];
    cloumn = [channel_h(1) zeros(1,length(signal_received)-1)];
    toep_h = toeplitz(cloumn,channel_h);
    % 信道的自相关矩阵，包括噪声的影响，（如1/SNR决定）
    half_mmse = (toep_h*toep_h'+1/SNR*eye(size(toep_h*toep_h',1)));
    %         half_mmse = gather(half_mmse);
    half_mmse_inv = inv(half_mmse);
    %         inv_half_mmse = gpuArray(inv_half_mmse);
    %         w = toep_h'/(toep_h*toep_h'+type{3}/type{2}*eye(size(toep_h*toep_h',1)));
    % MMSE准则下的权重矩阵
    w = toep_h'/half_mmse_inv;
    %         toc(tstart)
    x_hat = w*signal_received.';
    x_hat = x_hat.';
    % % % % % % % % % % % MMSE end % % % % % % % % % % %
    % 导频判断
    pilot = x_hat(4:pilot_length+3);
    pilot = pilot./norm(pilot,2)*sqrt(length(pilot));
    pilot(pilot <= 0) = 0;
    pilot(pilot > 0) = 1;
    zero = zeros(1,zero_length);
    % 信号判断
    signal_data = x_hat(pilot_length+14:end);
    signal_data = signal_data./norm(signal_data,2)*sqrt(length(signal_data))*sqrt(5); % Make the symbolic energy of x_hat = the symbolic energy of transmitted signal
    signal_data = (signal_data+3)/2;
    signal_data(signal_data <= 0.5) = 0;
    signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
    signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
    signal_data(2.5 < signal_data) = 3;
    x_demod = [ pilot, zero, signal_data ];

elseif strcmpi(type,'LS')
    % 等化器阶数
    equal_order = type{2};
    % 权重向量
    w = zeros(equal_order,times);
    % 用于最小二乘估计矩阵
    U = zeros(equal_order,pilot_length-equal_order+1);
    Err = zeros(1,times);
    % 确定窗口大小
    headwindow = equal_order-(fix(equal_order/2)+1);
    % 根据窗口大小和粗同步点调整接收信号，可能在前面添加零
    if headwindow*times > fin_syn_point-1
        signal_addzero = [ zeros(1,headwindow*times-(fin_syn_point-1)) , signal_current ];
    else
        start_point = fin_syn_point - headwindow*times;
        signal_addzero = signal_current(start_point:end);
    end
    % 对于每个符号
    for i = 0:(times-1)
        %指定偏移下采样序列的样本数
        pilot_received = downsample(signal_addzero(1:pilot_length*times),times,i);
        % 填充矩阵 U
        for j = equal_order:pilot_length
            % 对导频信号进行移位，每一次迭代都会将一个移位的导频信号作为列向量添加到 U 矩阵中。
            % 这样，U 矩阵的每一列都包含了导频信号的一个不同移位版本。
            U(:,j-equal_order+1) = pilot_received(j:-1:j-equal_order+1)';
        end
        % 提取原始导频信号
        pilot_ori = signal_origin(1:size(U,2))';
        % 计算均衡器权重
        w(:,i+1) = (U*U')\U*pilot_ori;
        % 计算并存储误差
        error = U'*w(:,i+1)-pilot_ori;
        Err(i+1) = norm(error,2);
    end
    % 选取最小化权重
    [~,phase] = min(Err);
     signal_received_ls = downsample(signal_addzero,times,phase-1);
    % 获取对应的权重
    W = w(:,phase);
    % 卷积
    x_hat = conv(signal_received_ls,W,'valid');

    % 决定导频
    pilot = x_hat(1:pilot_length);
    pilot(pilot <= 0) = 0;
    pilot(pilot > 0) = 1;
    zero = zeros(1,zero_length);
    % 决定信号
    signal_data = x_hat(pilot_length+zero_length+1:end);
    signal_data = (signal_data+3)/2;
    signal_data(signal_data <= 0.5) = 0;
    signal_data(0.5 < signal_data & signal_data <= 1.5) = 1;
    signal_data(1.5 < signal_data & signal_data <= 2.5) = 2;
    signal_data(2.5 < signal_data) = 3;
    x_demod = [ pilot, zero, signal_data ];

end

