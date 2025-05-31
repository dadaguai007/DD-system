function out = subband_qr_kalman(mic, spk, frame_size)
% 基于QR-平方根卡尔曼滤波的子带自适应回声消除器
% 输入:
%   mic - 麦克风信号(含回声)
%   spk - 参考信号(扬声器播放)
%   frame_size - 处理帧长(采样点数)
% 输出:
%   out - 回声消除后的信号

    % 初始化输出
    out_len = min(length(mic), length(spk));
    out = zeros(out_len, 1);
    out_num = floor(out_len/frame_size);  % 完整帧数
    st = init(frame_size);  % 初始化参数
    
    % 分帧处理
    for i = 1 : out_num
        % 获取当前帧
        mic_frame = mic((i-1)*frame_size+1 : i*frame_size);
        spk_frame = spk((i-1)*frame_size+1 : i*frame_size);
        
        % 核心处理
        [st, out_frame] = saf_process(st, mic_frame, spk_frame);
        
        % 输出拼接
        out(1+(i-1)*frame_size:i*frame_size) = out_frame;
    end
    
    %% ======== 初始化函数 ========
    function st = init(frame_size)
        % 基础参数
        st.frame_len = frame_size;       % 帧长度
        st.K = frame_size*2;             % FFT长度
        st.half_bin = st.K/2 + 1;        % 单边频谱点数
        st.win_len = st.K*2;             % 窗长度
        st.notch_radius = .982;          % 陷波滤波器半径
        st.notch_mem = zeros(2,1);       % 陷波滤波器状态
        
        % 加载预计算的窗函数
        win_st = load('win_para_update.mat');
        st.win_global = win_st.win;      % 分析/合成窗
        
        % 初始化窗缓冲区
        st.ana_win_echo = zeros(1, st.win_len); % 麦克风分析窗
        st.ana_win_far = zeros(1, st.win_len);  % 参考信号分析窗
        st.sys_win = zeros(1, st.win_len);      % 合成窗
        
        % 子带滤波器参数
        st.tap = 15;  % 每个子带滤波器阶数
        st.subband_in = zeros(st.half_bin, st.tap);   % 子带输入历史
        st.subband_adf = zeros(st.half_bin, st.tap);  % 子带滤波器系数
        
        % QR-平方根卡尔曼参数
        st.Ryu = ones(st.half_bin,st.tap,st.tap)*5;  % 协方差平方根矩阵
        st.w_cov = ones(st.half_bin, 1)*0.1;         % 过程噪声方差
        st.v_cov = ones(st.half_bin, 1)*0.001;       % 观测噪声方差
        st.gain = zeros(st.half_bin,st.tap);         % 增益矩阵
        
        % 非线性处理(NLP)参数
        st.Eh = zeros(st.half_bin,1);     % 误差功率平滑值
        st.Yh = zeros(st.half_bin,1);     % 回声功率平滑值
        st.est_ps = zeros(st.half_bin,1); % 当前回声功率估计
        st.spec_ave = 0.01;               % 功率谱平滑因子
        st.Pey = 0;                       % 误差-回声互相关
        st.Pyy = 0;                       % 回声自相关
        st.beta0 = 0.016;                 % 自适应步长基数
        st.beta_max = st.beta0/4;         % 最大步长
        st.min_leak = 0.005;              % 最小泄漏因子
        st.echo_noise_ps = 0;             % 回声+噪声功率估计
        st.adapt_cnt = 0;                 % 自适应计数器
        st.res_old_ps = 0;                % 残留功率记忆
        st.suppress_gain = 20;            % 回声抑制增益
        st.wiener_gain = zeros(st.half_bin,1); % 维纳增益
        st.gain_floor = ones(st.half_bin,1)*0.01; % 增益下限
    end

    %% ===== 直流陷波滤波器 =====
    function [out,mem] = filter_dc_notch16(in, radius, len, mem)
        % IIR陷波滤波器，去除直流分量
        % 输入:
        %   in - 输入信号
        %   radius - 陷波半径(接近1)
        %   len - 信号长度
        %   mem - 滤波器状态
        % 输出:
        %   out - 滤波后信号
        %   mem - 更新后的状态
        
        out = zeros(size(in));
        den2 = radius*radius + .7*(1-radius)*(1-radius); % 分母系数
        for ii = 1:len
            vin = in(ii);
            vout = mem(1) + vin;
            mem(1) = mem(2) + 2*(-vin + radius*vout);
            mem(2) = vin - (den2*vout);
            out(ii) = radius*vout;  % 滤波输出
        end
    end

    %% ===== 核心处理函数 =====
    function [st, out] = saf_process(st, mic_frame, spk_frame)
        N = st.frame_len; % 帧长
        K = st.K;         % FFT长度
        M = st.tap;       % 滤波器阶数
        
        % 1. 麦克风信号: 直流陷波
        [mic_in, st.notch_mem] = filter_dc_notch16(mic_frame, st.notch_radius, N, st.notch_mem);

        % 2. 麦克风信号: 分帧加窗
        st.ana_win_echo = [st.ana_win_echo(N+1:end), mic_in']; % 更新缓冲区
        ana_win_echo_windowed = st.win_global .* st.ana_win_echo; % 加窗
        ana_wined_echo = ana_win_echo_windowed(1:K) + ana_win_echo_windowed(K+1:2*K); % 分段相加
        fft_out_echo = fft(ana_wined_echo); % FFT变换

        % 3. 参考信号: 分帧加窗
        st.ana_win_far = [st.ana_win_far(N+1:end), spk_frame'];
        ana_win_far_windowed = st.win_global .* st.ana_win_far;
        ana_wined_far = ana_win_far_windowed(1:K) + ana_win_far_windowed(K+1:2*K);
        fft_out_far = fft(ana_wined_far, K); 

        % 4. 子带自适应滤波
        st.subband_in = [fft_out_far(1:st.half_bin)', st.subband_in(:,1:M-1)]; % 更新输入
        subband_adf_out = sum(st.subband_adf .* st.subband_in, 2); % 回声估计
        subband_adf_err = fft_out_echo(1:st.half_bin)' - subband_adf_out; % 计算误差
        
        % 5. QR-平方根卡尔曼更新 (核心改进)
        for j = 1:st.half_bin
            % 更新观测噪声方差
            st.v_cov(j) = 0.9*st.v_cov(j) + 0.1*(subband_adf_err(j)*conj(subband_adf_err(j)));
            
            % ==== 预测步骤 ====
            % 构建增广矩阵 [过程噪声; 状态协方差]
            A1 = [(eye(M).*sqrt(st.w_cov(j))), squeeze(st.Ryu(j,:,:))]';
            % QR分解获取预测后协方差
            [~, B1] = qr(A1);
            Rmu = B1(1:M,1:M)'; % 提取预测协方差平方根
            
            % ==== 更新步骤 ====
            % 构建观测增广矩阵 [观测噪声; 状态观测]
            A2 = [sqrt(st.v_cov(j)), zeros(1,M); 
                  (st.subband_in(j,:) * Rmu)', Rmu'];
            % QR分解计算卡尔曼增益
            [~, B2] = qr(A2);
            r12 = B2(1,1); % 标量参数
            g_r12 = B2(1, 2:end)'; % 增益向量
            gain = g_r12 ./ (r12 + 1e-14); % 计算卡尔曼增益
            
            % 更新滤波器系数
            phi = gain .* subband_adf_err(j);
            st.subband_adf(j,:) = st.subband_adf(j,:) + phi.';
            
            % 更新协方差平方根
            st.Ryu(j,:,:) = B2(2:end,2:end)';
            
            % 更新过程噪声方差
            st.w_cov(j) = 0.9*st.w_cov(j) + 0.1*(sqrt(phi'*phi)/M);
        end
        
        % 6. 非线性处理
        [st, nlpout] = nlpProcess(st, subband_adf_err, subband_adf_out);
        
        % 7. 合成输出
        % 构建共轭对称频谱
        ifft_in = [nlpout', fliplr(conj(nlpout(2:end-1)'))]; 
        fft_out = real(ifft(ifft_in)); % IFFT转时域
        
        % 加窗处理
        win_in = [fft_out, fft_out];
        comp_out = win_in .* st.win_global;
        
        % 重叠相加合成
        st.sys_win = st.sys_win + comp_out;
        out = st.sys_win(1:N); % 输出当前帧
        
        % 更新合成窗
        st.sys_win = [st.sys_win(N+1:end), zeros(1,N)];
        st.adapt_cnt = st.adapt_cnt + 1; % 更新计数器
    end

    %% ===== 非线性处理函数 =====
    function [st, nlp_out] = nlpProcess(st, error, est_echo)
        % 计算功率谱
        st.est_ps = abs(est_echo).^2; % 回声功率
        res_ps = abs(error).^2;       % 残留功率
        
        % 计算功率变化量
        Eh_curr = res_ps - st.Eh;
        Yh_curr = st.est_ps - st.Yh;
        
        % 更新平滑功率
        Pey = sum(Eh_curr.*Yh_curr); % 误差-回声互相关
        Pyy = sum(Yh_curr.*Yh_curr); % 回声自相关
        st.Eh = (1-st.spec_ave)*st.Eh + st.spec_ave*res_ps;
        st.Yh = (1-st.spec_ave)*st.Yh + st.spec_ave*st.est_ps;
        
        % 计算总功率
        Syy = sum(st.est_ps); % 回声总功率
        See = sum(res_ps);    % 残留总功率
        Pyy = sqrt(Pyy);      % 归一化因子
        Pey = Pey/(Pyy+1e-10); % 归一化互相关
        
        % 计算自适应步长
        tmp32 = st.beta0*Syy/See;
        alpha = min(tmp32, st.beta_max);
        
        % 更新相关统计量
        st.Pyy = (1-alpha)*st.Pyy + alpha*Pyy;
        st.Pey = (1-alpha)*st.Pey + alpha*Pey;
        st.Pyy = max(st.Pyy,1); % 数值保护
        st.Pey = max(st.Pey, st.Pyy*st.min_leak); % 最小泄漏
        st.Pey = min(st.Pey, st.Pyy);             % 最大泄漏
        
        % 计算泄漏因子
        leak = st.Pey/st.Pyy;
        if leak > 0.5
            leak = 1; % 泄漏过大时禁用抑制
        end
        
        % 计算残留回声功率
        residual_ps = leak*st.est_ps*st.suppress_gain;
        
        % 更新噪声功率估计
        if st.adapt_cnt == 0
            st.echo_noise_ps = residual_ps;
        else
            st.echo_noise_ps = max(0.85*st.echo_noise_ps, residual_ps);
        end
        st.echo_noise_ps = max(st.echo_noise_ps, 1e-10); % 防除零
        
        % 计算后验信噪比
        postser = res_ps./st.echo_noise_ps - 1;
        postser = min(postser, 100); % 限幅
        
        % 更新残留功率记忆
        if st.adapt_cnt == 0
            st.res_old_ps = res_ps;
        end
        
        % 计算先验信噪比
        prioriser = 0.5*max(0, postser) + 0.5*(st.res_old_ps./st.echo_noise_ps);
        prioriser = min(prioriser, 100); % 限幅
        
        % 计算维纳增益
        st.wiener_gain = prioriser./(prioriser+1);
        st.wiener_gain = max(st.wiener_gain, st.gain_floor); % 增益下限保护
        
        % 更新残留功率估计
        st.res_old_ps = 0.8*st.res_old_ps + 0.2*st.wiener_gain.*res_ps;
        
        % 应用增益抑制残留回声
        nlp_out = st.wiener_gain.*error;
    end
end
