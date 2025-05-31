clc;clear;close all;
% 基于增益分配的自适应滤波算法（DMSFLNN-GR）
% ============== 参数初始化 ==============
run_number = 3;   % 独立运行次数
N = 30000;        % 总样本数
SA_N = 8;         % 子带数量
M = 65;           % 滤波器长度
LLL = 40;         % 输入向量长度
P = 2;            % 三角函数扩展阶数
LL = 2*LLL*P + LLL + 1;  % 扩展特征维度

% 预分配内存
y = zeros(1,N); y1 = zeros(1,N); E1 = zeros(1,N);
Input = zeros(1,LL); HK = zeros(SA_N,M);
E2 = zeros(1,N); d_1 = zeros(1,N);

% ============== 加载滤波器组 ==============
load filter_bank_8.mat  % 加载8通道滤波器组
for m = 1:M
    HK(:,m) = hk(:,m);  % 存储分析滤波器
end

% ============== 主循环 ==============
for i = 1:run_number
    disp(['运行: ' num2str(i)]);
    TFLNN_w = zeros(LL,1);  % 自适应滤波器权重
    
    % ----- 信号生成 -----
    input1 = unifrnd(-1,1,1,N);  % 均匀白噪声
    input = filter(1,[1,-0.9],input1);  % 通过IIR滤波器
    
    % ----- 非线性系统（带记忆）-----
    for ii = 1:N
        % 分段非线性系统（随ii增大记忆增强）
        if ii==1
            y1(ii)=0.9*input(ii)+0.8*sin(pi*input(ii));
        elseif ii>15  % 完整记忆模式
            y1(ii)=0.9*input(ii)+0.8*sin(pi*input(ii))+...
                   0.5*cos(pi*input(ii-1))+0.1*sin(pi*input(ii-2))+...
                   0.08*cos(pi*input(ii-4))-0.3*cos(pi*input(ii-15));
        ... % 其他阶段省略
        end
    end
    
    d_n = awgn(y1,30);  % 添加30dB噪声
    
    % ----- 噪声功率估计 -----
    noise = var(d_n-y1);  % 主带噪声方差
    for ii = 1:SA_N
        noise_s(ii) = noise * norm(HK(ii,:),2)^2;  % 子带噪声功率
    end
    
    % ----- 自适应滤波变量初始化 -----
    flag = 0;  % 更新控制
    Cw = 100*ones(LL,1)/LL;  % 权重协方差矩阵（初始均匀分布）
    N_e = zeros(1,SA_N);     % 子带误差功率
    V_d = zeros(LL,1);       % 子带期望功率
    V_y = zeros(LL,1);       % 子带输出功率
    N_nois = zeros(SA_N,N);  % 时变噪声功率
    v2 = zeros(LL,1);
    % ===== 核心处理循环 =====
    for p = 3*LLL-1:N  % 滑动窗口处理
        
        % 获取期望片段
        D_n1 = y1(p:-1:p-M+1);  % 无噪声期望
        D_n = d_n(p:-1:p-M+1);  % 含噪声期望
        
        % ----- 三角函数扩展 -----
        S = zeros(LL,M);  % 特征矩阵初始化
        for pp = 0:(M-1)
            Input = input(p-pp:-1:p-pp-LLL+1);  % 当前输入向量
            
            % 计算P阶三角函数扩展
            for j=1:P
                for jj=1:LLL
                    B(j,jj) = sin(j*pi*Input(jj));  % j阶正弦
                    C(j,jj) = cos(j*pi*Input(jj));  % j阶余弦
                end
                D(j,:) = [B(j,:) C(j,:)];  % 组合扩展项
            end
            
            % 重组特征（P=2时4组）
            F = zeros(4,LLL);
            for q = 1:LLL
                F(1,q) = D(1,q);       % 1阶正弦
                F(2,q) = D(1,LLL+q);   % 1阶余弦
                F(3,q) = D(2,q);       % 2阶正弦
                F(4,q) = D(2,LLL+q);   % 2阶余弦
            end
            F1 = reshape(F,1,[]);      % 特征展开
            F2 = [1 Input F1]';        % 添加常数项
            S(:,pp+1) = F2;            % 存入特征矩阵
        end
        
        % ----- 子带分解 -----
        s_h = S * HK';  % 输入信号子带分解
        d_h = HK * D_n'; % 期望信号子带分解
        
        % ----- 增益分配自适应滤波 -----
        if flag == 0 || flag == SA_N  % 块更新触发
            flag = 1;  % 重置标志
            TFLNN_del_w = zeros(LL,1);  % 权重增量
            Cd = 0;  % 协方差修正项
            
            for kk = 1:SA_N  % 遍历子带
                % 子带输出与误差
                Y(kk) = TFLNN_w' * s_h(:,kk);
                GR_e(kk,p) = d_h(kk) - Y(kk);
                
                % 递归功率估计（指数平滑）
                ff = 1-1/(1*LL);  % 平滑因子
                N_e(kk) = ff*N_e(kk) + (1-ff)*GR_e(kk,p)^2;  % 误差功率
                V_d(kk) = ff*V_d(kk) + (1-ff)*d_h(kk)^2;      % 期望功率
                V_y(kk) = ff*V_y(kk) + (1-ff)*Y(kk)^2;       % 输出功率
                
                % 噪声功率跟踪（基于功率关系）
                N_nois(kk,p) = V_d(kk)*N_e(kk)/(N_e(kk)+V_y(kk));
                if N_nois(kk,p) <= 0  % 防止负功率
                    N_nois(kk,p) = -N_nois(kk,p);
                end
                
                % 卡尔曼增益计算
                g(:,kk) = Cw.*s_h(:,kk) / (sum(s_h(:,kk).*Cw.*s_h(:,kk)) + N_nois(kk,p));
                
                % 累积权重更新
                TFLNN_del_w = TFLNN_del_w + g(:,kk)*GR_e(kk,p);
                Cd = Cd + s_h(:,kk).*g(:,kk);  % 协方差修正
            end
            
            % ----- 权重协方差更新（三选一）-----
            % 方式1：混合更新（推荐）
            v2 = 0.985*v2 + 0.015*TFLNN_del_w.^2;  % 方差平滑
            Cw = Cw - Cd .* Cw + min(norm(TFLNN_del_w)^2/LL, v2);
            
            % 方式2：标量更新（简化）
            % Cw = Cw - Cd .* Cw + norm(TFLNN_del_w)^2/LL;
            
            % 方式3：指数平滑更新
            % v2 = 0.9*v2 + 0.1*TFLNN_del_w.^2;
            % Cw = Cw - Cd .* Cw + v2;
            
            % 更新权重
            TFLNN_w = TFLNN_w + TFLNN_del_w;
            
            % 主带误差计算
            d_1(p) = S(:,1)' * TFLNN_w;
            E1(p) = D_n1(1) - d_1(p);
        else
            % 非更新时刻保持权重
            flag = flag + 1;
            d_1(p) = S(:,1)' * TFLNN_w;
            E1(p) = D_n1(1) - d_1(p);
        end
    end
    E2 = E2 + E1.^2;  % 累积平方误差
end

% ============== 结果可视化 ==============
E1_vector1_new = E2/run_number;  % 平均EMSE
figure(3)
plot(10*log10(E1_vector1_new));  % dB单位
xlabel('迭代次数','fontsize',16);
ylabel('EMSE(dB)','fontsize',16);
legend('DMSFLNN-GR'); grid on;
