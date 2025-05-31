clc;clear;close all;
% close all; % 关闭所有图形窗口（可选）

% ============== 参数初始化 ==============
run_number = 3;   % 独立运行次数（用于结果平均）
N = 30000;        % 总样本点数
SA_N = 8;         % 子带数量
LLL = 40;         % 输入向量长度（时延点数）
P = 2;            % 三角函数扩展阶数
LL = 2*LLL*P + LLL + 1;  % 扩展后特征维度（1常量 + LLL原始输入 + 2P*LLL扩展项）
u = 1;            % 自适应步长

% 预分配内存
y = zeros(1,N); y1 = zeros(1,N); E3 = zeros(1,N);
Input = zeros(1,LL); HK = zeros(SA_N,LLL);
E2 = zeros(1,N); d_1 = zeros(1,N);

% ============== 加载滤波器组 ==============
load filter_bank_8.mat  % 载入8通道滤波器组（分析滤波器hk）
% 截取前LLL个滤波器系数
for m = 1:LLL
    HK(:,m) = hk(:,m);
end

% ============== 主循环（多次独立运行） ==============
for i = 1:run_number
    disp(['运行次数: ' num2str(i)]);
    TFLNN_w = zeros(LL,1);  % 自适应滤波器权重初始化
    
    % ----- 信号生成 -----
    input1 = unifrnd(-1,1,1,N);  % 均匀分布白噪声
    input = filter(1,[1,-0.1],input1);  % 通过IIR滤波器着色
    
    % 非线性系统（分段Sigmoid函数）
    for ii=1:N
        q(ii) = (2/3)*input(ii) - (3/10)*input(ii)^2;
        if q(ii)>0
            rho=4;  % 正输入区使用陡峭Sigmoid
        else
            rho=1/2;% 负输入区使用平缓Sigmoid
        end
        y1(ii) = 2*(1/(1+exp(-rho*q(ii))))-1;
    end
    d_n = awgn(y1,30);  % 添加30dB高斯白噪声
    
    % ----- 自适应处理 -----
    flag = 0;  % 块更新控制标志
    for p = LLL:N-(LLL-1)  % 滑动窗口处理
        
        % 获取当前期望片段（含噪声/无噪声）
        D_n1 = y1(p+(LLL-1):-1:p);  % 无噪声期望（逆序）
        D_n = d_n(p+(LLL-1):-1:p);  % 含噪声期望（逆序）
        
        % 预分配特征矩阵
        D = zeros(P,2*LLL); B = zeros(P,LLL); C = zeros(P,LLL);
        S = zeros(LL,LLL);  % 扩展特征矩阵（每列对应一个延迟点）
        
        % === 三角函数特征扩展 ===
        for pp = 0:(LLL-1)
            % 获取当前输入向量（长度为LLL）
            Input = input(p+pp:-1:p+pp-LLL+1);
            
            % 计算P阶三角函数扩展
            for j=1:P
                for jj=1:LLL
                    B(j,jj) = sin(j*pi*Input(jj));  % j阶正弦
                    C(j,jj) = cos(j*pi*Input(jj));  % j阶余弦
                end
                D(j,:) = [B(j,:) C(j,:)];  % 组合扩展项
            end
            
            % 重组扩展特征（按阶数排序）
            F = zeros(4,LLL);  % P=2时4组特征
            for q = 1:LLL
                F(1,q) = D(1,q);       % 一阶正弦
                F(2,q) = D(1,LLL+q);   % 一阶余弦
                F(3,q) = D(2,q);       % 二阶正弦
                F(4,q) = D(2,LLL+q);   % 二阶余弦
            end
            F1 = reshape(F,1,[]);      % 展开为行向量
            F2 = [1 Input F1]';        % 添加常数项和原始输入
            S(:,LLL-pp) = F2;          % 存储到特征矩阵
        end
        % 非线性因子维度：每次输入一组向量，得到更多数量的展开式；至于矩阵的大小，由子滤波器的长度所决定。
        % 即示例中出现 201*40； 是因为子滤波器的长度只有40 。
        % === 子带分解 === !!! important  输入信号 与  非线性因子相乘（非线性因子维度 与 滤波器长度 大小）
        s_h = S * HK';  % 输入信号子带分解（LL×SA_N）
        d_h = HK * D_n'; % 期望信号子带分解（SA_N×1）
        
        % === 自适应滤波与权重更新 ===
        if flag == 0 || flag == 8  % 块更新策略（每8样本更新）
            flag = 1;  % 重置标志
            TFLNN_del_w = zeros(LL,1);  % 权重增量初始化
            
            % 遍历所有子带
            for kk = 1:SA_N
                Y(kk) = TFLNN_w' * s_h(:,kk);  % 子带输出
                e(kk) = d_h(kk) - Y(kk);       % 子带误差
                % 归一化LMS更新
                TFLNN_del_w = TFLNN_del_w + (e(kk)*s_h(:,kk))/(s_h(:,kk)'*s_h(:,kk));
            end
            
            % 更新滤波器权重
            TFLNN_w = TFLNN_w + u * TFLNN_del_w;
        else
            flag = flag + 1;  % 更新计数器
        end
        
        % === 主带误差计算 ===
        d_1(p) = S(:,1)' * TFLNN_w;  % 主带输出
        E3(p) = D_n1(1) - d_1(p);    % 主带误差（使用无噪声期望）
    end
    
    % 累积平方误差
    E2 = E2 + E3.^2;
end

% ============== 结果可视化 ==============
non2_large = E2/run_number;  % 计算平均EMSE
figure(1)
plot(10*log10(non2_large));  % 绘制EMSE学习曲线(dB)
xlabel('迭代次数','fontsize',16);
ylabel('EMSE(dB)','fontsize',16);
legend('SAF-TFLNN');
grid on;
