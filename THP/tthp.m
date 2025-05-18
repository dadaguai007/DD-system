% 清除工作区变量和命令窗口
clc;
clear all;

% 系统参数设置
N = 1.6;                    % 模运算的模数
x = ones(8,1)*0.75;         % 生成8点原始信号序列，每个元素为0.75
syms z;                     % 声明符号变量z（用于Z域分析）
channel = 1 - z^(-1);       % 定义信道传递函数：H(z) = 1 - z⁻¹（典型ISI信道模型）


% x=N.*rand(10,1,'double')-N/2;
% syms z;
% channel=1-z^{-1}+3*z^{-2}+0.1*z^{-10};
%% ========== 无预编码的传输仿真 ==========
% 将时域信号转换为Z域表达式
z_x = 0;                    % 初始化Z变换结果
for i = 1:length(x)
    z_x = z_x + x(i)*z^(1-i); % 逐个元素累加计算Z变换：X(z) = Σx(n)z⁻ⁿ
end

% 信号通过信道传输
z_output = z_x * channel;   % 信道输出Y(z) = X(z)·H(z)
o = iztrans(z_output);      % 逆Z变换得到时域表达式
output_without_precoder = delta2sequence(o, length(x)); % 转换为离散序列

%% ========== THP预编码处理 ==========
% 信道逆预均衡
channel_inverse = 1/channel; % 计算信道逆：H⁻¹(z) = 1/(1 - z⁻¹)
z_pre_equalised = z_x * channel_inverse; % 预均衡：X'(z) = X(z)·H⁻¹(z)
o = iztrans(z_pre_equalised);          % 逆Z变换
pre_equalised = delta2sequence(o, length(x)); % 转换为离散序列

% 模运算控制信号幅度
pre_equalised_with_mod = modulo(pre_equalised, N); % 将信号限制在[-N/2, N/2)

% 将预编码信号重新转换为Z域
zz_pre_equalised = 0;
for i = 1:length(x)
    zz_pre_equalised = zz_pre_equalised + pre_equalised_with_mod(i)*z^(1-i);
end

% 预编码信号通过实际信道
z_output = zz_pre_equalised * channel; % Y(z) = X'(z)·H(z)
o = iztrans(z_output);                % 逆Z变换
output = delta2sequence(o, length(x)); % 接收端时域信号
output_with_mod = modulo(output, N);  % 接收端模运算恢复信号

%% ========== 结果可视化 ==========
n = [0:length(x)-1]; % 时间轴

% 图1：原始信号与无预编码接收信号对比
figure
stem(n, [x, output_without_precoder]);
title('发射与接收序列对比（无预编码）');
xlabel('时间索引 n');
legend({'发射序列','接收序列'}, 'Location','northwest');
% saveas(gcf, 'pic/original.png');

% 图2：预编码信号（无模运算）
figure
stem(n, [x, pre_equalised]);
title('预编码信号对比（未加模运算）');
xlabel('时间索引 n');
legend({'原始序列','预编码信号'}, 'Location','northwest');
% saveas(gcf, 'pic/pre_eq_without_mod.png');

% （其余绘图部分类似，略）

%% ========== 自定义函数 ==========
function mod = modulo(sequence, N)
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

function s = delta2sequence(expression, len)
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
