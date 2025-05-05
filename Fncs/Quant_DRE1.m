function [sigout] = Quant_DRE1(sigin,h,Qbit_num)
%% 数字分辨率增强器(Digital Resolution Enhancer)
%% 根据论文《Low-Resolution Digital Pre-Compensation Enabled by Digital Resolution Enhancer》实现
%% 作者：cc
%% 输入参数：
%%   sigin - 输入信号
%%   h - 系统冲激响应(用于误差滤波)
%%   Qbit_num - 量化比特数
%% 输出：
%%   sigout - 增强后的输出信号

% 确保输入信号为行向量
if size(sigin,1) > 1
    sigin = sigin.';
end

% 提取信号的实部(假设输入为实信号)
sigini=real(sigin);

% 计算信号动态范围
I_max = max(abs(sigini));  % 信号最大绝对值
I_min = -I_max;            % 信号最小值为最大值的负数

% ================== 量化参数初始化 ==================
step = (I_max-I_min)/Qbit_num;  % 计算量化步长

% 生成量化区间划分点(需要Qbit_num-1个阈值)
for ii = 1:Qbit_num-1
    ipartition(1,ii) = I_min+step*ii;  % 等间距划分量化区间
end

% 生成量化码本(各区间中心值)
for jj = 1:Qbit_num
    icodebook(1,jj)  = I_min + step/2 + step*(jj-1);  % 每个区间的中心值
end

% ================== 初始量化过程 ==================
% 执行量化操作，得到：
% Iorder - 量化级别(1到Qbit_num)
% iroundsig - 量化后的信号
[Iorder,iroundsig]=quantiz(sigini,ipartition,icodebook);  
Iorder = Iorder+1;  % 调整量化级别索引从1开始

% 计算初始量化误差(量化信号 - 原始信号)
I_error = iroundsig-sigini;

% ================== 误差动态调整 ==================
n=length(icodebook);    % 码本长度(等于Qbit_num)
m=length(h);            % 系统响应长度
nvtb=length(sigin);     % 输入信号长度

% 扩展误差数组用于滑动窗口处理(首尾相接避免越界)
I_error = [I_error,I_error(1:m)];  % 在末尾添加前m个误差
Iorder = [Iorder,Iorder(1:m)];     % 对应的量化级别扩展

% 初始化候选误差存储矩阵
Qro1 = zeros(3,m);  % 非边界情况的3种候选误差调整
Qsq1 = zeros(3,m);  % 对应的滤波后误差
Qro2 = zeros(2,m);  % 边界情况的2种候选误差调整
Qsq2 = zeros(2,m);  % 对应的滤波后误差

% 主处理循环(从m+1到nvtb+m)
for ii = m+1:nvtb+m
    % 判断当前量化级别是否在中间区域(非边界)
    if Iorder(ii) ~= 1 && Iorder(ii) ~= Qbit_num  
        %% 情况1：非边界量化级(有3种调整可能)
        
        % 生成3种候选误差调整方案：
        % 1. 最后一位减step
        % 2. 保持原误差
        % 3. 最后一位加step
        Qro1(1,:) = I_error(ii-m+1:ii) - [zeros(1,m-1),step];
        Qro1(2,:) = I_error(ii-m+1:ii);
        Qro1(3,:) = I_error(ii-m+1:ii) + [zeros(1,m-1),step];
        
        % 对每种候选误差进行系统响应滤波
        Qsq1(1,:) = conv(Qro1(1,:),h,'same');  % 使用'same'保持长度一致
        Qsq1(2,:) = conv(Qro1(2,:),h,'same');
        Qsq1(3,:) = conv(Qro1(3,:),h,'same');
        
        % 选择总误差最小的调整方案
        [~,num] = min(sum(abs(Qsq1),2));  % 计算各方案绝对误差和
        I_error(ii) = Qro1(num,m);       % 更新当前误差
        
    else
        %% 情况2：边界量化级(最高/最低级，只有2种调整可能)
        
        % 确定调整方向：最高级减step，最低级加step
        adjust_sign = sign(Iorder(ii)-n/2);  % 计算符号因子
        
        % 生成2种候选误差调整：
        % 1. 调整step
        % 2. 保持原误差
        Qro2(1,:) = I_error(ii-m+1:ii) - [zeros(1,m-1),step].*adjust_sign;
        Qro2(2,:) = I_error(ii-m+1:ii);
        
        % 滤波处理
        Qsq2(1,:) = conv(Qro2(1,:),h,'same');
        Qsq2(2,:) = conv(Qro2(2,:),h,'same');
        
        % 选择最优方案
        [~,num] = min(sum(abs(Qsq2),2));
        I_error(ii) = Qro2(num,m);  % 更新当前误差
    end
end

% ================== 合成输出信号 ==================
% 重组误差序列(去除扩展部分)
I_error1 = [I_error(end-m+1:end),I_error(m+1:nvtb)];

% 将优化后的误差叠加到原始信号
sigini1 = I_error1 + sigini;  % 增强后的信号 = 原始信号 + 优化误差

% 确保输出信号维度与输入一致
if size(sigout,2) > 1
    sigout = sigout.';
end

end
