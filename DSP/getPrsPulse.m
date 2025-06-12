function [xout,coeff]=getPrsPulse(SpS, Nsamples, reverse, alpha,taps)
%prs FTN格式
% 输入参数：
%   SpS      : 每个符号的样本数（默认值=2）
%   Nsamples : 输出脉冲的符号长度（默认值=16）
%   reverse  : 是否翻转脉冲 ['true'/'false']（默认='false'）
%   alpha    : 根升余弦滤波器的滚降因子（默认=0.01）
% 输出：
%   xout     : 双二进制脉冲成形滤波器系数
if nargin<1
    SpS=2;
end
if nargin<2
    Nsamples=16;
end
if nargin<3
    reverse='false';
end
if nargin<4
    alpha=0.01;
end

syms x;  % 定义符号变量x
y = (1 + x)^taps;  % 构造二项式表达式 (1+x)^N
coeff = double(flip(coeffs(expand(y))));  % 多项式展开 -> 提取系数 -> 反转顺序 -> 转为数值

K = length(coeff) - 1;                % 多项式阶数
max_delay = K;                         % 最大延迟符号数
% 扩展长度（避免边界效应）
N_extend = 2 * max_delay*SpS;   % 扩展的符号长度

% 扩展滤波器长度（避免边界效应）
N = Nsamples + N_extend;  % 实际生成长度 = 输出长度 + 两端各扩展SpS个样本

% 参数说明：
%   alpha : 滚降因子（控制带宽）
%   N     : 滤波器长度（符号数）
%   SpS   : 每个符号样本数
%   'sqrt': 生成根升余弦而非普通升余弦
p = rcosdesign(alpha,N,SpS,'sqrt');



% === 核心操作：构造prs脉冲 ===
% 时域实现：原始脉冲 + 延迟SpS样本的脉冲
% 向左滚动SpS 个符号

% 初始化输出脉冲
x = zeros(1, length(p));
% ===== 核心：多项式加权叠加 =====
for k = 0:K
    delay_samples = k * SpS;           % 延迟样本数
    weighted_pulse = coeff(k+1) * circshift(p, -delay_samples);
    x = x + weighted_pulse;            % 累加加权延迟脉冲
end


%翻转
% 当 reverse='true' 时生成匹配滤波器
% 脉冲翻转
if strcmp(reverse, 'true')
    x = flip(x);
end


% 截取有效部分 (保留中间Nsamples长度，去除扩展的边界样本)
start_idx = max_delay * SpS + 1;
end_idx = length(x) - max_delay * SpS;
xout = x(start_idx : end_idx);


end