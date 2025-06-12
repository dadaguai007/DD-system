function xout=Duob(SpS, Nsamples, reverse, alpha)
% Dubo FTN格式（双二进制编码脉冲生成器）
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

% 扩展滤波器长度（避免边界效应）
N = Nsamples + 2*SpS;  % 实际生成长度 = 输出长度 + 两端各扩展SpS个样本

% 参数说明：
%   alpha : 滚降因子（控制带宽）
%   N     : 滤波器长度（符号数）
%   SpS   : 每个符号样本数
%   'sqrt': 生成根升余弦而非普通升余弦

p = rcosdesign(alpha,N,SpS,'sqrt');



% === 核心操作：构造双二进制脉冲 ===
% 双二进制编码的数学形式：(1 + D) 滤波器
% 时域实现：原始脉冲 + 延迟SpS样本的脉冲
% 向左滚动SpS 个符号
% 双二进制编码的传递函数 1+D
x = p + circshift(p, -SpS);

%翻转
% 当 reverse='true' 时生成匹配滤波器
if strcmp(reverse,'true')
    indrev = length(p):-1:1;
    x = x(indrev);
end

% 截取有效部分（移除两端扩展）
xout = x(SpS+1 : end-SpS); 
% 保留中间Nsamples长度，去除扩展的边界样本

end