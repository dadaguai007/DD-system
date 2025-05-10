function [filtercoeff, wc] = prsFilter(D, span, sps)
% PRS_POLY 生成部分响应系统（Partial Response System）的脉冲成形滤波器系数
% 
% 输入参数:
%   D     - 设计参数，控制滤波器抽头的增益/延迟因子（影响频率响应）
%   span  - 滤波器的时间跨度（以符号数为单位）
%   sps   - 每个符号的采样数（Samples Per Symbol）
%
% 输出参数:
%   filtercoeff - 生成的FIR滤波器系数（行向量）
%   wc          - 滤波器的归一化截止频率

% ================== 1. 计算滤波器抽头数 ==================
taps = span * sps;  % 总抽头数 = 时间跨度 × 每符号采样数

% ========== 2. 生成二项式展开系数（组合数系数） ==========
syms x;  % 定义符号变量x
y = (1 + x)^taps;  % 构造二项式表达式 (1+x)^N
coeff = double(flip(coeffs(expand(y))));  % 多项式展开 -> 提取系数 -> 反转顺序 -> 转为数值
% 操作说明:
%   expand(y)    - 展开多项式（例如 (1+x)^3 → 1 + 3x + 3x² + x³）
%   coeffs(...)  - 提取系数（默认按降幂排列，例如 [x^3, x^2, x, 1] → [1, 3, 3, 1]）
%   flip(...)    - 反转系数顺序为升幂排列
%   double(...)  - 将符号数值转换为双精度浮点数

clear x y;  % 清除符号变量释放内存

% ========== 3. 构造滤波器系数（含D的幂次项） ==========
filtercoeff = 1;  % 初始化滤波器系数，首项为D^0=1
for i = 1:taps
    filtercoeff = [filtercoeff, D^i];  % 逐步追加D的幂次项：[1, D, D^2, ..., D^taps]
end

% ========== 4. 组合系数（二项式系数 × D幂次项） ==========
filtercoeff = filtercoeff .* coeff;  % 逐元素相乘，得到最终滤波器系数
% 数学形式：H(z) = Σ [C(k) * D^k] * z^{-k}, 其中C(k)为二项式系数

% ========== 5. 可选：归一化滤波器系数 ==========
% filtercoeff = filtercoeff / sum(filtercoeff);  % 注释掉的归一化操作（按需启用）

% ========== 6. 计算截止频率 ==========
wc = getcutoff(filtercoeff, 3);  

end
