function w = getCutoff(ch, ampdB)
% 计算滤波器的归一化截止频率
% 
% 输入参数:
%   ch     - 滤波器系数（行向量或列向量）
%   ampdB  - 截止频率处的衰减值（单位：dB，如3dB对应半功率点）
%
% 输出参数:
%   w      - 归一化截止频率（范围0~1，1对应π弧度/样本）

% ========== 1. 将分贝衰减转换为线性幅度比 ==========
amp = 10^(ampdB/10); 
% 示例：当ampdB=3时，计算10^(3/10)=2，对应幅度衰减至1/sqrt(2) ≈ 0.707
% 注意：此处使用功率比计算，20*log10(1/sqrt(amp)) = -ampdB/2 dB

% ========== 2. 计算滤波器频率响应 ==========
[H, f] = freqz(ch, 1, 2^12); 
% 参数说明:
%   ch    - 分子系数（FIR滤波器）
%   1     - 分母系数（FIR滤波器分母为1）
%   2^12  - 频率响应计算点数（4096点）
% 输出说明:
%   H - 复数频率响应（0~π范围内）
%   f - 对应的归一化角频率向量（0~π）

% ========== 3. 获取最大幅度值 ==========
amx = max(abs(H)); % 提取频率响应的最大幅度（通常为直流分量）

% ========== 4. 寻找截止频率点 ==========
[~, idx] = min(abs(abs(H)-amx/sqrt(amp)));
% 计算逻辑：
%   1. amx/sqrt(amp) 是目标幅度值（如3dB对应amx/sqrt(2)）
%   2. 寻找实际幅度响应最接近目标值的频率点
% 注意：当存在多个交点时，返回第一个满足条件的频率点

% ========== 5. 计算归一化截止频率 ==========
w = f(idx) / pi; 
% 将角频率从0~π范围归一化到0~1范围（1对应Nyquist频率）
end
