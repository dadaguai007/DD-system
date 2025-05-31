% Matlab function to generate a pseudo-QMF cosine modulated filterbank
% Parameters
%   nbands = n of subbands
%   f_len = filter length (must be a power of 2)
%   
% Aironi Carlo 2019
% --------------------------------------------------------------------

function [a_bank,s_bank] = build_pQMF(nbands,f_len)

% nbands=8; % 参数
% f_len=128;
% Lowpass prototype definition

stopedge = 1/nbands;
passedge = 1/(4*nbands);        % omegapass start value

tol = 0.00001;
step = 0.1 * passedge;
tcost = 0;
way = -1;
pcost = 10;
exit_flag = 0;
res_all = 4096;

while exit_flag == 0
    hopt = firpm(f_len-1,[0,passedge,stopedge,1],[1,1,0,0],[3,1]);
    
    H_all = fft(hopt,res_all);
    res_1b = floor(res_all/(2*nbands));         % n. of samples each band
    H_1b = zeros(res_1b,1);
    for k = 1:res_1b
        H_1b(k) = abs(H_all(res_1b-k+2))^2 + abs(H_all(k))^2; 
    end
    tcost = max(abs(H_1b - ones(max(size(H_1b)),1)));
    if tcost > pcost
        step = step/2;
        way = -way;
    end
    if abs(pcost - tcost) < tol
      exit_flag = 1;
    end
    pcost = tcost;
    passedge = passedge + way*step;
end

% cosine modulation

a_bank = zeros(nbands,f_len);       % analysis bank
s_bank = zeros(nbands,f_len);       % synthesis bank

for m = 1:nbands
    for n = 1:f_len
        arg_i = (pi/(nbands))*((m-1)+0.5)*((n-1)-((f_len-1)/2));
        theta_i = ((-1)^(m-1))*(pi/4);
        a_bank(m,n) = 2*hopt(n)*cos(arg_i + theta_i);
        s_bank(m,n) = 2*hopt(n)*cos(arg_i - theta_i);
    end
end


end

%% 中文注释

% 生成余弦调制伪QMF滤波器组
% 参数：
%   nbands = 子带数量
%   f_len = 滤波器长度（必须为2的幂）
% 输出：
%   a_bank = 分析滤波器组 (nbands × f_len)
%   s_bank = 合成滤波器组 (nbands × f_len)

% function [a_bank, s_bank] = build_pQMF(nbands, f_len)
% 
% % ============== 低通原型滤波器设计 ==============
% stopedge = 1/nbands;          % 阻带起始频率 (归一化)
% passedge = 1/(4*nbands);      % 初始通带截止频率 (归一化)
% 
% % 迭代优化参数设置
% tol = 0.00001;                % 收敛容差
% step = 0.1 * passedge;        % 频率调整步长
% tcost = 0;                    % 当前迭代代价
% way = -1;                     % 调整方向（初始负方向）
% pcost = 10;                   % 前次迭代代价（初始设为大值）
% exit_flag = 0;                % 循环退出标志
% res_all = 4096;               % FFT分辨率（用于频响计算）
% 
% % 原型滤波器迭代优化
% while exit_flag == 0
%     % 使用Parks-McClellan算法设计等波纹滤波器
%     hopt = firpm(f_len-1, [0, passedge, stopedge, 1], [1, 1, 0, 0], [3, 1]);
%     
%     % 计算全频带频响
%     H_all = fft(hopt, res_all);
%     
%     % 计算第一子带混叠误差
%     res_1b = floor(res_all/(2*nbands));  % 每个子带的频率采样点数
%     H_1b = zeros(res_1b, 1);             % 存储混叠误差
%     for k = 1:res_1b
%         % 计算镜像频率能量和（关键优化目标）
%         H_1b(k) = abs(H_all(res_1b - k + 2))^2 + abs(H_all(k))^2; 
%     end
%     
%     % 评估当前设计（最大偏离理想值程度）
%     tcost = max(abs(H_1b - ones(size(H_1b))));
%     
%     % 方向控制：当误差增大时反转调整方向
%     if tcost > pcost
%         step = step / 2;     % 步长减半
%         way = -way;          % 方向反转
%     end
%     
%     % 收敛检查：误差变化小于容差则终止
%     if abs(pcost - tcost) < tol
%         exit_flag = 1;      % 满足收敛条件
%     end
%     
%     pcost = tcost;                 % 更新前次代价
%     passedge = passedge + way*step;% 调整通带边缘频率
% end
% 
% % ============== 余弦调制生成滤波器组 ==============
% a_bank = zeros(nbands, f_len);  % 初始化分析滤波器组
% s_bank = zeros(nbands, f_len);  % 初始化合成滤波器组
% 
% % 对每个子带进行余弦调制
% for m = 1:nbands
%     for n = 1:f_len
%         % 计算调制相位参数
%         arg_i = (pi / nbands) * ((m-1) + 0.5) * (n - 1 - (f_len-1)/2);
%         
%         % 分析滤波器相位偏移 (-1)^(m-1)*π/4
%         theta_i = ((-1)^(m-1)) * (pi/4);
%         
%         % 生成分析滤波器系数
%         a_bank(m, n) = 2 * hopt(n) * cos(arg_i + theta_i);
%         
%         % 生成合成滤波器系数（相位符号相反）
%         s_bank(m, n) = 2 * hopt(n) * cos(arg_i - theta_i);
%     end
% end
% 
% end
