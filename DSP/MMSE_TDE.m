function [w,G]=MMSE_TDE(num_taps, channell)
% 生成G矩阵，基于信道响应和抽头数
G = getG(channell, num_taps);
E = getE(length(channell)+num_taps-1);
disp("G mat is: ")
disp(G)
disp("E mat is: ")
disp(E)
% 为了得到taps数的抽头数
% G为【num_taps+channel-1，num_taps】; 运算后为【num_taps+channel-1，num_taps+channel-1】
% 与taps 数不符
w = getMMSEW(G',E,1);
disp("w for 0 SNR is: ")
disp(w)
end

%% get e vector
function [E] = getE(Esize)
E = zeros(Esize,1);
E(1,1)=1;
end

%% (Toeplitz)矩阵
% 生成表示信道卷积操作的托普利兹(Toeplitz)矩阵
% 形成对角线排列的托普利兹结构
function [G] = getG(channel, num_tap)
% 长度为 L_h 的信道与长度为 L_w 的滤波器卷积，输出长度为 L_h + L_w - 1
% 滤波器抽头数量即为信号的长度

% 计算 G 矩阵的行数，等于信道冲激响应的长度加上均衡器的抽头数减1。
row_size = length(channel) + num_tap - 1;

col_size = num_tap;  % 矩阵列数 = 均衡器抽头数
G = zeros(row_size, col_size);  % 初始化全零矩阵
% 翻转信道响应，便于引用
% 卷积运算本质是"翻转+滑动"操作
% 矩阵构建时需要从后向前取信道系数
channel_flipped = flip(channel);
for row = 1:row_size  % 遍历每一行输出
    % 获取当前行需要填充的列范围
    [start_pt, end_pt] = get_idx(row, length(channel)-1, col_size-1);
    % 从翻转后的信道响应中提取相应样本
    channel_to_be_inserted = channel_flipped(end-(end_pt-start_pt):end);
    % 行索引过大时的边界处理
    %     边界处理原理：
    % 当行号超过均衡器长度时：
    %
    % 不再从信道尾部取系数
    % 改为从信道开头取系数（对应卷积的拖尾部分）
    % 数学意义：处理卷积的"不完全重叠"区域
    if row > col_size
        % 调整选取量：从信道开头取数据
        channel_to_be_inserted = channel_flipped(1:end_pt-start_pt+1);
    end
    % 将信道系数填充到矩阵对角线位置
    G(row, start_pt:end_pt) = channel_to_be_inserted;
end  % 结束行遍历
end
%%

% 计算G矩阵中当前行需要填充的列范围
function [start_pt, end_pt] = get_idx(row, N1plusN2, M1plusM2)
end_pt = row;  % 默认结束位置 = 当前行号
start_pt = 1;   % 默认起始位置 = 1
% 防止结束位置超过均衡器范围
% 当行号 > 均衡器长度时，填充列数不能超过均衡器维度
if row > M1plusM2 + 1
    end_pt = M1plusM2 + 1;  % 限制为最大列数
end
% 调整起始位置（处理行号超过信道长度的情况）
% 当行号 > 信道长度时，填充需要跳过前面的零值
if row > N1plusN2 + 1
    start_pt = row - N1plusN2;  % 计算新的起始位置
end

end
