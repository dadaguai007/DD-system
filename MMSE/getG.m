%% MMSE 均衡器，得到矩阵G进行均衡
% 时域MMSE均衡器
%% get G matrix
function [G] = getG(channel, num_tap)
    % 计算 G 矩阵的行数，等于信道冲激响应的长度加上均衡器的抽头数减1。
    row_size = length(channel) + num_tap -1 ; % N2+M2+1
    col_size = num_tap;
    G = zeros(row_size, col_size);
    % 翻转信道响应，便于引用
    channel_flipped=flip(channel);
    for row=1:row_size
        [start_pt, end_pt] = get_idx(row, length(channel)-1, col_size-1);
        % 从信道响应中提取相应样本
        channel_to_be_inserted = channel_flipped(end-(end_pt-start_pt):end);
        if row>col_size 
            % 行索引过大，调整选取量
            channel_to_be_inserted = channel_flipped(1:end_pt-start_pt+1);
        end
        G(row,start_pt:end_pt)=channel_to_be_inserted;
    end
end
%% get starting and ending point of the flowing channel in G
function [start_pt, end_pt] = get_idx(row, N1plusN2, M1plusM2)
    end_pt = row;
    start_pt = 1;
    if row>M1plusM2+1
        end_pt = M1plusM2+1;
    end
    if row>N1plusN2+1
        start_pt = row-N1plusN2;
    end
end


