function channel_h_coefficient = Get_channel_coefficient(signal_origin,signal_current,times,coar_syn_point,pilot_length,num_of_windows)
% coar_syn_point：粗同步点的位置。
% pilot_length：导频信号的长度。
% num_of_windows：窗口的数量。互相关计算的窗口数量，窗口通常用于定义信号的一个特定部分，该部分将被用来执行某些计算或分析。
% 找出pilot信号
signal_ori = signal_origin(1:pilot_length);
% 进行上采样
signal_ori2 = upsample(signal_ori,times);
if coar_syn_point*times-(num_of_windows*times-1) > 0
    %length(signal_rec2) = num_of_windows*times + length(signal_ori2)，取决于窗口长度和导频信号长度
    signal_rec2 = signal_current(coar_syn_point*times-(num_of_windows*times-1):coar_syn_point*times+(num_of_windows*times)+length(signal_ori2)-1);
    [r2,~] = xcorr(signal_rec2,signal_ori2,times*2*num_of_windows-1);
    r3 = abs(r2(2*num_of_windows*times-1+1:end));
    [~,index2] = max(r3);
    index3 = index2-1+coar_syn_point*times-(num_of_windows*times-1);
else
    signal_rec2 = signal_current(1:coar_syn_point*times+(num_of_windows*times-1)+length(signal_ori2)-1);
    [r2,~] = xcorr(signal_rec2,signal_ori2,coar_syn_point*times+(num_of_windows*times-1)-1);
    r3 = abs(r2(coar_syn_point*times+(num_of_windows*times-1):end));
    [~,index2] = max(r3);
    index3 = index2-1+1;
end
% 确定最终同步位置
fin_syn_point = index3;
signal_received = signal_current(fin_syn_point:times:end);
r_phase = mod(index2,times);        % Find the phase of maximum correlation point
if r_phase == 0
    r_phase = times;
end
r_for_judge = r3(r_phase:times:end);
[~,maxlocation] = max(r_for_judge);
if maxlocation >= 12 && maxlocation+19 <= length(r_for_judge)
    channel_h_coefficient = r_for_judge(maxlocation-11:maxlocation+19);
else
    channel_h_coefficient = r_for_judge(5:14);                          % Find 10 points with maximum correlation as h1 ... h10
end
end