function [res] = Build_DP_ABS_VolterraInput(input_X,W)
% Note：这种情况下，权重的抽头数量需要重新定义！！！
% 构建2阶截断的Volterra 输入
% input_X: 输入信号向量
% 只适用于二阶
res = [];
u=1;
% L=min(tapslen_1+P,length(input_X));
% 裁剪的乘法偏移W
for i = 1:W
    for j = i:length(input_X)
        res = [res abs(input_X(i)+input_X(j))];
        index_mat_nonlinear1(u, :) = [i, j];
        u=u+1;
    end
end
end