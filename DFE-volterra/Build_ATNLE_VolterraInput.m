function [res] = Build_ATNLE_VolterraInput(input_X,order)
% 构建 Volterra 输入
% input_X: 输入信号向量
% order: Volterra 级数
res = [];
u=1;
if order == 2
    for i = 1:length(input_X)
        res = [res abs(input_X(i))];
        index_mat_nonlinear1(u, :) = i;
        u=u+1;
    end
elseif order == 3
    for i = 1:length(input_X)
        res = [res abs(input_X(i))*input_X(i)];
        index_mat_nonlinear1(u, :) = i;
        u=u+1;
    end
end
end