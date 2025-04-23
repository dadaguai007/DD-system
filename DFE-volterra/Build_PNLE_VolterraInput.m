function [res] = Build_PNLE_VolterraInput(input_X,order)
% 构建 多项式Volterra 输入
% input_X: 输入信号向量
% order: Volterra 级数
res = [];
if order == 2
    for i = 1:length(input_X)
        res=[res input_X(i).^2];
    end
elseif order == 3
    for i = 1:length(input_X)
        res=[res input_X(i).^3];
    end
end
end