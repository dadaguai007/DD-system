function [res] = Build_ABF_Input(input_X,dn)
% 构建 Volterra 输入
% input_X: 输入信号向量
% order: Volterra 级数
res = [];
u=1;
% dn 值取不到当前的本身值
for i = 1:length(input_X)
    for j = 1:length(dn)-1
        res = [res input_X(i) * dn(j)];
        index_mat_nonlinear1(u, :) = [i, j];
        u=u+1;
    end
end


end








