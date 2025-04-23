function [res] = Build_Cut2nd_VolterraInput(input_X,P)
% Note：这种情况下，权重的抽头数量需要重新定义！！！！，再采用ABS的方式继续简化？
% 构建2阶截断的Volterra 输入
% input_X: 输入信号向量
% order: Volterra 级数
% 只适用于二阶
% P 为最大乘积间隔，tapslen_1为一阶的记忆长度，input_X为二阶记忆长度
res = [];
u=1;
% L=min(tapslen_1+P,length(input_X));

for i = 1:length(input_X)
    L=min(i+P,length(input_X));
    for j = i:L
        res = [res input_X(i)*input_X(j)];
        index_mat_nonlinear1(u, :) = [i, j];
        u=u+1;
    end
end

%     for i = 1:length(input_X)
%         for j = i:length(input_X)
%             for m = j:length(input_X)
%                 res = [res abs(input_X(i)+input_X(j)+input_X(m))];
%                 index_mat_nonlinear1(u, :) = [i, j,m];
%                 u=u+1;
%             end
%         end
%     end

end