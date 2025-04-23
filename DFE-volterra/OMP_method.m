function  [V_RF,V_D] =OMP_method(N_RF,AT)
% N_RF：射频链的数量。
% AT：初始的（基向量）集合，码本。
%初始化射频链矩阵V_RF为空，获取了码本的列数（即信号空间的维度）Ns，初始化了残差矩阵VRES为单位矩阵。
V_RF =[];
Ns = size(AT,2);
VRES  = eye(Ns);
for i = 1 : N_RF
    %当前的基向量AT，计算每个基向量与当前残差VRES的内积，然后取内积的对角元素，得到每个向量与残差的相关性。
    vi = AT'*VRES;
    vi = diag(vi*vi');
    %找到最大的值及其索引，这表示与当前残差最相关的基向量。
    [~,k] = max(vi);
    % 添加基向量
    V_RF =[V_RF AT(:,k)];
    % 移除已经选择的基向量
    AT(:,k) = [];
    V_D = inv((V_RF)'*V_RF)*V_RF'*VRES;% 使用pinv也行
    % 计算残差向量
    RES = eye(Ns)-V_RF*V_D;
    VRES = RES/norm(RES,'fro');
end
 
% 最终要得到V_D
end