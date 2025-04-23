function [hk] = l0_NLMS(par)
%L0_NLMS 此处显示有关此函数的摘要
%   此处显示详细说明
% 实现了一种结合了L0范数惩罚的非线性最小均方算法，用于更新权重矩阵hk，
% 其中f_beta函数用于引入非线性激活，以促进稀疏解。
% 这种算法在信号处理和机器学习中常用于自适应滤波和特征选择，特别是当需要稀疏解以提高模型的解释性时。
%% unpack parameters
uk=par.uk;
dk=par.dk;

gamma=par.gamma; % weight of l0-norm
mu=par.mu; % stepsize
beta=par.beta; % parameter of f_beta

delta=0.01;

%% iterations start
hk=zeros(size(uk,1),size(uk,2)+1);
for kk=1:size(uk,2)
    e_k=dk(kk)-uk(:,kk).'*hk(:,kk);
    hk(:,kk+1)=hk(:,kk)-mu*gamma*f_beta(hk(:,kk),beta)+mu*e_k*uk(:,kk)/(delta+norm(uk(:,kk))^2);
end
end
%一个规范化因子，用于防止除以零的情况发生，并提高数值稳定性
%标准最小均方（LMS）算法中的梯度下降部分，用于最小化误差e_k
function [y]=f_beta(x,beta)
y=(abs(x)<=1/beta).*(beta^2*x+sign(x).*beta);
end