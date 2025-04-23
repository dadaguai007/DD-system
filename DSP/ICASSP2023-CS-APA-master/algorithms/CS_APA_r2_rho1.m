function [hk] = CS_APA_r2_rho1(par)
%CS_APA_R2_RHO0 此处显示有关此函数的摘要
%   此处显示详细说明
% the proposed CS-APA with r=2, rho=1 and variable stepsize

%% unpack parameters
lmd=par.lmd; % the weight of the DC regularizer
eta=par.eta; % the threshold of the debiasing component
% hyperparameters for variable stepsize
C=par.C;
mu_max=par.mu_max;
alpha=par.alpha;

uk=par.uk;
dk=par.dk;

%% iterations start
delta=0.01;
% the soft shrinkage operator
soft_th=@(x,th)sign(x).*max(0,abs(x)-th);
% the derivative of the debiasing function phi_(eta,0)
phi_prime=@(x,eta)sign(x).*(abs(x)>=eta);

hk=zeros(size(uk,1),size(uk,2)+1);
Uk=zeros(size(uk,1),2);
dk_vec=zeros(2,1);
p=zeros(2,1);

for kk=1:size(uk,2)
    % update the Uk and dk_vec
    Uk=[uk(:,kk),Uk(:,1)];    
    dk_vec=[dk(kk);dk_vec(1:end-1)];    
    
    % update the variable stepsize mu
    e_k=dk_vec-Uk.'*hk(:,kk);
    Uk_pinvUkUk_ek=Uk*pinv(Uk.'*Uk+delta*eye(2))*e_k;    
    if kk==1
        p=Uk_pinvUkUk_ek;
    else
        p=alpha*p+(1-alpha)*Uk_pinvUkUk_ek;
    end
    mu=mu_max*norm(p)^2/(C+norm(p)^2);
    
    % update hk
    hk(:,kk+1)=soft_th(hk(:,kk)+mu*(Uk_pinvUkUk_ek+lmd*phi_prime(hk(:,kk),eta)),mu*lmd);    
end
end
%`CS_APA_r2_rho0` 函数实现的算法属于压缩感知（Compressed Sensing, CS）的范畴。
% 压缩感知是一种信号处理理论，它表明如果一个信号是稀疏的或者可以被表示为一个已知字典的稀疏线性组合，可以通过远少于传统香农采样定理要求的样本数量来准确地重建该信号。
% 压缩感知的概念体现在以下几个方面：
% 
% 1. **稀疏性假设**：算法假设目标信号（在这里由权重矩阵 `hk` 表示）是稀疏的，即信号中只有少数几个非零元素。
% 
% 2. **正则化**：通过正则化项（如 `lmd` 乘以去偏置分量）来促进稀疏解。正则化是压缩感知中的关键技术，用于平衡数据拟合和稀疏性约束。
% 
% 3. **迭代更新**：算法采用迭代方式更新权重向量，每次迭代都考虑了输入信号和输出信号，以及正则化和去偏置处理，这有助于逐步逼近稀疏解。
% 
% 4. **软阈值操作**：使用软阈值函数（`soft_th`）作为非线性激活函数，它在更新权重时进一步促进了稀疏性。
% 
% 5. **可变步长**：通过计算基于输入信号和误差的加权平均向量 `p` 并据此调整步长 `mu`，算法能够自适应地调整更新的幅度，这对于处理压缩感知问题中的非平稳信号尤其重要。
% 
% 6. **最小二乘问题的正则化解**：通过正则化的最小二乘问题来计算更新向量，这有助于在压缩感知中处理噪声和不完全数据。
% 
% 压缩感知算法通常用于信号恢复、图像处理、机器学习等领域，特别是在数据获取成本高或数据本身具有稀疏特性的场景中。`CS_APA_r2_rho0` 函数中的算法通过结合自适应步长、正则化和稀疏性促进技术，实现了压缩感知中信号的稀疏重建。