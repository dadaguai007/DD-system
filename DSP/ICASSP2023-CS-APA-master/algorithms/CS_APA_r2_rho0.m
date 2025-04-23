function [hk] = CS_APA_r2_rho0(par)
%CS_APA_R2_RHO0 此处显示有关此函数的摘要
%   此处显示详细说明
% the proposed CS-APA with r=2, rho=0 and variable stepsize

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
phi_prime=@(x,eta)sign(x).*min(1,abs(x)/eta);

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
    % p用来计算步长
    if kk==1
        p=Uk_pinvUkUk_ek;
        %这个结果是通过正则化的最小二乘问题计算得到的，表示在第一次迭代中，步长的调整完全基于当前的误差和输入信号。
    else
        p=alpha*p+(1-alpha)*Uk_pinvUkUk_ek;
        %alpha 的作用是给予之前迭代的 p 一定的权重，同时让新的 Uk_pinvUkUk_ek 也有贡献。
        % 这种加权平均是一种指数移动平均的形式，它可以平滑 p 的变化，使得步长的调整更加平稳。
    end
    % 可变步长mu的更新
    mu=mu_max*norm(p)^2/(C+norm(p)^2);
    
    % update hk
    hk(:,kk+1)=soft_th(hk(:,kk)+mu*(Uk_pinvUkUk_ek+lmd*phi_prime(hk(:,kk),eta)),mu*lmd);    
end
end
%

% Uk_pinvUkUk_ek=Uk*pinv(Uk.'*Uk+delta*eye(2))*e_k;
% 
% 
% 执行了以下步骤：
% 
% 1. `Uk.'*Uk`：首先计算矩阵 `Uk` 的转置与其自身的乘积。这通常被称为 Gram 矩阵，在信号处理和优化问题中，它表示输入信号的自相关性。
% 
% 2. `delta*eye(2)`：将一个很小的正数 `delta` 乘以一个 2x2 的单位矩阵。这个操作增加了一个正则化项，以确保矩阵是可逆的，即使在 `Uk` 的列不是完全线性独立的条件下。
% 
% 3. `pinv(Uk.'*Uk+delta*eye(2))`：计算上一步得到的矩阵的伪逆（Moore-Penrose 逆）。伪逆在矩阵可能不是满秩时仍然有效，并且可以用于求解线性方程组。
% 
% 4. `Uk*...*e_k`：将 `Uk` 与上一步计算得到的伪逆矩阵相乘，再与误差向量 `e_k` 相乘。这个操作本质上是在求解一个正则化的最小二乘问题，即：
% 
% 这个步骤的作用是计算一个正则化的最小二乘解，它代表了在正则化约束下，使得 `Uk` 的列向量线性组合最接近误差向量 `e_k` 的系数向量。在压缩感知和稀疏信号恢复中，这种方法可以用于估计信号的稀疏表示，同时避免过拟合。
% 
% 在 `CS_APA_r2_rho0` 函数的上下文中，`Uk_pinvUkUk_ek` 是 APA 更新向量的一部分，它将用于更新权重向量 `hk`，结合了输入信号的历史信息和当前的误差信息。
% 通过这种方式，算法在每一步迭代中都考虑到了信号的当前和历史信息，以期望找到稀疏且准确的信号表示。