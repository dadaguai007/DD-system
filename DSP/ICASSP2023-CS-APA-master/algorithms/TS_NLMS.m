function [hk] = TS_NLMS(par)
%TS_NLMS 此处显示有关此函数的摘要
%   此处显示详细说明
% 结合了软阈值操作的非线性最小均方算法，用于更新权重矩阵hk。
% 软阈值函数soft_th用于在权重更新中引入非线性和稀疏性，同时正则化项有助于防止过拟合。
%% unpack parameters
lmd=par.lmd;
uk=par.uk;
dk=par.dk;
b2=par.b2;

%% set stepsize
alpha=1;
delta=0.01;

%% iterations start
hk=zeros(size(uk,1),size(uk,2)+1);
%hk(:,1)=randn(size(uk,1),1);
for kk=1:size(uk,2)
    eps_k=uk(:,kk).'*hk(:,kk)-dk(kk);
    grad_fk_plus_lmd_g=uk(:,kk)*eps_k/(delta+norm(uk(:,kk))^2)+lmd*b2*(soft_th(hk(:,kk),1/b2)-hk(:,kk));
    hk(:,kk+1)=soft_th(hk(:,kk)-alpha*grad_fk_plus_lmd_g,alpha*lmd);
end
end

function [output]=soft_th(input, threshold)
output=sign(input).*max(0,abs(input)-threshold);
end
