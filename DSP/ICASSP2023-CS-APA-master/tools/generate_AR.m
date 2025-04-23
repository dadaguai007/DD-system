function [u] = generate_AR(gamma,K)
%GENERATE_AR 此处显示有关此函数的摘要
%   此处显示详细说明
% 生成一个自回归过程数列(AR)
% 其中序列的每个元素都是前一个元素的加权和加上一个新的随机扰动项。通过标准化，确保了序列具有单位方差。
% gamma参数控制了序列中的自相关程度，K参数定义了序列的长度。
u=zeros(K,1);
u(1)=randn;
for ii=2:K
    u(ii)=gamma*u(ii-1)+randn;
end
u=u/sqrt(var(u));
end

