function[y,e,time] = KRLS(nsig,d,tap_num,delay,sps)
if nargin < 5
    sps = 2;
end
a = 0.1; %tap_num越多，a越小
beta = 0.99;
lambda = 0.01;
u_all = buffer(nsig,tap_num,tap_num-sps,"nodelay");
%%迭代滤波%%
e=zeros(size(u_all,2),1);%误差矩阵
y=zeros(size(u_all,2)-1,1);%滤波输出
% z(1) = step*d(1);
c(:,1) = u_all(:,1);
q(1) = 1/(lambda*beta+1);
alpha(1) = q(1)*d(1);
tic;m=1;
for i=2:size(u_all,2)%开始迭代更新参数
    k_gauss = exp(-a*sum((c-repmat(u_all(:,i),1,i-1)).^2));
    y(m) = k_gauss*alpha;%输出
    e(m) = d(i+delay)-y(i-1);%误差
    z = q*k_gauss.';
    r = lambda*beta^i+1-z'*k_gauss.';
    q = r^-1*[q*r+z*z.' -z; -z.' 1];
    alpha = [alpha-z*r^-1*e(m); r^-1*e(m)];
    c = [c u_all(:,i)];
    m = m+1;
end
time = toc;
end