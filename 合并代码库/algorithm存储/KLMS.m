function[y,e,time] = KLMS(nsig,d,step,tap_num,delay,sps)
if nargin < 6
    sps = 2;
end
a = 0.3; %tap_num越多，a越小
nsig = nsig(:)';
d = d(:)';
u_all = buffer(nsig,tap_num,tap_num-sps,"nodelay");
%%迭代滤波%%
e=zeros(1,size(u_all,2));%误差矩阵
y=zeros(1,size(u_all,2));%滤波输出
z(1) = step*d(1);
c(:,1) = u_all(:,1);
tic
for i=2:size(u_all,2)%开始迭代更新参数
    k_gauss = exp(-a*sum((c-repmat(u_all(:,i),1,i-1)).^2));
    y(i-1)=z*k_gauss';%输出
    e(i-1)=d(i+delay)-y(i-1);%误差
    c = [c u_all(:,i)];
    z(i) = step*e(i-1);
end
time = toc;
end