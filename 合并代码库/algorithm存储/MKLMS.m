function[y,e,time] = MKLMS(nsig,d,step,tap_num,delay,sps)
if nargin < 6
    sps = 2;
end
a = [0.01 0.05]'; %tap_num越多，a越小
epsilon = 0.3;
nsig = nsig(:)';
d = d(:)';
u_all = buffer(nsig,tap_num,tap_num-sps,"nodelay");
%%迭代滤波%%
e=zeros(size(u_all,2),1);%误差矩阵
y=zeros(size(u_all,2)-1,1);%滤波输出
z = ones(size(a));
z = z*step*d(1);
c(:,1) = u_all(:,1);
m = 2;
tic
for i=2:size(u_all,2)%开始迭代更新参数
    k_gauss = exp(-a.*sum(abs(c-repmat(u_all(:,i),1,m-1)).^2));
    y(i-1)=trace(z*k_gauss');%sum(sum(z.*k_gauss));%输出
    e(i-1)=d(i+delay)-y(i-1);%误差
     [dis,idx] = min(sqrt(sum((c-repmat(u_all(:,i),1,m-1)).^2)));
    if dis <= epsilon
        z(:,idx) = z(:,idx)+step*e(i-1);
    else
        c = [c u_all(:,i)];
%         z = [z zeros(size(a))]+step*e(i-1)*[k_gauss ones(size(a))];
z(:,m) = ones(size(a))*step*e(i-1);
        m = m+1;
    end
end
time = toc;
end