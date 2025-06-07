function[y,e,time] = QCKRLS(nsig,d,tap_num,delay,sps)
if nargin < 5
    sps = 2;
end
a = 0.1; %tap_num越多，a越小
epsilon = 3;
beta = 0.999;
lambda = 0.01;
u_all = buffer(nsig,tap_num,tap_num-sps,"nodelay");
%%迭代滤波%%
e=zeros(size(u_all,2),1);%误差矩阵
y=zeros(size(u_all,2)-1,1);%滤波输出
c(:,1) = u_all(:,1);
da = 1;
q(1) = 1/(lambda*beta+1);
alpha(1) = q(1)*d(1);
tic;m=1;
for i=2:size(u_all,2)%开始迭代更新参数
    k_gauss = exp(-a*sum(abs(c-repmat(u_all(:,i),1,m)).^2));
    y(i-1) = k_gauss*alpha;%输出
    e(i-1) = d(i+delay)-y(i-1);%误差
    [dis,idx] = min(sqrt(sum(abs(c-repmat(u_all(:,i),1,m)).^2)));
    if dis <= epsilon
        xi = zeros(m,1);
        xi(idx) = 1;
        da = da+xi*xi';
        k = exp(-a*sum(abs(c-repmat(c(:,idx),1,m)).^2));
        q_new = q-q(:,idx)*(k*q)/(1+k*q(:,idx));
        alpha = alpha+q(:,idx)*(d(i+delay)-k*alpha)/(1+k*q(:,idx));
        q = q_new;
    else
        c = [c u_all(:,i)];
        z = q'*k_gauss.';
        zda = q*da*k_gauss.';
        r = lambda*beta^i+1-k_gauss*zda;
        q = r^-1*[q*r+zda*z.' -zda; -z.' 1];
        alpha = [alpha-zda*r^-1*e(m); r^-1*e(m)];
        da = [da zeros(size(da,1),1);zeros(1,size(da,2)) 1];
        m = m+1;
    end
    
end
time = toc;
end