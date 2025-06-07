function[y,e,w,time] = RGKLMK(nsig,d,step,tap_num,delay,sps,mf)
if nargin < 7
    mf = 'ook';
end
if strcmp(mf,'ook') == 1
    partition = 0;
    codebook = [-1,1];
else
    partition = [-2,0,2];
    codebook = [-3,-1,1,3];
end
if nargin < 6
    sps = 2;
end
d_len = length(d);
nsig = nsig(:)';
d = d(:)';
u = buffer(nsig,tap_num,tap_num-sps,"nodelay");
a = 0.1; %Gaussian kernel parameter
sigma = 0;
beta = 0.5;
p_len = 5;
p = [0 1:p_len];
c_p = sqrt((2*a).^p./factorial(p));
e=zeros(1,size(u,2));%误差矩阵
w=zeros(tap_num,p_len+1);%滤波器长度系数
y=zeros(1,size(u,2));
tic
for i = 1:size(u,2)%开始迭代更新参数
    x = u(:,i);%输入
    exp_x = exp(-a*x.^2);
    phi_x = c_p.*(x.^p.*exp_x);
    y(i) = sum(sum(w.*phi_x));%输出
    if i > d_len
        [~,d(i)] = quantiz(y(i), ...
            partition,codebook);
    end
    e(i) = d(i+delay)-y(i);%误差
    sigma = beta*sigma+(1-beta)*e(i)^2;
    w = w+2*step*(3*sigma-e(i)^2)*e(i)*phi_x;
end
time = toc;
end