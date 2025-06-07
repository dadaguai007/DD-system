function[y,e,w,time] = LMK(nsig,d,step,tap_num,delay,sps,mf)
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
%%迭代滤波%%
e=zeros(1,size(u,2));%误差矩阵
w=zeros(1,tap_num);%滤波器长度系数
y=zeros(1,size(u,2));
sigma = 0;
beta = 0.5;
tic
for i = 1:size(u,2)%开始迭代更新参数
    Xn = u(:,i);%输入
    y(i) = w*Xn;%输出
    if i > d_len
        [~,d(i)] = quantiz(y(i),partition,codebook);
    end
    e(i) = d(i+delay)-y(i);%误差
    sigma = beta*sigma+(1-beta)*e(i)^2;
    w = w+2*step*(3*sigma-e(i)^2)*e(i)*Xn';
end
time = toc;
end