function[y,e,w,time] = LMS(nsig,d,step,tap_num,delay,sps,mf)
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
u = buffer(nsig,tap_num,tap_num-sps,"nodelay");
%%迭代滤波%%
e=zeros(size(u,2),1);%误差矩阵
w=zeros(tap_num,1);%滤波器长度系数
y=zeros(size(u,2),1);
tic
for i=1:size(u,2)%开始迭代更新参数
    Xn=u(:,i);%输入
    y(i)=w'*Xn;%输出
    if i > d_len-delay
        [~,d(i)] = quantiz(y(i),partition,codebook);
    end
    e(i)=d(i+delay)-y(i);%误差
    w=w+2*step*e(i)'*Xn;
end
time = toc;
end