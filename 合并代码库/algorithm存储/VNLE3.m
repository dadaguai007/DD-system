function[y,e,w,time] = VNLE3(nsig,d,step,in_len,delay,sps,mf)
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
nsig_len = length(nsig);
d_len = length(d);
u = buffer(nsig,in_len,in_len-sps,"nodelay");

%%vnle-LMS迭代滤波%%
e=zeros(floor(nsig_len/2)-in_len,1);%误差矩阵
w=zeros((in_len^2+3*in_len)/2,1);%滤波器长度系数
index_mat=zeros((in_len^2+in_len)/2,2);
y=zeros(size(u,2),1);

m = 1;
for j = 1:in_len
    for k = j:in_len
        index_mat1(m,:) = [j,k];
        m = m+1;
    end
end
m = 1;
in_len_3 = floor(in_len/2);
for k = 1:in_len_3
    for j = k:in_len_3
        for q = j:in_len_3
            index_mat2(m,:) = [k,j,q];
            m = m+1;
        end
    end
end
w3 = zeros(m-1,1);
w = [w;w3];
tic
for i=1:size(u,2)%开始迭代更新参数
    x=u(:,i);%输入
    x1 = x(index_mat1(:,1)).*x(index_mat1(:,2));
    u1 = x(1:in_len_3);
    x2 = u1(index_mat2(:,1)).*u1(index_mat2(:,2)).*u1(index_mat2(:,3));
    xn = [x;x1;x2];
    y(i)=w'*xn;%输出
    if i > d_len-delay
        [~,d(i+delay)] = quantiz(y(i),partition,codebook);
    end
    e(i) = d(i+delay) - y(i);
    w=w+2*step*e(i)'*xn;
end
time = toc;
end