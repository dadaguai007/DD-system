function w=APfilter(input_X,dn,N,step)
% N为taps数
w=zeros(N,1);
% 构建R矩阵,能够满足M行
U=zeros(N,M);
% 滑动窗口函数
for i=1:M
    L=input_X(i:i+N);
    U(:,i)=L;
end
C = cov(U);
% 选取合适长度
dap=dn(1:M);
yap=U.'*w;
eap=dap-yap;
w=w+step*U*1/(U.'*U+C)*eap;

end