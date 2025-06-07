function[y,e,time] = DFE_KLMS(nsig,d,step,tap_num,delay,sps,mf)
if nargin < 7
    mf = 'ook';
end
if strcmp(mf,'ook') == 1
    partition = 0;
    code  book = [-1,1];
else
    partition = [-2,0,2];
    codebook = [-3,-1,1,3];
end
if nargin < 6
    sps = 2;
end
a1 = 0.1;
a2 = 5;
z1(1) = step(1)*d(1);
z2(1) = 1;
d_len = length(d);
nsig = nsig(:)';
d = d(:)';
u1 = buffer(nsig,tap_num(1),tap_num(1)-sps,"nodelay");
u2 = buffer([zeros(1,2*tap_num(2)) d(delay+2:end)],tap_num(2),tap_num(2)-1,"nodelay");
c1(:,1) = u1(:,1);
c2(:,1) = u2(:,1);
e=zeros(1,size(u1,2));%误差矩阵
y = zeros(1,size(u1,2)); %滤波输出
tic
for i = 2:size(u1,2)%d_len-delay%;开始迭代更新参数
    k_gauss1 = exp(-a1*sum((c1-repmat(u1(:,i),1,i-1)).^2));
    k_gauss2 = exp(-a2*sum((c2-repmat(u2(:,i),1,i-1)).^2));
    y(i-1)=z1*k_gauss1'+z2*k_gauss2';%输出
    e(i-1)=d(i+delay)-y(i-1);%误差
    c1 = [c1 u1(:,i)];
    c2 = [c2 u2(:,i+1)];
    z1(i) = step(1)*e(i-1);
    z2(i) = step(2)*e(i-1);
end
% x_dfe = u2(:,i);
% % y = y(1:i);
% for i = d_len-delay+1:size(u1,2)
%     k_gauss1 = exp(-a1*sum((c1-repmat(u1(:,i),1,i-1)).^2));
%     k_gauss2 = exp(-a2*sum((c2-repmat(x_dfe,1,i-1)).^2));
%     y(i-1)=z1*k_gauss1'+z2*k_gauss2';%输出
%     [~,d1] = quantiz(y(i-1),partition,codebook);
%     e(i-1)=d1-y(i-1);%误差
%     c1 = [c1 u1(:,i)];
%     c2 = [c2 x_dfe];
%     x_dfe = cat(1,x_dfe(2:end),d1);
%     z1(i) = step(1)*e(i-1);
%     z2(i) = step(2)*e(i-1);
% end
time = toc;
%         w = [w_ffe w_dfe];
end