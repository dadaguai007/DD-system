
function[y,e,m] = QCKLMS(input,d,step,tap_num,delay,sps)
if nargin < 6
    sps = 1;
end
a1 = 0.1;
epsilon = 3;
u1 = buffer(input(:,1),tap_num,tap_num-sps,"nodelay");
e=zeros(size(u1,2),1);
y=zeros(size(u1,2)-1,1);%滤波输出
c_x(:,1) = u1(:,1);
z(1) = step*d(1);
m = 2;
for i = 2:size(u1,2)%开始迭代更新参数
    xn = u1(:,i);
    k = exp(-a1*sum(abs(c_x-repmat(xn,1,m-1)).^2));
    y(i-1) = z*k.';
    e(i-1) = d(i+delay)-y(i-1);
    [dis,idx] = min(sqrt(sum(abs(c_x-repmat(u1(:,i),1,m-1)).^2)));
    if dis <= epsilon
        z(idx) = z(idx)+step*e(i-1);
    else
        z(m) = step*e(i-1);
        c_x = [c_x xn];
        m = m+1;
    end
end
end