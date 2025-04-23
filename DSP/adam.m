function[y,e,w] = adam(input,dn,tap_num)

nsig_len = length(input);
input = input(:);
dn = dn(:);
dn = dn(1:nsig_len-tap_num+1);
beta1 = 0.89;%0.89;
beta2 = 0.96;%0.951;
theta = 0.08;%0.9;
epsilon = 1e-6;%0.111;
w = zeros(tap_num,1); %滤波器长度系数
m = zeros(tap_num,1);
v = zeros(tap_num,1);
for j = 1:nsig_len-tap_num+1
R(j,:) = input(j+ tap_num - 1:-1:j);
end
y = zeros(1,nsig_len); %滤波输出
y(1:tap_num) = input(1:tap_num); %假设滤波器前几个点为输入值
    for i = (tap_num+1):nsig_len-tap_num+1 %开始迭代更新参数
        g = 2/(nsig_len-tap_num+1)*R'*(R*w-dn);
        m = beta1*m+(1-beta1)*g;
        v = beta2*v+(1-beta2)*g.^2;
        m1 = m/(1-beta1^(i-tap_num));
        v1 = v/(1-beta1^(i-tap_num));
        
        w = w - theta*m1./(sqrt(v1)+epsilon);
        y(i) = R(i,:)*w;
        e(i-tap_num) = dn(i) - y(i);
    end
end
