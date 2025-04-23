function W=PNLE_Adam(input_X,N,M,dn,w,step_list)
% N为抽头数，M为训练数,dn为期望向量,一阶矩和二阶矩必须输入为上一时刻的数值，初始为0.
m1=0;
v1=0;
m2=0;
v2=0;
m3=0;
v3=0;
phi=1e-8;
% 构建R矩阵
R=zeros(M-N+1,N);
% 滑动窗口函数
for i=1:M-N+1
    L=input_X(i:i+N);
    L=fliplr(L);
    R(i,:)=L;
end
% 选取合适长度
Y=dn(M-N+1);
scale=2/(M-N+1);
% Adam参数
beta1=0.9;
beta2=0.999;


% 随意设置一个误差初始值进入迭代
J=1;
while J>1e-5
    h1=w(1:N);
    h2=w(N+1:2*N);
    h3=w(2*N+1:end);
    E=R*h1+R.^2*h2+R.^3*h3-Y;
    % 误差
    J=scale*(E.')*E;
    %梯度矩阵
    G=zeros(3,1);
    for i=1:length(G)
        R_inv=(R.^i).';
        % H变化，梯度也会更新
        G(i)=scale*R_inv*E;
    end
    % 一阶矩和二阶矩估计
    m1=beta1*m1+(1-beta1)*G(1);
    v1=beta2*v1+(1-beta2)*G(1).^2;
    m_hat1=m1/(1-beta1);
    v_hat1=v1/(1-beta2);
    % 抽头向量更新
    h_new1=h1-step_list(1)*m_hat1/(sqrt(v_hat1)+phi);

    % 一阶矩和二阶矩估计
    m2=beta1*m2+(1-beta1)*G(2);
    v2=beta2*v2+(1-beta2)*G(2).^2;
    m_hat2=m2/(1-beta1);
    v_hat2=v2/(1-beta2);
    % 抽头向量更新
    h_new2=h2-step_list(2)*m_hat2/(sqrt(v_hat2)+phi);


    % 一阶矩和二阶矩估计
    m3=beta1*m3+(1-beta1)*G(3);
    v3=beta2*v3+(1-beta2)*G(3).^2;
    m_hat3=m3/(1-beta1);
    v_hat3=v3/(1-beta2);
    % 抽头向量更新
    h_new3=h3-step_list(3)*m_hat3/(sqrt(v_hat3)+phi);

    w=cat(1,h_new1,h_new2,h_new3);
end
W=w;

% 输出后，经历一个for循环，将所有的接收信号与w进行均衡