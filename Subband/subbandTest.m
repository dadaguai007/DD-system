clc;clear;close all;

% 输入:
%   N  : 每次运行的迭代次数，数据长度
%   Nw : 滤波器长度 (阶数+1)
%   dn : 正则化参数 (γ)
%   mu : 步长参数 (μ)
% 输出:
%   err: 平均MSE的对数值(dB)
M=4;
hk= firBank(M);
nc=normrnd(0,0.01,1,N*M);       % noise at system output
xin = color1(N*M,1);             % input signal
d=filter(B,A,xin(1:end));       % desired signal
d=d+nc;                         % unknown system output
L  = M;
u  = mu;

%% 归一化子带自适应滤波算法（NSAF）


% Analysis of desired and input signals
for k=1:M
    x_temp=filter(hk(k,:),1,d);
    dsb(k,:)=x_temp(find(mod((1:length(x_temp))-1,L)==0)); % 下采样
    x_temp=filter(hk(k,:),1,xin);
    xsb(k,:)=x_temp(find(mod((1:length(x_temp))-1,L)==0)); % 下采样
end

%initializing temporary values
w_ol =zeros(1,Nw+1);      % initial coefficient vector for subband m
for m=1:M
    x_ol(m,:)=[zeros(1,Nw+1)];    % initial input vector for subband m
end


temp_error = zeros(1,Nw+1);
%Adaption Process
for k=1:N

    for m=1:M

        x_ol(m,:)=[xsb(m,k) x_ol(m,1:Nw)]; % new input vector for subband m
        x_temp=shiftdim(x_ol(m,:),1);
        w_temp=shiftdim(w_ol,1);
        e_ol(m,l,k)=dsb(m,k)-w_temp'*x_temp;   % error at subband m

        denom=(gamma+x_ol(m,:)*x_temp);
        temp_error = temp_error+( (e_ol(m,l,k)*shiftdim(x_temp,1))/denom);
    end
    w_ol = w_ol + u*(temp_error);

end
% 每个滤波器的MSE
for m=1:M
    mse_ol(m,:)=mean(e_ol(m,:,:).^2,2);    % MSE at subband m
end



%% 比列更新算法
rho= 0.05;     % 比例更新保护阈值

for k = 1:N
   
    temp_error = zeros(1, Nw+1);  % 累积更新量

    for m = 1:M
        % 4a. 更新输入向量 (移位寄存器)
        x_ol(m, :) = [xsb(m, k), x_ol(m, 1:Nw)];

        % 4b. 计算子带误差
        e_ol(m, l, k) = dsb(m, k) - w_ol * x_ol(m, :)';

        % 4c. 计算比例矩阵G
        t = [];  % 存储比例因子
        for j = 1:(Nw+1)
            % 计算保护阈值: max(ρ*全局最大系数, |当前系数|)
            tem = max(rho * max([0.01, abs(w_ol)]), abs(w_ol(j)));
            t = [t, tem];
        end
        g = zeros(1, (Nw+1));
        for f = 1:(Nw+1)
            % 归一化比例因子 (总和=Nw+1)
            g(f) = ((Nw+1) * t(f) / sum(t));
        end
        G = diag(g);  % 构建对角比例矩阵

        % 4d. 计算正则化分母
        denom = (gamma/(Nw+1)) + x_ol(m, :)*G*x_ol(m, :)';

        % 4e. 累积子带更新量
        temp_error = temp_error + (e_ol(m, l, k) * (G * x_ol(m, :)') / denom);
    end

    % 4f. 更新滤波器系数
    w_ol = w_ol + u * temp_error;
end