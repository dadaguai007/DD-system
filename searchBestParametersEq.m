function [w_opt1,w_opt2]=searchBestParametersEq(x1_range,x2_range,matchOut,refEq,skip,label,M)

clockRecovery = DspSyncDecoding( ...
    [],...         % 接收信号的采样率
    [], ...        % 接收信号的波特率
    M,...         % 接收信号的格式
    [], ...     % 上采样率
    [],...         % 时钟信号的上采样率
    skip, ...      % 误码计算起始位置
    label,...      % 参考信号
    []);



% 均衡器参数
EQ=struct();
EQ.u=0.001;
EQ.sps=2;
EQ.lamda=0.9999;
EQ.delta=0.01;

% % 参数搜索范围
% x1_range = 1:2:51;
% x2_range = 1:2:49;

% 计算三维网格数据
[X1, X2] = meshgrid(x1_range, x2_range);
Z = zeros(size(X1));

% 扫描
for i = 1:size(X1, 1)
    for j = 1:size(X1, 2)
        % 参数 1
        w1=X1(i,j);
        % 参数2
        w2=X2(i,j);
        % 目标值-根据不同的均衡器进行改造；
        % 参考抽头设置为中间
        EQ.k1=w1;
        EQ.k2=w2;
        EQ.ref=round(w1/2);
        % FFE_LMS
        [yEq,~,~] = FFE_LMS(EQ, matchOut.', refEq.');
        % decode
        [~,berEq]=clockRecovery.PAM_ExecuteDecoding(yEq);
        Z(i,j)=berEq ;
    end
end

% 找到最佳参数
[Z1,x1_ind]=min(Z);
[fval,x2_ind]=min(Z1);

w_opt1=X1(1,x1_ind);
w_opt2=X2(x2_ind,1);


% 绘制三维图
figure;
surf(X1, X2, Z);
hold on;
plot3(w_opt1, w_opt2, fval, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
colorbar;
title('参数优化结果三维图');
xlabel('参数1');
ylabel('参数2');
zlabel('目标函数值');
grid on;


end