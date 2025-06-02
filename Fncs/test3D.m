clc;clear;close all;
% 定义目标函数（示例函数：Rastrigin函数）
objective = @(x) x(1)^2 + x(2)^2 - 10*(2 - cos(2*pi*x(1)) - cos(2*pi*x(2)));

% 初始猜测值
x0 = [0, 0];

% 参数搜索范围
x1_range = -5:0.1:5;
x2_range = -5:0.1:5;

% 计算三维网格数据
[X1, X2] = meshgrid(x1_range, x2_range);
Z = zeros(size(X1));

for i = 1:size(X1, 1)
    for j = 1:size(X1, 2)
        Z(i,j) = objective([X1(i,j), X2(i,j)]);
    end
end

% 使用fminsearch寻找最优参数
options = optimset('Display', 'iter');
[x_opt, fval] = fminsearch(objective, x0, options);

% 绘制三维图
figure;
surf(X1, X2, Z);
hold on;
plot3(x_opt(1), x_opt(2), fval, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
colorbar;
title('参数优化结果三维图');
xlabel('参数1');
ylabel('参数2');
zlabel('目标函数值');
grid on;
% 调整视角
view(45, 45)
% 输出最优结果
fprintf('最优参数值: 参数1 = %.4f, 参数2 = %.4f\n', x_opt(1), x_opt(2));
fprintf('最优目标函数值: %.4f\n', fval);