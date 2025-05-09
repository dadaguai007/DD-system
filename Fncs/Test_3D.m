clc;clear;close all;
% 定义x和y的范围
x = linspace(-2*pi, 2*pi, 50); % 从-2π到2π，共50个点
y = linspace(-2*pi, 2*pi, 50);

% 创建网格
[X, Y] = meshgrid(x, y);

% 计算z的值
Z = sin(X) .* cos(Y);

% 绘制曲面图
surf(X, Y, Z);

% 添加图形标题和轴标签
title('三维曲面图示例')
xlabel('X轴')
ylabel('Y轴')
zlabel('Z轴')

% 添加颜色条
colorbar

% 调整视角
view(45, 45)

% 保存图形
% saveas(gcf, '3DSurfacePlot.png');