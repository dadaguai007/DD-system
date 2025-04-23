% 定义拟合函数 A*sinh(B*X)
func = @(params, X) params(1) * sinh(params(2) * X);

% 输入数据点 (X_i, Y_i)，这里是示例数据，您可以替换为您的实际数据
X_data = (-10:10)';
Y_data = X_data.^3; % 示例数据，您需要替换为实际数据

% 初始参数的猜测值，注意这里的顺序和定义的函数参数一致
initial_guess = [1, 1];

% 使用 lsqcurvefit 进行拟合，得到最佳参数 A 和 B
params_best_fit = lsqcurvefit(func, initial_guess, X_data, Y_data);

% 提取最佳参数 A 和 B
A_best_fit = params_best_fit(1);
B_best_fit = params_best_fit(2);

% 打印最佳参数
disp(['Best-fit A: ', num2str(A_best_fit)]);
disp(['Best-fit B: ', num2str(B_best_fit)]);

% 绘制拟合曲线和原始数据
X_fit = linspace(-10, 10, 1000); % 生成拟合曲线的 X 值范围
Y_fit = func(params_best_fit, X_fit);

plot(X_data, Y_data, 'ro', 'MarkerSize', 8); % 原始数据
hold on;
plot(X_fit, Y_fit, 'b-', 'LineWidth', 2); % 拟合曲线
xlabel('X');
ylabel('Y');
legend('Data', 'Best-fit curve');
grid on;
hold off;