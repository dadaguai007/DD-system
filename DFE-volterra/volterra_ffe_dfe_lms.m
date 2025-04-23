function [y,e,w] = volterra_ffe_dfe_lms(sym_pam,ref_sym_pam,train_len,test_len,taps_list,step_len,delay)
% % 三阶Volterra级数LMS算法
% % sym_pam 滤波器输入信号,行向量
% % ref_sym_pam 期望信号，行向量
% % taps_list中包含的分别是一阶、二阶、三阶的记忆长度

%输入:
%   sym_pam - 滤波器的输入信号，应为行向量
%   ref_sym_pam - 期望信号，也是行向量
%   train_len - 训练序列的长度
%   test_len - 测试序列的长度
%   taps_list - 一个包含一阶、二阶、三阶记忆长度的数组，
%               以及对应的一阶、二阶、三阶反馈记忆长度
%   step_len - LMS算法的步长
%   delay - 期望信号的延迟，用于时间对齐

% 说明:
%   本函数实现了一个三阶Volterra级数的LMS算法，包括前馈和反馈部分。
%   它首先在训练阶段使用LMS算法调整权重，然后在测试阶段使用这些权重来产生输出信号。
%   训练和测试阶段都涉及构建Volterra输入，并根据LMS规则更新权重。
%
sym_pam = sym_pam(:).';
ref_sym_pam = ref_sym_pam(:).';

tapslen_1 = taps_list(1);
tapslen_2 = taps_list(2);
tapslen_3 = taps_list(3);

% 反馈记忆长度
fblen_1 = taps_list(4);
fblen_2 = taps_list(5);
fblen_3 = taps_list(6);

%初始化
% 具体来说，权重向量w的长度计算如下：
% 一阶前馈记忆长度：tapslen_1。
% 二阶前馈记忆长度：tapslen_2 * (tapslen_2 + 1) / 2。这是因为在二阶项中，每个元素与其他元素可以组合成一对，计算组合的数量时使用的是组合公式，即从tapslen_2个元素中选择2个元素的组合数。
% 三阶前馈记忆长度：tapslen_3 * (tapslen_3 + 1) * (tapslen_3 + 2) / 6。同样地，这是计算从tapslen_3个元素中选择3个元素的组合数。
% 一阶反馈记忆长度：fblen_1。
% 二阶反馈记忆长度：fblen_2 * (fblen_2 + 1) / 2。
% 三阶反馈记忆长度：fblen_3 * (fblen_3 + 1) * (fblen_3 + 2) / 6。
% 将这些长度相加，得到权重向量w的总长度。然后，使用zeros函数创建一个全零的行向量，其长度等于计算出的总长度，这就是Volterra滤波器的初始权重向量。

w = zeros(tapslen_1+tapslen_2*(tapslen_2+1)/2 + tapslen_3*(tapslen_3+1)*(tapslen_3+2)/6 +...
           fblen_1+fblen_2*(fblen_2+1)/2 + fblen_3*(fblen_3+1)*(fblen_3+2)/6 ,1);

fb = zeros(1,fblen_1);
t=zeros(tapslen_1,1);
t1=zeros(tapslen_1,1);
sps=2;
%% train 训练
for i_train = 1:train_len
    %构建volterra输入
    t=cat(1,t(sps+1:end),sym_pam(sps*i_train-sps+1:1:sps*i_train).');
    t1=cat(1,sym_pam(sps*i_train-sps+1:1:sps*i_train).',t1(1:end-sps));
    %一阶前馈输入
    x1 = sym_pam(i_train : i_train+tapslen_1-1);
    %二阶前馈输入
    x2 = x1(round((tapslen_1-tapslen_2)/2)+1 : end - fix((tapslen_1-tapslen_2)/2));
    x2_vol = step_len .* BuildVolterraInput(x2,2);
    %三阶前馈输入
    x3 = x1(round((tapslen_1-tapslen_3)/2)+1 : end - fix((tapslen_1-tapslen_3)/2));
    x3_vol = step_len .* step_len .* BuildVolterraInput(x3,3);
    %一阶反馈输入
    fb1_vol = fb(1:fblen_1);
    %二阶反馈输入
    fb2_vol = step_len .* BuildVolterraInput(fb(1:fblen_2),2);
    %三阶反馈输入
    fb3_vol = step_len .* step_len .* BuildVolterraInput(fb(1:fblen_3),3);
    %组合所有输入
    x_all = [x1 x2_vol x3_vol fb1_vol fb2_vol fb3_vol];
    
    e(i_train) = ref_sym_pam(i_train+delay) - x_all * w;
    
    %使用lms更新抽头
    w = w +  e(i_train) * step_len * x_all.';
    
    %反馈更新
    fb = [ref_sym_pam(i_train+delay) fb(1:end-1)];
    
end
% figure;plot(abs(e)) % 看误差曲线
% figure;plot(w) % 看抽头分布

fb = zeros(1,fblen_1);
%% test测试
for i_test = train_len+1:train_len+test_len
    %构建volterra输入
    x1 = sym_pam(i_test : i_test+tapslen_1-1);
    x2 = x1(round((tapslen_1-tapslen_2)/2)+1 : end - fix((tapslen_1-tapslen_2)/2));
    x3 = x1(round((tapslen_1-tapslen_3)/2)+1 : end - fix((tapslen_1-tapslen_3)/2));
    %二阶输入
    x2_vol = BuildVolterraInput(x2,2);
    %三阶输入
    x3_vol = BuildVolterraInput(x3,3);
     %一阶反馈输入
    fb1_vol = fb(1:fblen_1);
    %二阶反馈输入
    fb2_vol = BuildVolterraInput(fb(1:fblen_2),2);
    %三阶反馈输入
    fb3_vol = BuildVolterraInput(fb(1:fblen_3),3);
    %组合所有输入
    x_all = [x1 x2_vol x3_vol fb1_vol fb2_vol fb3_vol];

    y(i_test-train_len) =  x_all * w;
    
    %反馈更新
    if  fblen_1 ~= 0
        fb = [pammod(pamdemod(y(i_test-train_len),4),4) fb(1:end-1)];
    end
end

end
