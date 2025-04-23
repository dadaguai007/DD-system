% Volterra DFE LMS
function [x1_all,x2_vol_all,x3_vol_all,...
    fb1_vol_all,fb2_vol_all,fb3_vol_all]=Build_Volterra_Dictionary(xn,dn,sps,ref,taps_list)

% 输入向量一定小于等于参考向量

%三阶volterra
% 转换为列
xn = xn(:);
dn = dn(:);
% 序列长度
L1 = length(xn);
% 期望输入
L2 = length(dn);


% 前馈输入的抽头
% 一阶
tapslen_1=taps_list(1);
% 二阶
tapslen_2=taps_list(2);
% 三阶
tapslen_3=taps_list(3);
%反馈输入抽头
fblen_1=taps_list(4);
fblen_2=taps_list(5);
fblen_3=taps_list(6);

% 每个长度都是不同的：整个数组长度，从数组抽取两个元素的长度，从数组中抽取三个元素的长度

% 参考抽头 % 参考抽头应该存在两种方法进行定义
xn = cat(1,xn(ref+1:end),dn(end-ref+1:end));
% 反馈抽头
fb = zeros(fblen_1,1);
% 前馈抽头
x1=zeros(tapslen_1,1);

% 输出
x1_all=[];
x2_vol_all=[];
x3_vol_all=[];
fb2_vol_all=[];
fb2_vol_all=[];
fb3_vol_all=[];
for idx = 1:L1
    %构建volterra输入 （步长可舍去）
    %一阶前馈输入
    x1 = cat(1,x1(sps+1:end),xn(sps*idx-sps+1:1:sps*idx));
    %二阶前馈输入 % 从x1中确定好抽头的数量
    x2 = x1(round((tapslen_1-tapslen_2)/2)+1 : end - fix((tapslen_1-tapslen_2)/2));
    % 二阶核
    x2_vol =   BuildVolterraInput(x2,2);
    %三阶前馈输入
    x3 = x1(round((tapslen_1-tapslen_3)/2)+1 : end - fix((tapslen_1-tapslen_3)/2));
    % 三阶核
    x3_vol = BuildVolterraInput(x3,3);


    % 反馈输入
    %一阶反馈输入
    fb1_vol = fb(1:fblen_1);
    %二阶反馈输入 （步长可舍去）
    fb2_vol =  BuildVolterraInput(fb(1:fblen_2),2);
    %三阶反馈输入
    fb3_vol =  BuildVolterraInput(fb(1:fblen_3),3);
    %组合所有输入
%     x_all = [x1 x2_vol x3_vol fb1_vol fb2_vol fb3_vol];
    fb_new = dn(idx);
    x1_all = [x1_all;x1];
    x2_vol_all = [x2_vol_all;x2_vol];
    x3_vol_all = [x3_vol_all;x3_vol];
    fb1_vol_all = [fb1_vol_all;fb1_vol];
    fb2_vol_all = [fb2_vol_all;fb2_vol];
    fb3_vol_all = [fb3_vol_all;fb3_vol];
    %反馈更新
    fb = [fb(2:end) fb_new];

end


end