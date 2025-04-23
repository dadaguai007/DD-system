%数字时钟恢复
data = resample(data,osr*fb,fs);  %四倍上采样，为了满足后续时钟恢复频率
f_clk_rec = cr(data.',fb,osr,0); % 注意：f_clk_rec是对fb的估计，不是对fb*osr的估计！
t_cur = 1/osr/fb*(0:length(data)-1); % 注意：插值前data1的采样率是osr*fb
t_new = 0:1/f_clk_rec/osr:t_cur(end); % 注意：插值后data2的采样率是osr*f_clk_rec
data1 = interp1(t_cur,data,t_new,'spline');
data1 = resample(data1,sps,osr);
x=pnorm(data1);
% 同步
[xn,ref2,ff] = sync(x,ref2);

% xn = xn(1:length(ref2)*5);
