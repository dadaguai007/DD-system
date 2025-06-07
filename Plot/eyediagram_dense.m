function [fh, scr_out] = eyediagram_dense(varargin)
% EYEDIAGRAM_DENSE 生成高密度眼图
% 输入参数：
%   signal    : 输入信号向量
%   npoints   : 每段眼图的点数（周期长度）
% 可选参数：
%   "time", time_vector  : 提供时间向量生成时间轴标签
%   "midbit", true/false : 是否在比特中间标记时间轴（默认false）
%   "histogram", true/false : 是否显示直方图（默认false）
%   "xresolution", 值     : 设置x轴分辨率（像素）
%   "yresolution", 值     : 设置y轴分辨率（像素）
% 输出：
%   fh        : 图形句柄
%   scr_out   : 结构体包含处理数据
%     .x        : x轴坐标向量
%     .y        : y轴坐标向量
%     .image    : 处理后的眼图矩阵
%     .raw_data : 原始信号分箱计数矩阵

% To Do ： x轴和y轴的添加实现；直方图位置可控

%% 输入验证
switch nargin
    case 0
        help eyediagram_dense;
        error('未提供输入参数');
    case 1
        error('未提供每段点数');
    case {3,5,7,9,11}
        error('参数数量不足');
end
if nargin>12, error('参数过多'); 
end

%% 初始化参数
signal = varargin{1};
rep    = varargin{2};       % 每段点数
time   = [];                % 时间向量初始化
midbit = false;             % 中比特标记开关
plothist = false;           % 直方图开关
screen_size = [0,0];        % 分辨率初始化

% 确保信号为行向量
if size(signal,1) > size(signal,2)
    signal = signal'; 
end

% 解析可选参数
% % 遍历从第3个开始的输入参数 (前两个固定为signal和npoints)
% for i = 3:nargin
%     % 检查是否为"time"参数
%     if strcmpi("time", varargin{i})
%         tempstore = varargin{i+1};  % 获取时间向量
%         % 确保时间向量为行向量 (与signal方向一致)
%         if size(tempstore,1) > size(tempstore,2)
%             tempstore = tempstore';
%         end
%         % 验证时间向量与信号尺寸匹配
%         if sum(size(tempstore) ~= size(signal))
%             error('时间向量必须与信号尺寸相同');
%         end
%         time = tempstore;  % 存储有效时间向量
%         clear tempstore    % 清理临时变量
%     
%     % 检查是否为"midbit"参数 (比特中心标记)
%     elseif strcmpi("midbit", varargin{i})
%         % 验证输入为布尔值或数值型 (0/1)
%         if ~(isa(varargin{i+1},"logical") || isa(varargin{i+1},"double"))
%             error('midbit参数需为布尔值');
%         end
%         midbit = varargin{i+1};  % 存储设置值
%     
%     % 检查是否为"histogram"参数 (直方图显示)
%     elseif strcmpi("histogram", varargin{i})
%         % 验证输入为布尔值或数值型 (0/1)
%         if ~(isa(varargin{i+1},"logical") || isa(varargin{i+1},"double"))
%             error('histogram参数需为布尔值');
%         end
%         plothist = varargin{i+1};  % 存储设置值
%     
%     % 检查是否为"xresolution"参数 (x轴分辨率)
%     elseif strcmpi("xresolution", varargin{i})
%         % 验证输入为数值型
%         if ~(isa(varargin{i+1},"double"))
%             error('x分辨率需为数值');
%         end
%         % 验证输入为标量
%         if length(varargin{i+1}) > 1
%             error('x分辨率需为标量');
%         end
%         % 验证正值
%         if varargin{i+1} < 1
%             error('x分辨率需为正数');
%         end
%         screen_size(1) = varargin{i+1};  % 存储x分辨率
%     
%     % 检查是否为"yresolution"参数 (y轴分辨率)
%     elseif strcmpi("yresolution", varargin{i})
%         % 验证输入为数值型
%         if ~(isa(varargin{i+1},"double"))
%             error('y分辨率需为数值');
%         end
%         % 验证正值
%         if varargin{i+1} < 1
%             error('y分辨率需为正数');
%         end
%         screen_size(2) = varargin{i+1};  % 存储y分辨率
%     end
% end

% 解析可选参数
for i = 3:2:nargin
    key = varargin{i}; value = varargin{i+1};
    if strcmpi("time", key)
        if size(value,1)>size(value,2), value = value'; end
        if ~isequal(size(value),size(signal))
            error('时间向量与信号尺寸不匹配');
        end
        time = value;
    elseif strcmpi("midbit", key)
        midbit = logical(value);  % 中比特标记
    elseif strcmpi("histogram", key)
        plothist = logical(value);% 直方图开关
    elseif strcmpi("xresolution", key)
        screen_size(1) = value;   % x分辨率
    elseif strcmpi("yresolution", key)
        screen_size(2) = value;   % y分辨率
    end
end

%% 分辨率设置
if sum(screen_size)==0
    screen_size = round(500*[1,0.776]); % 默认分辨率500x388
elseif screen_size(2)==0
    screen_size(2) = round(screen_size(1)*0.776); % 按宽高比计算y
elseif screen_size(1)==0
    screen_size(1) = round(screen_size(2)*1.2886); % 按宽高比计算x
end
screen = zeros(screen_size); % 初始化眼图矩阵

%% 信号预处理
% 截断信号使总长度可被分段整除
num2delete = rep*mod(numel(signal)/rep,1);
signal = signal(round(num2delete)+1:end); 
% 分段重组：每行代表一个眼图周期
sig_sliced = reshape(signal, rep, [])'; 

%% 坐标轴设置
% y轴（电压）分箱边界
yscale = linspace(min(signal), max(signal), screen_size(2));
% x轴（时间）刻度（若有时间向量）
if isempty(time)
    xscale = [];
else
    xscale = linspace(0, time(rep)-time(1), screen_size(1));
end
scr_out.x = xscale;  % 存储输出
scr_out.y = yscale;

%% 电压分箱处理
for i = 1:size(sig_sliced,1)   % 遍历每个周期
    slice = sig_sliced(i,:);
    % 将当前周期插值到目标x分辨率
    sliceinterp = interp1(linspace(0,1,rep), slice, linspace(0,1,screen_size(1)));
    for j = 1:numel(sliceinterp)  % 遍历周期内每个点
        [~, bin] = min(abs(sliceinterp(j)-yscale)); % 寻找最近电压分箱
        screen(j, bin) = screen(j, bin) + 1;        % 计数累加
    end
end

%% 图像后处理
screen = flipud(screen');      % 转置并垂直翻转
scr_out.raw_data = screen;     % 存储原始数据
image = screen/max(screen(:)); % 归一化
image = sqrt(image);           % Gamma校正(γ=0.5)
scr_out.image = image/max(image(:)); % 再次归一化

%% 绘图
fh = figure();
hEye = axes(fh);
imshow(image);  % 显示眼图
hold on;

%% 自动标记比特0/1电平（y轴）
try
    % 通过直方图分析寻找逻辑电平
    [counts, edges] = histcounts(signal, floor(screen_size(2)/2));
    bin_centers = edges(1:end-1) + diff(edges)/2;
    
    % 使用二阶导数定位直方图峰值
    [~, idx1] = min(gradient(gradient(lowpass(counts,0.1))));
    level1 = bin_centers(idx1);
    
    % 在未分析区域寻找第二峰值
    half = floor(screen_size(2)/4);
    if idx1 < half
        counts2 = counts(half:end);
        offset = half;
    else
        counts2 = counts(1:half);
        offset = 0;
    end
    [~, idx2] = min(gradient(gradient(lowpass(counts2,0.1))));
    level2 = bin_centers(idx2+offset);
    
    % 按中心距离排序
    [~, order] = sort(abs([idx1, idx2+offset]-screen_size(2)/2));
    levels = [level1, level2];
catch
    warning('y轴标记失败'); 
end

%% 自动标记过零点（x轴）
try
    [~, filter_coeff] = lowpass(sig_sliced(1,:), 0.1);%使用lowpass函数对切片进行低通滤波，以平滑直方图并减小噪声
    filter_radius = ceil(numel(filter_coeff)/2);  % 滤波器半宽
    zero_crossings = [];
    
    % 检测每个周期的过零点
    for i = 1:size(sig_sliced,1)
        if max(sig_sliced(i,:))<=0 || min(sig_sliced(i,:))>=0, continue; end
        
        % 边界扩展滤波
        padded = [repmat(sig_sliced(i,1),1,filter_radius), ...
                 sig_sliced(i,:), ...
                 repmat(sig_sliced(i,end),1,filter_radius)];
        filtered = filter(filter_coeff,1,padded);
        filtered = filtered(2*filter_radius+1:end); % 去除边界效应
        
        % 根据起始极性检测过零点
        if sig_sliced(i,1) > 0
            cross_idx = find(filtered<0,1);
        else
            cross_idx = find(filtered>0,1);
        end
        zero_crossings(end+1) = cross_idx;
    end
    
    % 聚类分析过零点
    zero_crossings = sort(zero_crossings);
    gaps = find(diff(zero_crossings) > 25); % 寻找显著间隔
    clusters = [1, gaps, numel(zero_crossings)];
    mean_crossings = arrayfun(@(k) mean(zero_crossings(clusters(k)+1:clusters(k+1)-1)),...
        1:numel(clusters)-1);
    
    % 转换为图像坐标
    scale_factor = screen_size(1)/rep;
    img_crossings = mean_crossings * scale_factor;
    step_size = mean(diff(img_crossings));
    
    % 中比特调整
    if midbit
        img_crossings = img_crossings + step_size/2;
        img_crossings = [img_crossings(1)-2*step_size, ...
                         img_crossings(1)-step_size, ...
                         img_crossings, ...
                         img_crossings(end)+step_size];
        img_crossings = img_crossings(img_crossings>0 & img_crossings<screen_size(1));
    end
    
    % 设置x轴刻度和标签
    xticks(round(img_crossings));
    if ~isempty(time)
        time_vals = time(round(img_crossings/scale_factor));
        [unit_prefix, scale] = get_prefix(max(time_vals));
        xticklabels(round((time_vals-min(time_vals))*10^-scale,1));
        xlabel(['时间 [', unit_prefix, 's]']);
    else
        xticklabels(0:numel(img_crossings)-1);
        xlabel('比特周期');
    end
catch
    warning('x轴标记失败'); 
end

%% 设置y轴标签
try
    max_voltage = max(abs([min(signal),max(signal)]));
    [unit_prefix, scale] = get_prefix(max_voltage);
    yticks(sort(abs([idx1, idx2+offset]-screen_size(2)/2)*2)); % 按中心距离排序
    yticklabels(round(levels(order)*10^-scale,1));
    ylabel(strcat("Voltage [",unit_prefix,"V]"));
catch
    warning('y轴标记失败'); 
end

%% 直方图绘制
if plothist
    hXHist = axes('Position',[0.1 0.85 0.8 0.1]); % 时间直方图位置
    histogram(zero_crossings, 'BinWidth', 5, ...
             'FaceColor','k', 'EdgeColor','none');
    axis tight; set(hXHist,'Visible','off');
    
    hYHist = axes('Position',[0.85 0.1 0.1 0.8]); % 电压直方图位置
    histogram(signal, 'BinWidth', diff(yscale(1:2)));
    view(90,-90);  % 旋转90度
    axis tight; set(hYHist,'Visible','off');
end

%% 单位前缀辅助函数
function [prefix, power] = get_prefix(value)
    exponents = -12:3:12; % 指数范围
    prefixes = ["p","n","µ","m","","k","M","G","T"]; % 单位前缀
    exp_idx = find(value >= 10.^exponents, 1, 'last'); % 匹配最佳前缀
    if isempty(exp_idx)
        exp_idx=1; 
    end
    prefix = prefixes(exp_idx);
    power = exponents(exp_idx);
end

end
