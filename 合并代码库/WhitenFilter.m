function FullResultMat=WhitenFilter(OriData,LayerNum)
% LayerNum通道数
for i=1:LayerNum  % 遍历每个通道
    ChosenData = OriData(i,:); 
    
    % 步骤1: 计算自相关函数
    [Corr,lg] = xcorr(ChosenData); 
    Corr(lg<0) = [];  % 保留非负延迟部分
    
    % 步骤2: Levinson-Durbin算法求反射系数
    [A,E,K] = levinson(Corr,N);  % K=反射系数向量
    
    % 步骤3: 格型滤波器白化处理
    [WhiteningResult,~] = latcfilt(K, ChosenData); 
    FullResultMat(i,:) = WhiteningResult; % 存储结果
end
