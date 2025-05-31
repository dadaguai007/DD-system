%% Equalizers 时域MMSE均衡器 ， 得到权重向量
function [w]=getMMSEW(G,E,snr)
% 计算得到时域权重
    I = eye(size(G,1));
    w = (G*G' + I./snr)\(G*E);
end