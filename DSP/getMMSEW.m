%% Equalizers 时域MMSE均衡器 ， 得到权重向量
function [w]=getMMSEW(G,E,snr)
% 计算得到时域权重
% G为信道卷积矩阵，SNR为噪声量，E为单位脉冲
% G*G'等效为接收信号的自相关矩阵；G*E等效为接收信号与发射信号的互相关向量，中间为噪声分量
% 需要噪声为白
    I = eye(size(G,1));
    w = (G*G' + I./snr)\(G*E);
end