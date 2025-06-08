%% get e vector
% 用于噪声消除以及MMSE滤波中
function [E] = getE(Esize)
    E = zeros(Esize,1);
    E(1,1)=1;
end