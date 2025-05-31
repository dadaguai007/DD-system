function w=MMSE_TDE(num_taps, channell)
    % 生成G矩阵，基于信道响应和抽头数
    G = getG(channell, num_taps);
    E = getE(length(channell)+num_taps-1);
    disp("G mat is: ")
    disp(G)
    disp("E mat is: ")
    disp(E)
    w = getMMSEW(G',E,1);
    disp("w for 0 SNR is: ")
    disp(w)
end

%% get e vector
function [E] = getE(Esize)
    E = zeros(Esize,1);
    E(1,1)=1;
end