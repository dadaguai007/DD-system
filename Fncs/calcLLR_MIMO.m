function LLRs = calcLLR_MIMO(rxSymb, sigma2, K,states,statesIndex, bitMap, px)
% LLR calculation

M = length(constSymb);
b = log2(M);

% 设置为空，不断累加
LLRs = [];
% 状态数
num_candidates = size(states, 1);
% 状态向量
Gamma = zeros(1,num_candidates);
% 此方法为app方法，求解LLR(最精确的方法)


for i = 1:length(rxSymb)
    for j = 1:num_candidates
        states_candidate = states(j,:);
        Gamma(j)=sum(exp(-abs(yEq(i) - states_candidate.').^2 / sigma2) * px(statesIndex(j,:)));
    end
    prob=Gamma;
    % prob = exp(-abs(rxSymb(i) - states).^2 / sigma2) .* px;
    % prob =max( (-abs(rxSymb(i) - constSymb).^2 / sigma2) + log(px));
    
    for index=1:K

        for indBit = b*(index-1)+1:b*index
            p0 = sum(prob(bitMap(:, indBit) == 0));
            p1 = sum(prob(bitMap(:, indBit) == 1));
            % 索引，一个符号对应log2(M),出去后进行转置。
            LLR((i-1)* b+ indBit) = log(p0) - log(p1);
        end
    end
    LLRs=[LLRs,LLR];
end
end
