function [ out ] = QAMdecision( in, BitPerSymbol )

% BitPerSymbol: modulation format
[N1,N2] = size(in);
out = zeros(N1,N2);
if BitPerSymbol == 2
    Symbol = [-1+1i, -1-1i, 1+1i, 1-1i];
    else if BitPerSymbol == 4
            Symbol = [-3+3*1i,-3+1i,-3-3*1i,-3-1i,-1+3*1i,-1+1i,-1-3*1i,-1-1i,3+3*1i,3+1i,3-3*1i,3-1i,1+3*1i,1+1i,1-3*1i,1-1i];
        end
end

if BitPerSymbol == 1.5
    Symbol = [-1-1i, -1-0i, 0-1i, -1+1i,1+0i,1+1i,1-1i,0+1i];
end


for ii = 1:N1
    for jj = 1:N2
        Dis = abs(in(ii,jj)-Symbol).^2;
        [~, min_ind]=min(Dis);
        out(ii,jj) = Symbol(min_ind);
    end
end

end

