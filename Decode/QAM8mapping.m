function [ Symbol ] = QAM8mapping( Bit, BitPerSymbol )

% Symbol: output QAM symbol sequence
% Bit: input bit sequence
% BitPerSymbol: modulation format
N = length(Bit);
if BitPerSymbol == 3
      Symbol = zeros(1,floor(N/3));
    for ii = 1:floor(N/3)
        switch 4*Bit(3*ii-2)+2*Bit(3*ii-1)+1*Bit(3*ii)
            case 0
                Symbol(ii) = -1-1i;
            case 1
                Symbol(ii) = -1-0i;
            case 2
                Symbol(ii) = 0-1i;
            case 3
                Symbol(ii) = -1+1i;
            case 4
                 Symbol(ii) = 1+0i;
            case 5
                Symbol(ii) = 1+1i;
            case 6
                Symbol(ii) = 1-1i;
            case 7
                Symbol(ii) = 0+1i;
        end
    end
end
end