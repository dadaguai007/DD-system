function [ Symbol ] = QAMmapping( Bit, BitPerSymbol )

% Symbol: output QAM symbol sequence
% Bit: input bit sequence
% BitPerSymbol: modulation format
N = length(Bit);
if BitPerSymbol == 1.5
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


if BitPerSymbol == 2
    Symbol = zeros(1,floor(N/2));
    for ii = 1:floor(N/2)
        switch 2*Bit(2*ii-1)+1*Bit(2*ii)
            case 0
                Symbol(ii) = -1+1i;
            case 1
                Symbol(ii) = -1-1i;
            case 2
                Symbol(ii) = 1+1i;
            case 3
                Symbol(ii) = 1-1i;
        end
    end
    
else if BitPerSymbol == 4
        Symbol = zeros(1,floor(N/4));
        for ii = 1:floor(N/4)
            switch 8*Bit(4*ii-3)+4*Bit(4*ii-2)+2*Bit(4*ii-1)+1*Bit(4*ii)
                case 0
                    Symbol(ii) = -3+3*1i;
                case 1
                    Symbol(ii) = -3+1i;
                case 2
                    Symbol(ii) = -3-3*1i;
                case 3
                    Symbol(ii) = -3-1i;
                case 4
                    Symbol(ii) = -1+3*1i;
                case 5
                    Symbol(ii) = -1+1i;
                case 6
                    Symbol(ii) = -1-3*1i;
                case 7
                    Symbol(ii) = -1-1i;
                case 8
                    Symbol(ii) = 3+3*1i;
                case 9
                    Symbol(ii) = 3+1i;
                case 10
                    Symbol(ii) = 3-3*1i;
                case 11
                    Symbol(ii) = 3-1i;
                case 12
                    Symbol(ii) = 1+3*1i;
                case 13
                    Symbol(ii) = 1+1i;
                case 14
                    Symbol(ii) = 1-3*1i;
                case 15
                    Symbol(ii) = 1-1i;
            end
        end
    end
end

end

