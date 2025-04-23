function [ RxBit ] = QAM8Demapping( RxSym, BitPerSymbol )

RxBit=[];k= 1;
if BitPerSymbol == 3
    for i=1:length(RxSym)
        switch RxSym(i)
            case -1-1i
                RxBit(k) =0;
                RxBit(k+1)=0;
                RxBit(k+2)=0;
            case -1-0i
                RxBit(k) =0;
                RxBit(k+1)=0;
                RxBit(k+2)=1;
            case 0-1i
                RxBit(k) =0;
                RxBit(k+1)=1;
                RxBit(k+2)=0;
            case -1+1i
                RxBit(k) =0;
                RxBit(k+1)=1;
                RxBit(k+2)=1;
            case 1+0i
                RxBit(k) =1;
                RxBit(k+1)=0;
                RxBit(k+2)=0;
            case 1+1i
                RxBit(k) =1;
                RxBit(k+1)=0;
                RxBit(k+2)=1;
            case 1-1i
                RxBit(k) =1;
                RxBit(k+1)=1;
                RxBit(k+2)=0;
            case  0+1i
                RxBit(k) =1;
                RxBit(k+1)=1;
                RxBit(k+2)=1;

        end
        k=k+3;

    end

end
