function [pam16demod] = pam16demod(RxSym)
   pam16demod=zeros(1,length(RxSym)*4);
   k = 1;
for i=1:length(RxSym)
   
    switch RxSym(i)
        case 0
            pam16demod(k) =0;
            pam16demod(k+1)=0;
            pam16demod(k+2)=0;
             pam16demod(k+3)=0;
        case 1
            pam16demod(k) =0;
            pam16demod(k+1)=0;
            pam16demod(k+2)=0;
            pam16demod(k+3)=1;
        case 2
            pam16demod(k) =0;
            pam16demod(k+1)=0;
            pam16demod(k+2)=1;
            pam16demod(k+3)=1;
        case 3
            pam16demod(k) =0;
            pam16demod(k+1)=0;
            pam16demod(k+2)=1;
            pam16demod(k+3)=0;
        case 4
            pam16demod(k) =0;
            pam16demod(k+1)=1;
            pam16demod(k+2)=1;
            pam16demod(k+3)=0;
        case 5
            pam16demod(k) =0;
            pam16demod(k+1)=1;
            pam16demod(k+2)=1;
            pam16demod(k+3)=1;
          case 6
            pam16demod(k) =0;
            pam16demod(k+1)=1;
            pam16demod(k+2)=0;
            pam16demod(k+3)=1;
          case 7
            pam16demod(k) =0;
            pam16demod(k+1)=1;
            pam16demod(k+2)=0;
            pam16demod(k+3)=0;
         case 8
            pam16demod(k) =1;
            pam16demod(k+1)=1;
            pam16demod(k+2)=0;
            pam16demod(k+3)=0;
        case 9
            pam16demod(k) =1;
            pam16demod(k+1)=1;
            pam16demod(k+2)=0;
            pam16demod(k+3)=1;
        case 10
            pam16demod(k) =1;
            pam16demod(k+1)=1;
            pam16demod(k+2)=1;
            pam16demod(k+3)=1;
        case 11
            pam16demod(k) =1;
            pam16demod(k+1)=1;
            pam16demod(k+2)=1;
            pam16demod(k+3)=0;
        case 12
            pam16demod(k) =1;
            pam16demod(k+1)=0;
            pam16demod(k+2)=1;
            pam16demod(k+3)=0;
        case 13
            pam16demod(k) =1;
            pam16demod(k+1)=0;
            pam16demod(k+2)=1;
            pam16demod(k+3)=1;
          case 14
            pam16demod(k) =1;
            pam16demod(k+1)=0;
            pam16demod(k+2)=0;
            pam16demod(k+3)=1;
          case 15
            pam16demod(k) =1;
            pam16demod(k+1)=0;
            pam16demod(k+2)=0;
            pam16demod(k+3)=0;

    end
    k=k+4;
end

     end
   
      