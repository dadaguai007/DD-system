function [pam4demod] = pam4demod(RxSym)
   pam4demod=zeros(1,length(RxSym)*2);
   k = 1;
for i=1:length(RxSym)
   
    switch RxSym(i)
        case 0
            pam4demod(k) =0;
            pam4demod(k+1)=0;
        case 1
             pam4demod(k) = 0;
              pam4demod(k+1) = 1;
        case 2
             pam4demod(k)=1;
              pam4demod(k+1)=1;
        case 3
             pam4demod(k) = 1;
              pam4demod(k+1) = 0;
    end
    k=k+2;
    
     
end

     end
   
      