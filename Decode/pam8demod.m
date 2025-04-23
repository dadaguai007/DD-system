function [pam8demod] = pam8demod(RxSym)
   pam8demod=zeros(1,length(RxSym)*3);
   k = 1;
for i=1:length(RxSym)
   
    switch RxSym(i)
        case 0
            pam8demod(k) =0;
            pam8demod(k+1)=0;
            pam8demod(k+2)=0;
        case 1
            pam8demod(k) =0;
            pam8demod(k+1)=0;
            pam8demod(k+2)=1;
        case 2
            pam8demod(k) =0;
            pam8demod(k+1)=1;
            pam8demod(k+2)=1;
        case 3
            pam8demod(k) =0;
            pam8demod(k+1)=1;
            pam8demod(k+2)=0;
        case 4
            pam8demod(k) =1;
            pam8demod(k+1)=1;
            pam8demod(k+2)=1;
        case 5
            pam8demod(k) =1;
            pam8demod(k+1)=1;
            pam8demod(k+2)=0;
          case 6
            pam8demod(k) =1;
            pam8demod(k+1)=0;
            pam8demod(k+2)=0;
          case 7
            pam8demod(k) =1;
            pam8demod(k+1)=0;
            pam8demod(k+2)=1;

    end
    k=k+3;
    
     
end

     end
   
      