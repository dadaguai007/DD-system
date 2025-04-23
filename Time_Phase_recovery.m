[xn,P]=Quick_Syn(data,ref2,1/fs,1/fb/sps);
xn(1:P(1))=[];
%     xn = xn(1:length(ref2)*5);
xn=pnorm(xn);