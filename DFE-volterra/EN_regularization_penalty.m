% L1和L2的正则化准则
function e=EN_regularization_penalty(yd,y,w,lamba,alera)

d=1/2*abs(yd-y).^2;
L=lamba*sum(alera*abs(w)+(1-alera)*abs(w).^2);
e=d+L;

end


