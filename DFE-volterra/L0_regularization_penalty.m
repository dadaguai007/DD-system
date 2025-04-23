% L0正则化代价函数
function e=L0_regularization_penalty(yd,y,w,lamba)
loc=find(abs(w)>0);
L0=length(loc);
e=1/2*abs(yd-y).^2+lamba*sum(L0);


end
