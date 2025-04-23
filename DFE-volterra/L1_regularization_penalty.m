% L1正则化代价函数
function e=L1_regularization_penalty(yd,y,w,lamba)

e=1/2*abs(yd-y).^2+lamba*sum(abs(w));


end
