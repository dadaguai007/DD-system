% TPCA
function w_var=weight_variance(w,Y,lamba)

weight=(w*w.')*lamba;
total_Y=Y*(Y.');
w_var=weight/total_Y;

end

